#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 (AGPL-3.0).
# 
# The full copyright notice, including terms governing use, modification,
# and redistribution, is contained in the files LICENSE and COPYRIGHT,
# which can be found at the root of the source code distribution tree:
# - LICENSE <https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE>
# - COPYRIGHT <https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT>
#
"""
Contains the Output class definition

Please note that this module is private. The Output class is
available in the main ``galaxia_ananke`` namespace - use that instead.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Union, Tuple, List, Dict, Iterable
from numpy.typing import NDArray, ArrayLike
from warnings import warn
from functools import cached_property
import concurrent.futures
import pathos
import gc
import os
import pathlib
import itertools
import numpy as np
import h5py as h5
import ebf
import vaex
import pandas as pd
from astropy import units, coordinates
import vaex.dataframe

from .__metadata__ import *
from ._constants import *
from ._templates import *
from ._defaults import *
from .utils import classproperty, CallableDFtoNone, CallableDFtoInt, RecordingDataFrame, State, common_entries, metadata_dicts_to_consolidated_json
from .photometry.PhotoSystem import PhotoSystem
from . import Input

if TYPE_CHECKING:
    from . import Survey

__all__ = ['Output']


def shift_g_lon(lon): # restrict longitude values to be within (-180,180)
    return -((-lon+180)%360-180)


def _flush_extra_columns_to_hdf5(vaex_df: vaex.DataFrame, hdf5_file: pathlib.Path, with_columns: Optional[Iterable] = (), verbose: bool = True) -> None:  # temporary until vaex supports it
    _temp = vaex.open(hdf5_file)
    old_column_names = set(_temp.column_names)
    _temp.close()
    extra_columns = [k for k in set(vaex_df.column_names)-old_column_names if not k.startswith('__')]
    with_columns = list(set(with_columns) - set(extra_columns))
    with h5.File(hdf5_file, 'r+') as f5:
        for k in extra_columns:
            f5.create_dataset(name=k, data=vaex_df[k].to_numpy())
        if verbose and extra_columns:
            print(f"Exported the following quantities to {hdf5_file}")
            print(extra_columns)
        for k in with_columns:
            f5[k][...] = vaex_df[k].to_numpy()
        if verbose and len(with_columns):
            print(f"Overwritten the following quantities to {hdf5_file}")
            print(with_columns)


def _decorate_post_processing(pp: CallableDFtoNone, hdf5_path_input: bool = False, flush_with_columns: Optional[Iterable] = (), max_thread_workers: int = None, verbose: bool = True) -> CallableDFtoNone:
    def new_pp(*args) -> None:
        first_arg = args[0]
        if not isinstance(first_arg, list):
            first_arg = [first_arg]
        for _temp in first_arg:
            if hdf5_path_input:
                hdf5_file: pathlib.Path = _temp
                old_vaex_main_executor = vaex.dataframe.main_executor
                vaex.dataframe.main_executor = vaex.execution.ExecutorLocal(vaex.multithreading.ThreadPoolIndex(max_workers=max_thread_workers))
                vaex_df: vaex.DataFrame = vaex.open(hdf5_file)
            else:
                vaex_df: vaex.DataFrame = _temp
            pp(vaex_df, *args[1:])
            if hdf5_path_input:
                _flush_extra_columns_to_hdf5(vaex_df, hdf5_file, with_columns=flush_with_columns, verbose=verbose)
                vaex_df.close()
                vaex.dataframe.main_executor = old_vaex_main_executor
            gc.collect()
    return new_pp


class Output:
    _position_prop = (('px', 'py', 'pz'), "Position coordinates in $kpc$")
    _velocity_prop = (('vx', 'vy', 'vz'), "Velocity coordinates in $km/s$")
    _celestial_prop = (('ra', 'dec'), "Celestial equatorial coordinates in $degrees$")
    _galactic_prop = (('glon', 'glat'), "Celestial galactic coordinates in $degrees$")
    _distance_prop = ('rad', "Distance in $kpc$")
    _modulus_prop = ('dmod', "Distance modulus in magnitude units")
    _trgbmass_prop = ('mtip', "Tip of the Red Giant Branch stellar mass in solar masses")
    _currentmass_prop = ('mact', "Current stellar mass in solar masses")
    _zamsmass_prop = ('smass', "Zero Age Main Sequence stellar mass in solar masses")
    _age_prop = ('age', "Stellar ages in years and decimal logarithmic scale")
    _surfacegravity_prop = ('grav', "Surface gravity in CGS units and decimal logarithmic scale")
    _metallicity_prop = ('feh', "Stellar metallicity $[Fe/H]$ in $dex$ relative to solar")
    _temperature_prop = ('teff', "Surface temperature in Kelvin and decimal logarithmic scale")
    _luminosity_prop = ('lum', "Stellar luminosity in solar luminosities and decimal logarithmic scale")
    _parentindex_prop = Input._parentindex_prop
    _partitionindex_prop = Input._partitionindex_prop
    _satindex_prop = Input._populationindex_prop
    _satindex = 'satid'
    _particleflag_prop = ('partid', "Flag = 1 if star not at center of its parent particle")
    _parallax_prop = ('pi', "Parallax in milliarcseconds")
    _propermotion_prop = (('mura', 'mudec'), "Equatorial proper motions in milliarcseconds per year")
    _galacticpropermotion_prop = (('mul', 'mub'), "Galactic proper motions in milliarcseconds per year")
    _radialvelocity_prop = ('vr', "Radial velocity in $km/s$")
    _vaex_under_list = ['_repr_html_']
    def __init__(self, survey: Survey) -> None:
        """
            Driver to exploit the output of Galaxia.
            
            Call signature::
                output = Output(survey)
            
            Parameters
            ----------
            survey : :obj:`Survey`
                Survey object that returned this output.
            
            Notes
            -----
            An Output object almost behaves as a vaex DataFrame, also please
            consult ``vaex`` online tutorials for more hands-on information:
            
                https://vaex.io/docs/tutorial.html
            
            The DataFrame represents the catalogue with columns corresponding
            to properties of the stars from the synthetic stellar population
            it simulates.
            
            .. warning:: When generated directly by ``galaxia_ananke``, the
                         catalogue properties reflect directly the quantities
                         as computed by Galaxia. However the catalogue can be
                         modified/amended by applying post-processing routines
                         using the method
                         ``apply_post_process_pipeline_and_flush``. Also if
                         such ``Output`` object was generated by other software
                         than ``galaxia_ananke``, post-processing may have been
                         applied: also please refer to that software
                         documentation for a more complete overview of the
                         catalogue.
            
            The catalogue properties include the photometric magnitudes per
            filter, with each filter identified by a key in the following
            lowercase format:
            
                ``photosys_filtername``
            
            where

            * ``photo_sys`` corresponds to the chosen photometric system
            * ``filtername`` corresponds to a filter name of that system
            
            As an example, the photometry in filters ``gbp``, ``grp`` & ``g``
            of the Gaia DR2 system identified as ``GAIA__DR2`` are respectively
            under keys ``gaia__dr2_gbp``, ``gaia__dr2_grp`` & ``gaia__dr2_g``.

            With those are also always included the following properties:
            {_output_properties}
            
            Additionally, depending on what optional properties were provided
            with the input particle data, the output can also include the
            following properties:
            {_optional_properties}
        """
        self.__survey: Survey = survey
        self.__vaex = None
        self.__vaex_per_partition = None
        self.__path = None
        self.__make_state()
        if not self.caching:
            self.__clear_ebfs(force=True)
        self._max_pp_workers = 1
        self._pp_auto_flush = True

    class _State(State):
        pass

    @classproperty
    def _export_properties(cls):
        return {
            cls._position_prop,
            cls._velocity_prop,
            cls._distance_prop,
            cls._modulus_prop,
            cls._trgbmass_prop,
            cls._currentmass_prop,
            cls._zamsmass_prop,
            cls._age_prop,
            cls._surfacegravity_prop,
            cls._metallicity_prop,
            cls._temperature_prop,
            cls._luminosity_prop,
            cls._parentindex_prop,
            cls._particleflag_prop,
            cls._partitionindex_prop
            }
    
    @classproperty
    def _postprocess_properties(cls):
        return {
            cls._celestial_prop,
            cls._galactic_prop,
            cls._parallax_prop,
            cls._propermotion_prop,
            cls._galacticpropermotion_prop,
            cls._radialvelocity_prop
            }
    
    @classproperty
    def _all_optional_properties(cls):
        return Input._optional_properties \
              - {cls._parentindex_prop, cls._partitionindex_prop, cls._satindex_prop} \
                | {(cls._satindex, cls._satindex_prop[1])}
    
    @classproperty
    def _export_keys(cls):
        return tuple(sum([(_p[0],) if isinstance(_p[0], str) else _p[0] for _p in cls._export_properties], ()))
    
    @classproperty
    def _postprocess_keys(cls):
        return tuple(sum([(_p[0],) if isinstance(_p[0], str) else _p[0] for _p in cls._postprocess_properties], ()))
    
    @classproperty
    def _pos(cls):
        return cls._position_prop[0]
    
    @classproperty
    def _vel(cls):
        return cls._velocity_prop[0]
    
    @classproperty
    def _cel(cls):
        return cls._celestial_prop[0]
    
    @classproperty
    def _gal(cls):
        return cls._galactic_prop[0]
    
    @classproperty
    def _rad(cls):
        return cls._distance_prop[0]
    
    @classproperty
    def _dmod(cls):
        return cls._modulus_prop[0]
    
    @classproperty
    def _mtip(cls):
        return cls._trgbmass_prop[0]
    
    @classproperty
    def _mact(cls):
        return cls._currentmass_prop[0]
    
    @classproperty
    def _mini(cls):
        return cls._zamsmass_prop[0]
    
    @classproperty
    def _age(cls):
        return cls._age_prop[0]
    
    @classproperty
    def _grav(cls):
        return cls._surfacegravity_prop[0]
    
    @classproperty
    def _feh(cls):
        return cls._metallicity_prop[0]
    
    @classproperty
    def _teff(cls):
        return cls._temperature_prop[0]
    
    @classproperty
    def _lum(cls):
        return cls._luminosity_prop[0]
    
    @classproperty
    def _parentid(cls):
        return cls._parentindex_prop[0]
    
    @classproperty
    def _partitionid(cls):
        return cls._partitionindex_prop[0]
        
    @classproperty
    def _partid(cls):
        return cls._particleflag_prop[0]
    
    @classproperty
    def _pi(cls):
        return cls._parallax_prop[0]
    
    @classproperty
    def _mu(cls):
        return cls._propermotion_prop[0]
    
    @classproperty
    def _mugal(cls):
        return cls._galacticpropermotion_prop[0]
    
    @classproperty
    def _vr(cls):
        return cls._radialvelocity_prop[0]
    
    def __dir__(self):
        return sorted({i for i in self.__vaex.__dir__() if not i.startswith('_')}.union(
            super(Output, self).__dir__()).union(
            self._vaex_under_list if self.__vaex is not None else []))

    def __repr__(self) -> str:
        return repr(self._vaex)
    
    def __getitem__(self, item: str):
        return self._vaex[item]
    
    def __setitem__(self, item: str, value):
        self._vaex[item] = value
    
    def __getattr__(self, item: str):
        if (item in self.__vaex.__dir__() and not item.startswith('_')) or (item in self._vaex_under_list and self.__vaex is not None):
            return getattr(self.__vaex, item)
        else:
            return self.__getattribute__(item)

    @classmethod
    def _compile_export_mag_names(cls, photosystems: list[PhotoSystem]) -> Tuple[str]:
        return tuple(itertools.chain.from_iterable([photosystem.to_export_keys for photosystem in photosystems]))
    
    @classmethod
    def _make_export_keys(cls, photosystems: list[PhotoSystem], extra_keys=()) -> Tuple[str]:
        return tuple(set(cls._export_keys).union(extra_keys).union(cls._compile_export_mag_names(photosystems)))

    @classmethod
    def _make_catalogue_keys(cls, photosystems: list[PhotoSystem], extra_keys=()) -> Tuple[str]:
        return cls._make_export_keys(photosystems, extra_keys=cls._postprocess_keys+extra_keys)

    def _make_input_optional_keys(self) -> Tuple[str]:
        return tuple(k if k != self._satindex_prop[0] else self._satindex for k in self.survey.input.optional_keys())

    def _redefine_partitions_in_ebfs(self, partitioning_rule: Optional[CallableDFtoInt]) -> None:
        if partitioning_rule is not None:  # TODO test instead partitioning_rule?
            ebfs: List[pathlib.Path] = self._ebfs
            export_keys: Tuple[str] = self.export_keys
            # with concurrent.futures.ThreadPoolExecutor(max_workers=self._max_pp_workers) as executor:  # credit to https://www.squash.io/how-to-parallelize-a-simple-python-loop/
            #     # Submit tasks to the executor
            #     futures = [executor.submit(self._singlethread_redefine_partitions, ebf_file, partitioning_rule, export_keys)
            #                for ebf_file in ebfs]
            #     # Collect the results
            #     _ = [future.result() for future in concurrent.futures.as_completed(futures)]
            with pathos.pools.ProcessPool(self._max_pp_workers) as executor:  # credit to https://github.com/uqfoundation/pathos/issues/158#issuecomment-449636971
                # Submit tasks to the executor
                futures = [executor.apipe(self._singlethread_redefine_partitions, ebf_file, partitioning_rule, export_keys)
                        for ebf_file in ebfs]
                # Collect the results
                _ = [future.get() for future in futures]

    @classmethod
    def _singlethread_redefine_partitions(cls, ebf_file: pathlib.Path, partitioning_rule: CallableDFtoInt, export_keys: Tuple[str]) -> None:
        dummy_df = RecordingDataFrame([], columns=export_keys, dtype=float)
        _ = partitioning_rule(dummy_df)
        ebf_df = pd.DataFrame({key: ebf.read(str(ebf_file), f"/{key}") for key in dummy_df.record_of_all_used_keys})
        new_partition_id = partitioning_rule(ebf_df)
        del ebf_df
        ebf.update_ind(str(ebf_file), f'/{cls._partitionid}', new_partition_id)
        partition_id_sorter = np.argsort(new_partition_id)
        del new_partition_id
        gc.collect()
        for key in ebf._EbfMap.keys(str(ebf_file)):
            if key not in [b'/log', b'/center'] and not key.startswith(b'/.'):
                ebf.update_ind(str(ebf_file), key, ebf.read(str(ebf_file), key)[partition_id_sorter])
        # Update ebf file log to flag that a modification was made
        full_log = str(ebf.read(str(ebf_file), '/log')).rstrip('\n')
        ebf._EbfTable.remove(str(ebf_file), '/log')  # TODO update when a more dedicated routine is available
        ebf.write(
            str(ebf_file), '/log',
            f"{full_log}\n"
            f"# Modified by Python's {NAME} v{__version__},\n"
            f"# software available at <{__url__}>.\n",
            'a'
        )
        
    def _ebf_to_hdf5(self) -> None:
        ebfs: List[pathlib.Path] = self._ebfs
        export_keys: Tuple[str] = self.export_keys
        partition_ids: List[int] = list(self._hdf5s.keys())
        # with concurrent.futures.ThreadPoolExecutor(max_workers=self._max_pp_workers) as executor:  # credit to https://www.squash.io/how-to-parallelize-a-simple-python-loop/
        #     # Submit tasks to the executor
        #     futures = [executor.submit(self._singlethread_ebf_to_hdf5, i, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs, ebfs, export_keys)
        #                for i, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs in common_entries(self._hdf5s, self.__ebfs_part_slices, self.__ebfs_part_lengths)]
        #     # Collect the results
        #     _ = [future.result() for future in concurrent.futures.as_completed(futures)]
        with pathos.pools.ProcessPool(self._max_pp_workers) as executor:  # credit to https://github.com/uqfoundation/pathos/issues/158#issuecomment-449636971
            # Submit tasks to the executor
            futures = [executor.apipe(self._singlethread_ebf_to_hdf5, i, partition_id, partition_ids, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs, ebfs, export_keys, self.parameters, self.verbose)
                       for i, (partition_id, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs) in enumerate(common_entries(self._hdf5s, self.__ebfs_part_slices, self.__ebfs_part_lengths))]
            # Collect the results
            _ = [future.get() for future in futures]
        if not(self._pp_auto_flush):
            self.__reload_vaex()

    @classmethod
    def _singlethread_ebf_to_hdf5(cls, i: int, partition_id: int,
                                  partition_ids: List[int],
                                  hdf5_file: pathlib.Path,
                                  part_slices_in_ebfs: Dict[str, List[slice]],
                                  part_lengths_in_ebfs: Dict[str, int],
                                  ebfs: List[pathlib.Path], export_keys: Tuple[str],
                                  parameters: Dict[str, Union[str,float,int]],
                                  verbose: bool = True) -> None:
        header_ebf: str          = str(ebfs[0].resolve())
        ebfs: List[pathlib.Path] = [ebf_path for ebf_path in ebfs if ebf_path.name in part_lengths_in_ebfs]
        n_ebfs: int              = len(ebfs)
        if ebfs:
            i_ebf: int               = i % n_ebfs
            ebfs: List[pathlib.Path] = ebfs[i_ebf:]+ebfs[:i_ebf]
            header_ebf: str          = str(ebfs[0].resolve())
        #
        data_length: int         = sum(part_lengths_in_ebfs.values())
        ebfs_slices: Dict[str, slice] = {
            ebf_path.name: slice(bounds[0],bounds[1])
            for ebf_path, bounds
            in zip(ebfs,
                   np.repeat(
                       np.cumsum(
                           [0]+[part_lengths_in_ebfs[ebf_path.name] for ebf_path in ebfs]
                           ),
                       [1]+(n_ebfs-1)*[2]+[1]
                       ).reshape((n_ebfs,2)
                       ) if ebfs else []
                   )
            }
        with h5.File(hdf5_file, 'w') as f5:
            f5datasets = {name: f5.create_dataset(name=name,
                                                  shape=(data_length,),
                                                  dtype=ebf.getHeader(header_ebf, f"/{name}").get_dtype())
                            for name in export_keys}
            ebf_logs = []
            for ebf_path in ebfs:
                ebf_name: str            = ebf_path.name
                ebf_str: str             = str(ebf_path.resolve())
                ebf_log: Dict[str, str]  = cls.__parse_galaxia_ebf_log(ebf_str)
                f5data_slice: slice      = ebfs_slices[ebf_name]
                part_slices: List[slice] = part_slices_in_ebfs[ebf_name]
                for name in export_keys:
                    head = f5data_slice.start
                    for p_slice in part_slices:
                        f5datasets[name][head:(head:=head+p_slice.stop-p_slice.start)] = ebf.read(
                            ebf_str, f"/{name}", begin=p_slice.start, end=p_slice.stop
                            )
                if verbose:
                    print(f"Exported the following quantities from {ebf_path} to {hdf5_file} for partition {partition_id}")
                    print(list(f5.keys()))
                ebf_logs.append(ebf_log)
            f5.attrs["Galaxia_ebf_output_logs"] = metadata_dicts_to_consolidated_json(ebf_logs, metadata_name = "Galaxia_ebf_output_logs")
            f5.attrs[f"{NAME}_header"] = (
                f"File generated by Python module {NAME} v{__version__}, "
                f"software available at <{__url__}> under {__license__}."
            )
            f5.attrs[f"cat_partition_id"] = partition_id
            f5.attrs[f"all_partition_ids"] = partition_ids
            for k, v in parameters.items():
                f5.attrs[f"{PARAMETER_ATTRIBUTE_PREFIX}{k}"] = v if isinstance(v, (bool, int, float, str)) else str(v)  # TODO consider Iterable case?

    @classmethod
    def __parse_galaxia_ebf_log(cls, ebf_str: str) -> Dict[str, str]:
        log_string: str = str(ebf.read(ebf_str, '/log'))
        log_header, log_body = log_string.split('\n# <parameterfile>\n')
        log_body, log_footer = log_body.split('\n# </parameterfile>\n')
        log_body: str = log_body.replace('\n# ','\n')
        log_dict: Dict[str, str] = {parts[0]: ' '.join(parts[1:]) for parts in (line.split() for line in log_body.splitlines()) if len(parts)>0}
        return {"log_header": log_header, **log_dict, "log_footer": log_footer}

    def read_galaxia_output(self, partitioning_rule: Optional[CallableDFtoInt], max_pp_workers: int, pp_auto_flush: bool) -> None:
        self._max_pp_workers = max_pp_workers
        self._pp_auto_flush = pp_auto_flush
        self.check_state_before_running(description="redefine_partitions_in_ebfs")(self._redefine_partitions_in_ebfs)(partitioning_rule)
        self.check_state_before_running(description="convert_ebf_to_hdf5", level=1)(self._ebf_to_hdf5)()

    ### DEFINING POST PROCESSING PIPELINES BELOW # TODO consider a PostProcess class that runs postprocess pipeline at __call__ and holds flush_with_columns

    @classmethod
    def __pp_convert_cartesian_to_galactic(cls, df: pd.DataFrame) -> None:
        """
        converts positions & velocities from mock catalog Cartesian coordinates (relative to solar position) 
        into Galactic coordinates, assuming Sun is on -x axis (use rotateStars)
        """
        GC = coordinates.Galactic(u = df[cls._pos[0]].to_numpy()*units.kpc,
                                  v = df[cls._pos[1]].to_numpy()*units.kpc,
                                  w = df[cls._pos[2]].to_numpy()*units.kpc,
                                  U = df[cls._vel[0]].to_numpy()*units.km/units.s,
                                  V = df[cls._vel[1]].to_numpy()*units.km/units.s,
                                  W = df[cls._vel[2]].to_numpy()*units.km/units.s,
                                  representation_type = coordinates.CartesianRepresentation,
                                  differential_type   = coordinates.CartesianDifferential)
        df[cls._gal[0]] = shift_g_lon(GC.spherical.lon.value)
        df[cls._gal[1]] = GC.spherical.lat.value
        df[cls._rad]    = GC.spherical.distance.value
        ####################################
        if GC.shape[0]:
            df[cls._mugal[0]] = GC.sphericalcoslat.differentials['s'].d_lon_coslat.value
            df[cls._mugal[1]] = GC.sphericalcoslat.differentials['s'].d_lat.value
            df[cls._vr]       = GC.sphericalcoslat.differentials['s'].d_distance.value
        else:
            df[cls._mugal[0]] = df[cls._mugal[1]] = df[cls._vr] = df[cls._pos[0]].to_numpy()+df[cls._vel[0]].to_numpy()

    @classmethod
    def __pp_convert_galactic_to_icrs(cls, df: pd.DataFrame) -> None:
        """
        converts PMs in galactic coordinates (mulcosb, mub) in arcsec/yr (as output by Galaxia)
        to ra/dec in mas/yr (units of output catalog)
        """
        GC = coordinates.Galactic(l               = df[cls._gal[0]].to_numpy()*units.degree,
                                  b               = df[cls._gal[1]].to_numpy()*units.degree,
                                  distance        = df[cls._rad].to_numpy()*units.kpc,
                                  pm_l_cosb       = df[cls._mugal[0]].to_numpy()*units.mas/units.yr,
                                  pm_b            = df[cls._mugal[1]].to_numpy()*units.mas/units.yr,
                                  radial_velocity = df[cls._vr].to_numpy()*units.km/units.s)
        GC_icrs = GC.transform_to(coordinates.ICRS())
        df[cls._cel[0]] = GC_icrs.ra.value
        df[cls._cel[1]] = GC_icrs.dec.value
        ####################################
        if GC.shape[0]:
            df[cls._mu[0]]  = GC_icrs.pm_ra_cosdec.to(units.mas/units.yr).value
            df[cls._mu[1]]  = GC_icrs.pm_dec.to(units.mas/units.yr).value
        else:
            df[cls._mu[0]] = df[cls._mu[1]] = df[cls._gal[0]].to_numpy()+df[cls._mugal[0]].to_numpy()
    
    @classmethod
    def __pp_convert_icrs_to_galactic(cls, df: pd.DataFrame) -> None:
        """
        converts PMs from ICRS coordinates (muacosd, mudec) to Galactic (mul, mub)
        input and output in mas/yr for PMs and degrees for positions
        also exports the galactic lat and longitude
        """
        IC = coordinates.ICRS(ra           = df[cls._cel[0]].to_numpy()*units.degree,
                              dec          = df[cls._cel[1]].to_numpy()*units.degree,
                              pm_ra_cosdec = df[cls._mu[0]].to_numpy()*units.mas/units.yr,
                              pm_dec       = df[cls._mu[1]].to_numpy()*units.mas/units.yr)

        IC_gal = IC.transform_to(coordinates.Galactic())
        df[cls._gal[0]]   = shift_g_lon(IC_gal.l.value)
        df[cls._gal[1]]   = IC_gal.b.value
        ####################################
        if IC.shape[0]:
            df[cls._mugal[0]] = IC_gal.pm_l_cosb.to(units.mas/units.yr).value
            df[cls._mugal[1]] = IC_gal.pm_b.to(units.mas/units.yr).value
        else:
            df[cls._mugal[0]] = df[cls._mugal[1]] = df[cls._cel[0]].to_numpy()+df[cls._mu[0]].to_numpy()

    @classmethod
    def __pp_last_conversions(cls, df: pd.DataFrame) -> None:
        df[cls._pi]   = 1.0/df[cls._rad]  # parallax in mas (from distance in kpc)
        df[cls._teff] = 10**df[cls._teff]  #Galaxia returns log10(teff/K)
        df[cls._lum]  = 10**df[cls._lum]  #Galaxia returns log10(lum/lsun)

    @property
    def _max_pp_workers(self) -> int:
        return self.__max_pp_workers
    
    @_max_pp_workers.setter
    def _max_pp_workers(self, value: int) -> None:
        self.__max_pp_workers: int = value

    @property
    def _pp_auto_flush(self) -> bool:
        return self.__pp_auto_flush
    
    @_pp_auto_flush.setter
    def _pp_auto_flush(self, value: bool) -> None:
        self.__pp_auto_flush: bool = value

    @staticmethod
    def __consolidate_partitions_per_process(partitions: NDArray, lengths: NDArray, max_workers: int) -> Dict[int, List[int]]:
        df = pd.DataFrame({'partition': partitions, 'length': lengths, 'process_id': 0*lengths})
        df.set_index('partition', drop=True, inplace=True)
        df.sort_values('length', inplace=True)
        df['length_cumsum'] = df.length.cumsum()
        df['length_cumsum_norm'] = (
            df.length_cumsum/df.length_cumsum.iloc[-1]
            if df.length_cumsum.iloc[-1]
            else np.linspace(0,1,df.shape[0]+1)[1:]
            )
        df['process_id'] = np.ceil(df.length_cumsum_norm*max_workers).astype('int')-1
        unique = df.process_id.unique()
        df['process_id'] = df.process_id.map(dict(zip(unique, range(len(unique)))))
        df['temp'] = df.length * (-1)**(df.process_id)
        df.sort_values(['process_id','temp'], inplace=True)
        return df.groupby('process_id').groups

    def apply_post_process_pipeline_and_flush(self, post_process: CallableDFtoNone, *args, flush_with_columns=(), hold_flush: bool = False, hold_reload: bool = False, consolidate_partitions_per_process: bool = False) -> None:
        """
            Apply a given post processing routine to the catalogue

            Parameters
            ----------
            post_process : callable
                Post processing pipeline to apply to the catalogue. This must
                be defined as a callable that returns nothing, and take only
                positional arguments, the first of which being the DataFrame
                representing the catalogue.

            \*args : callable args
                Any other positinoal arguments that should be passed to the
                ``post_process`` callable pipeline, in the order they should
                be passed.

            flush_with_columns : iterable
                If given an iterable structure of existing column keys, the
                flushing done after application of the post-processing
                will also overwrite those in the backend file with their
                current in-memory values. Default to an empty tuple. 

            hold_flush : bool
                Flag to hold the flushing from being done after application of
                the post-processing. Default to False.

            hold_reload : bool
                Flag to hold the reload from being done after application of
                the post-processing and flushing. Default to False.
        """
        # post_process(self._vaex, *args)
        vaex_df_or_hdf5_or_list_s = self._hdf5s if self._pp_auto_flush else self._vaex_per_partition
        if consolidate_partitions_per_process:
            partitions_per_process: Dict[int, List[int]] = self.__consolidate_partitions_per_process(
                np.array(list(self.__partitions_lengths.keys())),
                np.array(list(self.__partitions_lengths.values())),
                self._max_pp_workers
                )
            vaex_df_or_hdf5_or_list_s = {process: [vaex_df_or_hdf5_or_list_s[i] for i in partitions]
                                         for process, partitions in partitions_per_process.items()}
        if self._pp_auto_flush:
            with pathos.pools.ProcessPool(self._max_pp_workers) as executor:  # credit to https://github.com/uqfoundation/pathos/issues/158#issuecomment-449636971
                # Submit tasks to the executor
                futures = [executor.apipe(_decorate_post_processing(post_process,
                                                                    self._pp_auto_flush,
                                                                    flush_with_columns=flush_with_columns,
                                                                    max_thread_workers=int(np.ceil(os.cpu_count()/self._max_pp_workers)),
                                                                    verbose=self.verbose),
                                        vaex_df_or_hdf5_or_list, *args)
                        for vaex_df_or_hdf5_or_list in vaex_df_or_hdf5_or_list_s.values()]
                # Collect the results
                _ = [future.get() for future in futures]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self._max_pp_workers) as executor:  # credit to https://www.squash.io/how-to-parallelize-a-simple-python-loop/
                # Submit tasks to the executor
                futures = [executor.submit(_decorate_post_processing(post_process,
                                                                    self._pp_auto_flush,
                                                                    flush_with_columns=flush_with_columns,
                                                                    verbose=self.verbose),
                                        vaex_df_or_hdf5, *args)
                        for vaex_df_or_hdf5 in vaex_df_or_hdf5_or_list_s]
                # Collect the results
                _ = [future.result() for future in concurrent.futures.as_completed(futures)]
        if not(hold_flush) and not(self._pp_auto_flush):
            self.flush_extra_columns_to_hdf5(with_columns=flush_with_columns)
        if not(hold_reload):
            self.__reload_vaex()

    def post_process_output(self) -> None:
        self.check_state_before_running(description="pp_cartesian_to_galactic")(self._pp_convert_cartesian_to_galactic)(hold_reload=True)
        self.check_state_before_running(description="pp_galactic_to_icrs", level=1)(self._pp_convert_galactic_to_icrs)(hold_reload=True)
        self.check_state_before_running(description="pp_last_conversions", level=1)(self._pp_last_conversions)(hold_reload=True)
        self.__reload_vaex()

    def _pp_convert_cartesian_to_galactic(self, **kwargs) -> None:
        pipeline_name = "convert_cartesian_to_galactic"
        if self.verbose:
            print(f"Running {pipeline_name} post-processing pipeline")
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_cartesian_to_galactic, flush_with_columns=self._gal+(self._rad,), **kwargs)

    def _pp_convert_galactic_to_icrs(self, **kwargs) -> None:
        pipeline_name = "convert_galactic_to_icrs"
        if self.verbose:
            print(f"Running {pipeline_name} post-processing pipeline")
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_galactic_to_icrs, flush_with_columns=self._cel, **kwargs)
    
    def _pp_convert_icrs_to_galactic(self, **kwargs) -> None:
        pipeline_name = "convert_icrs_to_galactic"
        if self.verbose:
            print(f"Running {pipeline_name} post-processing pipeline")
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_icrs_to_galactic, flush_with_columns=self._gal+self._mugal, **kwargs)

    def _pp_last_conversions(self, **kwargs) -> None:
        pipeline_name = "last_conversions"
        if self.verbose:
            print(f"Running {pipeline_name} post-processing pipeline")
        self.apply_post_process_pipeline_and_flush(self.__pp_last_conversions, flush_with_columns=(self._teff, self._lum), **kwargs)

    def __name_with_ext(self, ext):
        name_base = self._file_base
        return name_base.parent / f"{name_base.name}{ext}"
    
    def save(self, path):  # TODO Gotta update this
        """
            Save output to new path

            .. danger:: currently not implemented
        """
        raise NotImplementedError
        old_path = self._path
        self.__path = pathlib.Path(path)
        self._vaex.close()
        old_path.rename(self._path)
        self.__vaex = vaex.open(self._path)

    def __make_state(self):
        self.__state: Output._State = self._State(self.__name_with_ext('.dummy'), self)

    @property
    def _state(self) -> Output._State:
        return self.__state

    @property
    def check_state_before_running(self):
        return self._state.check_state_file_before_running

    @property
    def caching(self):
        return self.survey.caching

    @property
    def verbose(self):
        return self.survey.verbose

    @property
    def survey(self):
        return self.__survey
    
    @property
    def photosystems(self):
        return self.survey.photosystems

    @property
    def export_keys(self) -> Tuple[str]:
        return self._make_export_keys(self.photosystems, extra_keys=self._make_input_optional_keys())
    
    @property
    def catalogue_keys(self) -> Tuple[str]:
        return self._make_catalogue_keys(self.photosystems, extra_keys=self._make_input_optional_keys())
    
    @property
    def output_dir(self):
        return pathlib.Path(self.parameters[FTTAGS.output_dir])

    @property
    def output_name(self):
        return f"{self.survey.surveyname_hash}.{self.survey.inputname_hash}"

    @property
    def rsun_skycoord(self):
        _temp = [self.parameters[k] for k in FTTAGS.rSun]
        return coordinates.SkyCoord(u=_temp[0], v=_temp[1], w=_temp[2], unit='kpc', representation_type='cartesian', frame='galactic')

    @property
    def parameters(self) -> Dict[str, Union[str,float,int]]:
        return self.survey.parameters
    
    @property
    def parameter_mag_color_names(self) -> str:
        return self.parameters[FTTAGS.mag_color_names]
    
    @property
    def parameter_magnitude_name(self) -> str:
        return self.parameter_mag_color_names.split(',')[0]

    @property
    def parameter_abs_mag_hi(self) -> str:
        return self.parameters[FTTAGS.abs_mag_lim_hi]

    @property
    def parameter_app_mag_hi(self) -> str:
        return self.parameters[FTTAGS.app_mag_lim_hi]

    @cached_property
    def __ebfs_partitions(self) -> Dict[int, Dict[str, NDArray]]:
        return_dict = {}
        for ebf_path in self._ebfs:
            if ebf_path.exists():
                ebf_name: str = ebf_path.name
                ebf_str: str = str(ebf_path.resolve())
                indices_generator = (item
                                     for item in pd.DataFrame(
                                         ebf.read(ebf_str, f"/{self._partitionid}")
                                         ).groupby([0]).indices.items() # groupby ignores ids that don't appear
                                     if item[0]>=0)  # hard-coded condition to only take partition id >= 0
                for i, ind in indices_generator:
                    if i in return_dict:
                        return_dict[i][ebf_name] = ind
                    else:
                        return_dict[i] = {ebf_name: ind}
            else:
                raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        if not return_dict:
            return_dict[-1] = {}
        return return_dict

    @cached_property
    def __ebfs_part_lengths(self) -> Dict[int, Dict[str, int]]:
        return {i: {ebf_name: len(indices)
                    for ebf_name,indices in indices_per_ebf.items()}
                for i,indices_per_ebf in self.__ebfs_partitions.items()}

    @cached_property
    def __partitions_lengths(self) -> Dict[int, int]:
        return {i: sum(part_lengths_in_ebfs.values()) for i, part_lengths_in_ebfs in self.__ebfs_part_lengths.items()}

    @cached_property
    def __ebfs_part_slices(self) -> Dict[int, Dict[str, List[slice]]]:
        return {i: {ebf_name: [slice(start, stop)
                    for start, stop in zip(
                        [indices[0]]+indices[where_slice_change+1].tolist(),
                        (indices[where_slice_change]+1).tolist() + [indices[-1]+1]
                        )]
                    for ebf_name,indices in indices_per_ebf.items()
                    if (where_slice_change:=np.where(np.diff(indices)>1)[0]) is not None}
                for i,indices_per_ebf in self.__ebfs_partitions.items()}

    @cached_property
    def __vaex_partitions(self) -> Dict[int, NDArray]:
        return self._vaex[self._partitionid].to_pandas_series().to_frame().groupby([0]).indices

    @cached_property
    def __vaex_partition_slices(self) -> Dict[int, slice]:
        return {i: slice(indices[0], indices[-1]+1) for i,indices in self.__vaex_partitions.items()}

    @property
    def _vaex(self):
        if self.__vaex is None:
            raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        else:
            return self.__vaex

    @property
    def _vaex_per_partition(self):
        if self.__vaex_per_partition is None:
            raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        else:
            return self.__vaex_per_partition

    @property
    def _path(self):
        raise NotImplementedError
        if self.__path is None:
            return self.__hdf5
        else:
            return self.__path

    @property
    def _file_base(self):
        return self.output_dir / self.output_name

    @property
    def _ebf_glob_pattern(self):
        return self.__name_with_ext('.*.ebf')
    
    @property
    def _ebf_glob(self):
        _temp = self._ebf_glob_pattern
        return _temp.parent.glob(_temp.name)
    
    @cached_property
    def _ebfs(self):
        ebfs = list(self._ebf_glob)
        if not ebfs:
            raise RuntimeError("Galaxia didn't return any ebf output files")
        return ebfs
    
    def __clear_ebfs(self, force: bool = False) -> None:
        for ebf in self._ebf_glob:
            if True if force else input(f"You are about to remove {ebf}, do I proceed? [y/N] ") == 'y':
                ebf.unlink()
    
    @property
    def _hdf5_glob_pattern(self):
        return self.__name_with_ext('.*.h5')
    
    @cached_property
    def _hdf5s(self):
        pattern = self._hdf5_glob_pattern
        partitions = self.__ebfs_partitions
        length_tags = len(str(max(tuple(partitions.keys())+(0,))))
        return {i: pattern.parent / pattern.name.replace('*',f"{i:0{length_tags}d}") for i in partitions.keys()}

    @classmethod
    def __singlethread_flush_extra_columns_to_hdf5(cls, vaex_df: vaex.DataFrame, hdf5_file: pathlib.Path, with_columns: Optional[Iterable] = (), verbose: bool = True) -> None:  # temporary until vaex supports it
        _flush_extra_columns_to_hdf5(vaex_df, hdf5_file, with_columns, verbose)

    def flush_extra_columns_to_hdf5(self, with_columns: Optional[Iterable] = ()) -> None:  # temporary until vaex supports it
        """
            Flush the dataframe new columns to its backend memory-mapped file

            Parameters
            ----------
            with_columns : iterable
                If given an iterable structure of existing column keys, the
                flushing will also overwrite those in the backend file with
                their current in-memory values. Default to an empty tuple.
        """
        # with pathos.pools.ProcessPool(self.__max_pp_workers) as executor:  # credit to https://github.com/uqfoundation/pathos/issues/158#issuecomment-449636971
        #     # Submit tasks to the executor
        #     futures = [executor.apipe(_flush_extra_columns_to_hdf5, vaex_df, hdf5_file, with_columns, self.verbose)
        #                for _, hdf5_file, vaex_df in common_entries(self._hdf5s, self._vaex_per_partition)]
        #     # Collect the results
        #     _ = [future.get() for future in futures]
        with concurrent.futures.ThreadPoolExecutor(max_workers=self._max_pp_workers) as executor:  # credit to https://www.squash.io/how-to-parallelize-a-simple-python-loop/
            # Submit tasks to the executor
            futures = [executor.submit(_flush_extra_columns_to_hdf5, vaex_df, hdf5_file, with_columns, self.verbose)
                       for _, hdf5_file, vaex_df in common_entries(self._hdf5s, self._vaex_per_partition)]
            # Collect the results
            _ = [future.result() for future in concurrent.futures.as_completed(futures)]
        # self.__reload_vaex()
        # gc.collect()

    def __reload_vaex(self) -> None:
        if self.__vaex is not None:
            self.__vaex.close()
        if self._hdf5s.values():
            self.__vaex = vaex.open_many(map(str,self._hdf5s.values()))
        else:
            raise RuntimeError("Corrupted HDF5 internal dictionary")
        if self.__vaex_per_partition is not None and not self._pp_auto_flush:
            for i in self.__vaex_per_partition:
                self.__vaex_per_partition[i].close()
        self.__vaex_per_partition = {i: vaex.open(str(hdf5_file)) for i, hdf5_file in self._hdf5s.items()}
        gc.collect()


Output.__init__.__doc__ = Output.__init__.__doc__.format(_output_properties=''.join(
                                                [f"\n            * {desc} via key{'' if isinstance(key, str) else 's'} ``{str(key).replace(chr(39),'')}``"
                                                    for key, desc in Output._export_properties.union(Output._postprocess_properties)]),
                                                         _optional_properties=''.join(
                                                [f"\n            * {desc} via key{'' if isinstance(key, str) else 's'} ``{str(key).replace(chr(39),'')}``"
                                                    for key, desc in Output._all_optional_properties]))


if __name__ == '__main__':
    raise NotImplementedError()
