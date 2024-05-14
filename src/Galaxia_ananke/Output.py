#!/usr/bin/env python
"""
Contains the Output class definition

Please note that this module is private. The Output class is
available in the main ``Galaxia`` namespace - use that instead.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Tuple, List, Dict
from numpy.typing import NDArray, ArrayLike
from warnings import warn
from functools import cached_property
import concurrent.futures
import pathlib
import itertools
import numpy as np
import h5py as h5
import ebf
import vaex
import pandas as pd
from astropy import units, coordinates
from astropy.utils import classproperty

from .constants import *
from .templates import *
from .defaults import *
from .utils import CallableDFtoNone, common_entries
from .photometry.PhotoSystem import PhotoSystem
from . import Input

if TYPE_CHECKING:
    from . import Survey

__all__ = ['Output']


def shift_g_lon(lon): # restrict longitude values to be within (-180,180)
    return -((-lon+180)%360-180)


class Output:
    _position_prop = (('px', 'py', 'pz'), "Position coordinates in kpc")
    _velocity_prop = (('vx', 'vy', 'vz'), "Velocity coordinates in km/s")
    _celestial_prop = (('ra', 'dec'), "Celestial equatorial coordinates in degrees")
    _galactic_prop = (('glon', 'glat'), "Celestial galactic coordinates in degrees")
    _distance_prop = ('rad', "Distance in kpc")
    _modulus_prop = ('dmod', "Distance modulus in magnitude units")
    _trgbmass_prop = ('mtip', "Tip of the Red Giant Branch stellar mass in solar masses")
    _currentmass_prop = ('mact', "Current stellar mass in solar masses")
    _zamsmass_prop = ('smass', "Zero Age Main Sequence stellar mass in solar masses")
    _age_prop = ('age', "Stellar ages in years and decimal logarithmic scale")
    _surfacegravity_prop = ('grav', "Surface gravity in CGS units and decimal logarithmic scale")
    _metallicity_prop = ('feh', "Stellar metallicity [Fe/H] in dex relative to solar")
    _temperature_prop = ('teff', "Surface temperature in Kelvin and decimal logarithmic scale")
    _luminosity_prop = ('lum', "Stellar luminosity in solar luminosities and decimal logarithmic scale")
    _parentindex_prop = Input._parentindex_prop
    _partitionindex_prop = Input._partitionindex_prop
    _particleflag_prop = ('partid', "Flag = 1 if star not at center of its parent particle")
    _parallax_prop = ('pi', "Parallax in milliarcseconds")
    _propermotion_prop = (('mura', 'mudec'), "Equatorial proper motions in milliarcseconds per year")
    _galacticpropermotion_prop = (('mul', 'mub'), "Galactic proper motions in milliarcseconds per year")
    _radialvelocity_prop = ('vr', "Radial velocity in km/s")
    _vaex_under_list = ['_repr_html_']
    def __init__(self, survey: Survey, parameters: dict) -> None:
        """
            Driver to exploit the output of Galaxia.

            Call signature::

                output = Output(survey, parameters)

            Parameters
            ----------
            survey : :obj:`Survey`
                Survey object that returned this output.
            
            parameters : dict
                Dictionary all of parameters passed by Survey that were used
                to generate this output.
            
            Notes
            -----
            An Output object almost behaves as a vaex DataFrame, also please
            consult vaex online tutorials for more hands-on information:
                https://vaex.io/docs/tutorial.html

            The DataFrame represents the catalogue with columns corresponding
            to properties of the synthetic stars. Those include the photometric
            magnitudes per filter, with each filter identified by a key in the
            lowercase format "photosys_filtername" where photo_sys corresponds
            to the photometric system and filtername corresponds to a filter
            name of that system. With those are also always included the
            following properties: {_output_properties}
            
            Additionally, depending on what optional properties were provided
            with the input particle data, the output can also include the
            following properties: {_optional_properties}
        """
        self.__survey = survey
        self.__parameters = parameters
        self.__vaex = None
        self.__vaex_per_partition = None
        self.__path = None
        self.__clear_ebfs(force=True)

    @classproperty
    def _export_properties(cls):
        return {
            cls._position_prop,
            cls._velocity_prop,
            cls._celestial_prop,
            cls._galactic_prop,
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
            cls._parallax_prop,
            cls._propermotion_prop,
            cls._galacticpropermotion_prop,
            cls._radialvelocity_prop
            }
    
    @classproperty
    def _all_optional_properties(cls):
        return Input._optional_properties - {cls._parentindex_prop, cls._partitionindex_prop}
    
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
        return cls._surfacegravity_prop[0]
    
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
    
    def __getitem__(self, item):
        return self._vaex[item]
    
    def __setitem__(self, item, value):
        self._vaex[item] = value
    
    def __getattr__(self, item):
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
        return tuple(k if k != 'id' else 'satid' for k in self.survey.input.optional_keys())

    def __ebf_to_hdf5_older(self):
        warn('This method is deprecated and does nothing at this time, this will be removed in future versions', DeprecationWarning, stacklevel=2)
        return
        hdf5_file = self.__hdf5
        with h5.File(hdf5_file, 'w') as f5:
            for k in self.export_keys:
                # print(f"Exporting {k}...")
                f5.create_dataset(name=k, data=ebf.read(str(self.__ebf), f"/{k}"))
            print(f"Exported the following quantities to {hdf5_file}")
            print(list(f5.keys()))
        self.__vaex = vaex.open(hdf5_file)

    def __ebf_to_hdf5_old(self):
        warn('This method is deprecated and does nothing at this time, this will be removed in future versions', DeprecationWarning, stacklevel=2)
        return
        for i, hdf5_file, partition_slices, partition_indices in common_entries(self._hdf5s, self.__ebf_part_slices, self.__ebf_partitions):
            with h5.File(hdf5_file, 'w') as f5:
                for k in self.export_keys:
                    # print(f"Exporting {k}...")
                    # f5.create_dataset(name=k, data=ebf.read_ind(str(self._ebf), f"/{k}", partition_indices))
                    data = np.zeros(partition_indices.shape[0],
                                    dtype=ebf.read_ind(str(self.__ebf), f"/{k}", [0]).dtype)
                    head = 0
                    for p_slice in partition_slices:
                        data[head:(head:=head+p_slice.stop-p_slice.start)] = ebf.read(str(self.__ebf),
                                                                                      f"/{k}",
                                                                                      begin=p_slice.start,
                                                                                      end=p_slice.stop)
                    f5.create_dataset(name=k,
                                      data=data)
                print(f"Exported the following quantities to {hdf5_file} for partition {i}")
                print(list(f5.keys()))
        self.__vaex = vaex.open_many(map(str,self._hdf5s.values()))

    def _ebf_to_hdf5(self) -> None:
        ebfs: List[pathlib.Path] = self._ebfs
        export_keys: Tuple[str] = self.export_keys
        with concurrent.futures.ThreadPoolExecutor() as executor:  # credit to https://www.squash.io/how-to-parallelize-a-simple-python-loop/
            # Submit tasks to the executor
            futures = [executor.submit(self.__singlethread_ebf_to_hdf5, i, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs, ebfs, export_keys)
                       for i, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs in common_entries(self._hdf5s, self.__ebfs_part_slices, self.__ebfs_part_lengths)]
            # Collect the results
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
        self.__vaex = vaex.open_many(map(str,self._hdf5s.values()))
        self.__vaex_per_partition = [vaex.open(str(hdf5_file)) for hdf5_file in self._hdf5s.values()]

    @classmethod
    def __singlethread_ebf_to_hdf5(cls, i: int, hdf5_file: pathlib.Path,
                                   part_slices_in_ebfs: Dict[str, List[slice]],
                                   part_lengths_in_ebfs: Dict[str, int],
                                   ebfs: List[pathlib.Path], export_keys: Tuple[str]) -> None:
        n_ebfs: int              = len(ebfs)
        data_length: int         = sum(part_lengths_in_ebfs.values())
        ebfs_slices: Dict[slice] = {ebf_path.name: slice(bounds[0],bounds[1])
                                    for ebf_path, bounds in zip(ebfs, 
                                                                np.repeat(np.cumsum(
                                                                    [0]+[part_lengths_in_ebfs[ebf_path.name] for ebf_path in ebfs]
                                                                            ),
                                                                            [1]+(n_ebfs-1)*[2]+[1]
                                                                            ).reshape((n_ebfs,2)))}
        ebf_sorter: NDArray      = (i + np.arange(n_ebfs)) % n_ebfs
        i_ebf: int               = ebf_sorter[0]
        first_ebf_str: str       = str(ebfs[i_ebf].resolve())
        with h5.File(hdf5_file, 'w') as f5:
            f5datasets = {name: f5.create_dataset(name=name,
                                                    shape=(data_length,),
                                                    dtype=ebf.read_ind(first_ebf_str, f"/{name}", [0]).dtype)
                            for name in export_keys}
            for ebf_path in ebfs[i_ebf:]+ebfs[:i_ebf]:
                ebf_name: str            = ebf_path.name
                ebf_str: str             = str(ebf_path.resolve())
                f5data_slice: slice      = ebfs_slices[ebf_name]
                part_slices: List[slice] = part_slices_in_ebfs[ebf_name]
                for name in export_keys:
                    head = f5data_slice.start
                    for p_slice in part_slices:
                        f5datasets[name][head:(head:=head+p_slice.stop-p_slice.start)] = ebf.read(
                            ebf_str, f"/{name}", begin=p_slice.start, end=p_slice.stop
                            )
                print(f"Exported the following quantities from {ebf_path} to {hdf5_file} for partition {i}")
                print(list(f5.keys()))

    ### DEFINING POST PROCESSING PIPELINES BELOW # TODO consider a PostProcess class that runs postprocess pipeline at __call__ and holds flush_with_columns

    @classmethod
    def __pp_convert_cartesian_to_galactic(cls, df: pd.DataFrame) -> None:
        """
        converts positions & velocities from mock catalog Cartesian coordinates (relative to solar position) 
        into Galactic coordinates, assuming Sun is on -x axis (use rotateStars)
        """
        gc = coordinates.Galactic(u = df[cls._pos[0]].to_numpy()*units.kpc,
                                  v = df[cls._pos[1]].to_numpy()*units.kpc,
                                  w = df[cls._pos[2]].to_numpy()*units.kpc,
                                  U = df[cls._vel[0]].to_numpy()*units.km/units.s,
                                  V = df[cls._vel[1]].to_numpy()*units.km/units.s,
                                  W = df[cls._vel[2]].to_numpy()*units.km/units.s,
                                  representation_type = coordinates.CartesianRepresentation,
                                  differential_type   = coordinates.CartesianDifferential)
        df[cls._gal[0]] = shift_g_lon(gc.spherical.lon.value)
        df[cls._gal[1]] = gc.spherical.lat.value
        df[cls._rad]    = gc.spherical.distance.value
        ####################################
        df[cls._mugal[0]] = gc.sphericalcoslat.differentials['s'].d_lon_coslat.value
        df[cls._mugal[1]] = gc.sphericalcoslat.differentials['s'].d_lat.value
        df[cls._vr]       = gc.sphericalcoslat.differentials['s'].d_distance.value

    @classmethod
    def __pp_convert_galactic_to_icrs(cls, df: pd.DataFrame) -> None:
        """
        converts PMs in galactic coordinates (mulcosb, mub) in arcsec/yr (as output by Galaxia)
        to ra/dec in mas/yr (units of output catalog)
        """
        c = coordinates.Galactic(l               = df[cls._gal[0]].to_numpy()*units.degree,
                                 b               = df[cls._gal[1]].to_numpy()*units.degree,
                                 distance        = df[cls._rad].to_numpy()*units.kpc,
                                 pm_l_cosb       = df[cls._mugal[0]].to_numpy()*units.mas/units.yr,
                                 pm_b            = df[cls._mugal[1]].to_numpy()*units.mas/units.yr,
                                 radial_velocity = df[cls._vr].to_numpy()*units.km/units.s)
        c_icrs = c.transform_to(coordinates.ICRS())
        df[cls._cel[0]] = c_icrs.ra.value
        df[cls._cel[1]] = c_icrs.dec.value
        ####################################
        df[cls._mu[0]]  = c_icrs.pm_ra_cosdec.to(units.mas/units.yr).value
        df[cls._mu[1]]  = c_icrs.pm_dec.to(units.mas/units.yr).value
    
    @classmethod
    def __pp_convert_icrs_to_galactic(cls, df: pd.DataFrame) -> None:
        """
        converts PMs from ICRS coordinates (muacosd, mudec) to Galactic (mul, mub)
        input and output in mas/yr for PMs and degrees for positions
        also exports the galactic lat and longitude
        """
        c = coordinates.ICRS(ra           = df[cls._cel[0]].to_numpy()*units.degree,
                             dec          = df[cls._cel[1]].to_numpy()*units.degree,
                             pm_ra_cosdec = df[cls._mu[0]].to_numpy()*units.mas/units.yr,
                             pm_dec       = df[cls._mu[1]].to_numpy()*units.mas/units.yr)

        c_gal = c.transform_to(coordinates.Galactic())
        df[cls._gal[0]]   = shift_g_lon(c_gal.l.value)
        df[cls._gal[1]]   = c_gal.b.value
        df[cls._mugal[0]] = c_gal.pm_l_cosb.to(units.mas/units.yr).value
        df[cls._mugal[1]] = c_gal.pm_b.to(units.mas/units.yr).value

    @classmethod
    def __pp_last_conversions(cls, df: pd.DataFrame) -> None:
        df[cls._pi]   = 1.0/df[cls._rad]  # parallax in mas (from distance in kpc)
        df[cls._teff] = 10**df[cls._teff]  #Galaxia returns log10(teff/K)
        df[cls._lum]  = 10**df[cls._lum]  #Galaxia returns log10(lum/lsun)

    def apply_post_process_pipeline_and_flush(self, post_process: CallableDFtoNone, *args, flush_with_columns=(), hold_flush: bool = False):
        post_process(self._vaex, *args)
        if not(hold_flush):
            self.flush_extra_columns_to_hdf5(with_columns=flush_with_columns)

    def _post_process(self) -> None:
        self._pp_convert_cartesian_to_galactic()
        self._pp_convert_galactic_to_icrs()
        self._pp_last_conversions()

    def _pp_convert_cartesian_to_galactic(self) -> None:
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_cartesian_to_galactic, flush_with_columns=self._gal+(self._rad,))

    def _pp_convert_galactic_to_icrs(self) -> None:
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_galactic_to_icrs, flush_with_columns=self._cel)
    
    def _pp_convert_icrs_to_galactic(self) -> None:
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_icrs_to_galactic, flush_with_columns=self._gal+self._mugal)

    def _pp_last_conversions(self) -> None:
        self.apply_post_process_pipeline_and_flush(self.__pp_last_conversions, flush_with_columns=(self._teff, self._lum))

    def __name_with_ext(self, ext):
        name_base = self._file_base
        return name_base.parent / f"{name_base.name}{ext}"
    
    def save(self, path):  # TODO Gotta update this
        """
            Save output to new path
        """
        raise NotImplementedError
        old_path = self._path
        self.__path = pathlib.Path(path)
        self._vaex.close()
        old_path.rename(self._path)
        self.__vaex = vaex.open(self._path)

    @property
    def survey(self):
        return self.__survey
    
    @property
    def photosystems(self):
        return self.survey.photosystems

    @property
    def isochrones(self):
        warn('This property will be deprecated, please use instead property photosystems', DeprecationWarning, stacklevel=2)
        return self.photosystems

    @property
    def export_keys(self) -> Tuple[str]:
        return self._make_export_keys(self.photosystems, extra_keys=self._make_input_optional_keys())
    
    @property
    def catalogue_keys(self) -> Tuple[str]:
        return self._make_catalogue_keys(self.photosystems, extra_keys=self._make_input_optional_keys())
    
    @property
    def output_dir(self):
        return pathlib.Path(self._parameters[FTTAGS.output_dir])

    @property
    def output_name(self):
        return f"{self.survey.surveyname}.{self.survey.inputname}"

    @property
    def rsun_skycoord(self):
        _temp = [self._parameters[k] for k in FTTAGS.rSun]
        return coordinates.SkyCoord(u=_temp[0], v=_temp[1], w=_temp[2], unit='kpc', representation_type='cartesian', frame='galactic')

    @property
    def _parameters(self):
        return self.__parameters
    
    @cached_property
    def __ebf_partitions(self) -> Dict[int, NDArray]:
        warn('This property is deprecated as remnant of the single ebf output implementation, this will be removed in future versions', DeprecationWarning, stacklevel=2)
        if self.__ebf.exists():
            return pd.DataFrame(ebf.read(str(self.__ebf), f"/{self._partitionid}")).groupby([0]).indices
        else:
            raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")

    @cached_property
    def __ebfs_partitions(self) -> Dict[int, Dict[str, NDArray]]:
        return_dict = {}
        for ebf_path in self._ebfs:
            if ebf_path.exists():
                ebf_name: str = ebf_path.name
                ebf_str: str = str(ebf_path.resolve())
                for i, ind in pd.DataFrame(ebf.read(ebf_str, f"/{self._partitionid}")).groupby([0]).indices.items():
                    if i in return_dict:
                        return_dict[i][ebf_name] = ind
                    else:
                        return_dict[i] = {ebf_name: ind}
            else:
                raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        return return_dict

    @property
    def __ebfs_part_lengths(self) -> Dict[int, Dict[str, int]]:
        return {i: {ebf_name: len(indices)
                    for ebf_name,indices in indices_per_ebf.items()}
                for i,indices_per_ebf in self.__ebfs_partitions.items()}

    @cached_property
    def __ebf_part_slices(self) -> Dict[int, List[slice]]:
        warn('This property is deprecated as remnant of the single ebf output implementation, this will be removed in future versions', DeprecationWarning, stacklevel=2)
        return {i: [slice(start, stop)
                    for start, stop in zip(
                        [indices[0]]+indices[where_slice_change+1].tolist(),
                        (indices[where_slice_change]+1).tolist() + [indices[-1]+1]
                        )]
                for i,indices in self.__ebf_partitions.items()
                if (where_slice_change:=np.where(np.diff(indices)>1)[0]) is not None}

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
    
    @cached_property
    def __ebf(self):
        warn('This property is deprecated as remnant of the single ebf/single hdf5 output implementation, this will be removed in future versions', DeprecationWarning, stacklevel=2)
        return next(self._ebf_glob)
    
    @property
    def _ebf_glob_pattern(self):
        return self.__name_with_ext('.*.ebf')
    
    @property
    def _ebf_glob(self):
        _temp = self._ebf_glob_pattern
        return _temp.parent.glob(_temp.name)
    
    @cached_property
    def _ebfs(self):
        return list(self._ebf_glob)
    
    def __clear_ebfs(self, force: bool = False) -> None:
        for ebf in self._ebf_glob:
            if True if force else input(f"You are about to remove {ebf}, do I proceed? [y/N] ") == 'y':
                ebf.unlink()

    @cached_property
    def __hdf5(self):
        warn('This property is deprecated as remnant of the single ebf/single hdf5 output implementation, this will be removed in future versions', DeprecationWarning, stacklevel=2)
        return self.__name_with_ext('.h5')
    
    @property
    def _hdf5_glob_pattern(self):
        return self.__name_with_ext('.*.h5')
    
    @cached_property
    def _hdf5s(self):
        pattern = self._hdf5_glob_pattern
        partitions = self.__ebfs_partitions
        length_tags = len(str(max(partitions.keys())))
        return {i: pattern.parent / pattern.name.replace('*',f"{i:0{length_tags}d}") for i in partitions.keys()}
    
    def __flush_extra_columns_to_hdf5_old(self, with_columns=()):  # temporary until vaex supports it
        warn('This method is deprecated and does nothing at this time, this will be removed in future versions', DeprecationWarning, stacklevel=2)
        return
        hdf5_file = self.__hdf5
        old_column_names = set(vaex.open(hdf5_file).column_names)
        with h5.File(hdf5_file, 'r+') as f5:
            extra_columns = [k for k in set(self.column_names)-old_column_names if not k.startswith('__')]
            for k in extra_columns:
                f5.create_dataset(name=k, data=self[k].to_numpy())
            if extra_columns:
                print(f"Exported the following quantities to {hdf5_file}")
                print(extra_columns)
            for k in with_columns:
                f5[k][...] = self[k].to_numpy()
            if with_columns:
                print(f"Overwritten the following quantities to {hdf5_file}")
                print(with_columns)
        self.__vaex = vaex.open(hdf5_file)

    def flush_extra_columns_to_hdf5(self, with_columns=()):  # temporary until vaex supports it
        old_column_names = set(vaex.open(str(self._hdf5s[0])).column_names)
        extra_columns = [k for k in set(self.column_names)-old_column_names if not k.startswith('__')]
        for i, hdf5_file, vaex_slice in common_entries(self._hdf5s, self.__vaex_partition_slices):
            with h5.File(hdf5_file, 'r+') as f5:
                for k in extra_columns:
                    f5.create_dataset(name=k, data=self[vaex_slice][k].to_numpy())
                if extra_columns:
                    print(f"Exported the following quantities to {hdf5_file} for partition {i}")
                    print(extra_columns)
                for k in with_columns:
                    f5[k][...] = self[vaex_slice][k].to_numpy()
                if with_columns:
                    print(f"Overwritten the following quantities to {hdf5_file} for partition {i}")
                    print(with_columns)
        self.__vaex = vaex.open_many(map(str,self._hdf5s.values()))


Output.__init__.__doc__ = Output.__init__.__doc__.format(_output_properties=''.join(
                                                [f"\n                 -{desc} via key `{str(key)}`"
                                                    for key, desc in Output._export_properties.union(Output._postprocess_properties)]),
                                                         _optional_properties=''.join(
                                                [f"\n                 -{desc} via key `{str(key)}`"
                                                    for key, desc in Output._all_optional_properties]))


if __name__ == '__main__':
    pass
