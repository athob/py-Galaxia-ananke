#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later).
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
from typing import TYPE_CHECKING, Any, Callable, Hashable, Optional, Union, Tuple, List, Set, Dict, Iterable, Iterator, cast
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
import vaex.progress

from .__metadata__ import *
from ._constants import *
from ._templates import *
from ._defaults import *
from .utils import classproperty, format_docstring, in_description_order, CallableDFtoNone, CallableDFtoInt, RecordingDataFrame, State, common_entries, metadata_dicts_to_consolidated_json
from .photometry.PhotoSystem import PhotoSystem
from .Input import Input

if TYPE_CHECKING:
    from .Survey import Survey

__all__ = ['Output']

_Property = Tuple[Union[str, Tuple[str, ...]], str]  # the key, or the keys, of a property, and its description


def shift_g_lon(lon: Any) -> Any: # restrict longitude values to be within (-180,180)
    """
    Shift longitudes to within -180 (excluded) and 180 (included) degrees.

    Parameters
    ----------
    lon : float or numpy.ndarray
        The longitudes in degrees, for instance between 0 and 360.

    Returns
    -------
    lon : float or numpy.ndarray
        The same angles, in the interval: 190 gives -170, and 180 stays 180.
    """
    return -((-lon+180)%360-180)


def _flush_extra_columns_to_hdf5(vaex_df: vaex.DataFrame, hdf5_file: pathlib.Path,
                                 with_columns: Iterable[str] = (),
                                 update_metadata: Optional[Dict[str, Union[str,float,int]]] = None,
                                 verbose: bool = True) -> None:  # temporary until vaex supports it
    """
    Write the new columns of a table in the HDF5 file of its partition, and optionally overwrite the ones that are there.

    This is temporary, until vaex supports it. A column is new if the file does not have it, and the columns whose names
    start with two underscores are not written.

    Parameters
    ----------
    vaex_df : vaex.dataframe.DataFrame
        The table of the partition, with its new columns.
    hdf5_file : pathlib.Path
        The file, which is opened to write in.
    with_columns : iterable of str, optional
        Existing columns to overwrite with the values that the table has now. Default to none.
    update_metadata : dict, optional
        Attributes to set in the file, by name: a value that is not a boolean, a number or a string is a function, which is
        given the attributes of the file and gives the value. Default to none.
    verbose : bool, optional
        Whether to print which columns are written and which are overwritten. Default to True.
    """
    _temp = vaex.open(hdf5_file, group=STARCATALOG_GROUP)
    old_column_names = set(_temp.column_names)
    _temp.close()
    extra_columns = [k for k in set(vaex_df.column_names)-old_column_names if not k.startswith('__')]
    with_columns = list(set(with_columns) - set(extra_columns))
    with h5.File(hdf5_file, 'r+') as f5:
        for k in extra_columns:
            f5[STARCATALOG_GROUP].create_dataset(name=k, data=vaex_df[k].to_numpy())
        if verbose and extra_columns:
            print(f"Exported the following quantities to {hdf5_file}")
            print(extra_columns)
        for k in with_columns:
            f5[STARCATALOG_GROUP][k][...] = vaex_df[k].to_numpy()
        if verbose and len(with_columns):
            print(f"Overwritten the following quantities to {hdf5_file}")
            print(with_columns)
        if isinstance(update_metadata, dict):
            for k, v in update_metadata.items():
                f5.attrs[k] = v if isinstance(v, (bool, int, float, str)) else v(f5.attrs)


def _decorate_post_processing(pp: CallableDFtoNone, hdf5_path_input: bool = False, flush_with_columns: Iterable[str] = (),
                              max_thread_workers: Optional[int] = None,
                              update_metadata: Optional[Dict[str, Union[str,float,int]]] = None,
                              verbose: bool = True) -> Callable[..., None]:
    """
    Make the function that a worker runs to apply a pipeline to partitions and flush what the pipeline made.

    The function takes a partition, or a list of them, and then the other arguments of the pipeline. A partition is the path
    of its HDF5 file if ``hdf5_path_input`` is set, which the function opens, with its own pool of threads for vaex, and
    closes after the new columns are written in the file. Otherwise it is the table of the partition, which the pipeline
    changes.

    Parameters
    ----------
    pp : callable
        The pipeline, which takes the table and the other arguments, and gives nothing.
    hdf5_path_input : bool, optional
        Whether the partitions are given as the paths of their HDF5 files, which are flushed after the pipeline. Default to
        False.
    flush_with_columns : iterable of str, optional
        The existing columns that are overwritten in the files when they are flushed. Default to none.
    max_thread_workers : int, optional
        The number of threads that vaex uses for a file that the function opens. Default to the number of vaex.
    update_metadata : dict, optional
        Attributes to set in the files when they are flushed (see ``_flush_extra_columns_to_hdf5``). Default to none.
    verbose : bool, optional
        Whether to print what is flushed. Default to True.

    Returns
    -------
    new_pp : callable
        The function.
    """
    def new_pp(*args: Any) -> None:
        first_arg = args[0]
        if not isinstance(first_arg, list):
            first_arg = [first_arg]
        for _temp in first_arg:
            if hdf5_path_input:
                hdf5_file: pathlib.Path = _temp
                old_vaex_main_executor = vaex.dataframe.main_executor
                vaex.dataframe.main_executor = vaex.execution.ExecutorLocal(vaex.multithreading.ThreadPoolIndex(max_workers=max_thread_workers))
                vaex_df: vaex.DataFrame = vaex.open(hdf5_file, group=STARCATALOG_GROUP)
            else:
                vaex_df = _temp
            pp(vaex_df, *args[1:])
            if hdf5_path_input:
                _flush_extra_columns_to_hdf5(vaex_df, hdf5_file, with_columns=flush_with_columns, update_metadata=update_metadata, verbose=verbose)
                vaex_df.close()
                vaex.dataframe.main_executor = old_vaex_main_executor
            gc.collect()
    return new_pp


class Output:
    """
    The catalogue of a survey, as the HDF5 files that hold it, with what makes them from the files of the backend and
    post-processes them.

    It wraps a vaex DataFrame, whose public names it has too. A survey makes it: the documentation of its constructor says
    what it holds.
    """
    _position_prop = (('px', 'py', 'pz'), "Position coordinates in $kpc$")
    _velocity_prop = (('vx', 'vy', 'vz'), "Velocity coordinates in $km/s$")
    _celestial_prop = (('ra', 'dec'), "Celestial equatorial coordinates in $degrees$")
    _galactic_prop = (('glon', 'glat'), "Celestial galactic coordinates in $degrees$")
    _distance_prop = ('rad', "Distance in $kpc$")
    _modulus_prop = ('dmod', "Distance modulus in magnitude units")
    _trgbmass_prop = ('mtip', "Tip of the Red Giant Branch stellar mass in solar masses")
    _currentmass_prop = ('mact', "Current stellar mass in solar masses")
    _initialmass_prop = ('minit', "Zero Age Main Sequence stellar mass in solar masses")
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
            When a pipeline run executes, it persists the resulting mock star
            catalogue to the filesystem as a set of one or more HDF5 files
            rather than holding it entirely in memory; the `Output` object
            provides the lightweight, on-disk interface to these files. Acting
            as a lazy proxy, it internally wraps a ``vaex.DataFrame`` and
            exposes a nearly identical API, facilitating efficient out-of-core
            manipulation of catalogues that may exceed available system RAM.
            For practical usage and advanced operations, consult the official
            ``vaex`` online tutorials:
            
                https://vaex.io/docs/tutorial.html
            
            The wrapped DataFrame represents the catalogue with columns
            corresponding to properties of the stars from the synthetic
            stellar population it simulates.
            
            .. warning:: When generated directly by ``galaxia_ananke``, the
                         catalogue properties reflect the quantities as
                         computed by Galaxia. However, the catalogue can be
                         modified or extended via post-processing routines
                         using the
                         ``apply_post_process_pipeline_and_flush`` method. If
                         the `Output` object originates from other software,
                         post-processing may already have been applied;
                         consult that software's documentation for a complete
                         overview of the catalogue contents.
            
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
            These magnitudes are absolute.

            With those are also always included the following properties:
            {_output_properties}
            
            Additionally, depending on what optional properties were provided
            with the input particle data, the output can also include the
            following properties:
            {_optional_properties}
            
            Internally, the catalogue is stored across one or more HDF5 files
            (single by default, splittable via a ``partitionid`` key of the
            input particles or a custom ``partitioning_rule``; a star whose
            partition is negative is left out of the catalogue). Each file
            contains three root groups: ``{STARCATALOG_GROUP}``, with datasets
            for individual catalogue columns, ``{PARENTPARTICLE_GROUP}``,
            holding parent-particle properties and within-file and cross-file
            demographic statistics, and ``{PARTITIONING_GROUP}``, whose dataset
            ``partition_id`` lists the partitions of the catalogue.
            Root-level HDF5 attributes additionally log pipeline runtime
            metadata (e.g., software versions and run parameters).

            These attributes record where a file comes from: ``Header`` names
            the version of this package that made it, ``Cat_partition_id`` is
            the partition that the file holds, ``All_partition_ids`` says where
            the numbers of all of them are, and the parameters and the metadata
            of the survey are under the prefixes ``SurveyParam_``, ``Input_``
            (the name and the hash of the input) and ``Survey_`` (its name, its
            hash and its sampling factor). ``Galaxia_ebf_output_logs`` is JSON
            text that holds the logs that the backend wrote in its ``.ebf``
            files, one for each file, in a compact form: the items that are the
            same in all of them are under ``common_values``, and the others are
            named in ``varying_keys`` and have their value for each file in
            ``variations``. Among those items, each photometric system of the
            survey has one named ``isochrone_info_`` followed by the name of
            the system: the raw descriptor of its isochrones, with their
            magnitude system, and the raw header of their files, which record
            the choices made when the isochrones were produced.
        """
        self.__survey: Survey = survey
        self.__vaex: Any = None
        self.__path: Optional[pathlib.Path] = None
        self.__make_state()
        if not self.caching:
            self.__clear_ebfs(force=True)
        self._max_pp_workers = 1
        self._pp_auto_flush = True

    class _State(State):
        pass

    @classproperty
    def _export_properties(cls) -> Set[_Property]:
        """
        The properties of the catalogue that come from the files of the backend,
        as pairs of the key, or the keys, and a description.
        """
        return {
            cls._position_prop,
            cls._velocity_prop,
            cls._distance_prop,
            cls._modulus_prop,
            cls._trgbmass_prop,
            cls._currentmass_prop,
            cls._initialmass_prop,
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
    def _postprocess_properties(cls) -> Set[_Property]:
        """
        The properties that the post-processing pipelines add to the catalogue,
        as pairs of the key, or the keys, and a description.
        """
        return {
            cls._celestial_prop,
            cls._galactic_prop,
            cls._parallax_prop,
            cls._propermotion_prop,
            cls._galacticpropermotion_prop,
            cls._radialvelocity_prop
            }
    
    @classproperty
    def _all_optional_properties(cls) -> Set[Tuple[str, str]]:
        """
        The properties that the input may have and the catalogue then has, as pairs of the key and a description.

        They are the optional properties of the input without the indices of the parents, of the partitions and of the
        populations, the last of which is given the key ``satid``.
        """
        return Input._optional_properties \
              - {cls._parentindex_prop, cls._partitionindex_prop, cls._satindex_prop} \
                | {(cls._satindex, cls._satindex_prop[1])}
    
    @classproperty
    def _export_keys(cls) -> Tuple[str, ...]:
        """The keys of the columns that come from the files of the backend."""
        return tuple(sum([(_p[0],) if isinstance(_p[0], str) else _p[0] for _p in cls._export_properties], ()))
    
    @classproperty
    def _postprocess_keys(cls) -> Tuple[str, ...]:
        """The keys of the columns that the post-processing pipelines add."""
        return tuple(sum([(_p[0],) if isinstance(_p[0], str) else _p[0] for _p in cls._postprocess_properties], ()))
    
    @classproperty
    def _pos(cls) -> Tuple[str, str, str]:
        """The keys of the three components of the positions in the catalogue."""
        return cls._position_prop[0]
    
    @classproperty
    def _vel(cls) -> Tuple[str, str, str]:
        """The keys of the three components of the velocities in the catalogue."""
        return cls._velocity_prop[0]
    
    @classproperty
    def _cel(cls) -> Tuple[str, str]:
        """The keys of the right ascension and the declination in the catalogue."""
        return cls._celestial_prop[0]
    
    @classproperty
    def _gal(cls) -> Tuple[str, str]:
        """The keys of the Galactic longitude and latitude in the catalogue."""
        return cls._galactic_prop[0]
    
    @classproperty
    def _rad(cls) -> str:
        """The key of the distances in the catalogue."""
        return cls._distance_prop[0]
    
    @classproperty
    def _dmod(cls) -> str:
        """The key of the distance moduli in the catalogue."""
        return cls._modulus_prop[0]
    
    @classproperty
    def _mtip(cls) -> str:
        """The key of the masses at the tip of the red giant branch in the catalogue."""
        return cls._trgbmass_prop[0]
    
    @classproperty
    def _mact(cls) -> str:
        """The key of the current masses in the catalogue."""
        return cls._currentmass_prop[0]
    
    @classproperty
    def _minit(cls) -> str:
        """The key of the initial masses in the catalogue."""
        return cls._initialmass_prop[0]
    
    @classproperty
    def _age(cls) -> str:
        """The key of the ages in the catalogue."""
        return cls._age_prop[0]
    
    @classproperty
    def _grav(cls) -> str:
        """The key of the surface gravities in the catalogue."""
        return cls._surfacegravity_prop[0]
    
    @classproperty
    def _feh(cls) -> str:
        """The key of the metallicities in the catalogue."""
        return cls._metallicity_prop[0]
    
    @classproperty
    def _teff(cls) -> str:
        """The key of the temperatures in the catalogue."""
        return cls._temperature_prop[0]
    
    @classproperty
    def _lum(cls) -> str:
        """The key of the luminosities in the catalogue."""
        return cls._luminosity_prop[0]
    
    @classproperty
    def _parentid(cls) -> str:
        """The key of the indices of the parent particles in the catalogue."""
        return cls._parentindex_prop[0]
    
    @classproperty
    def _partitionid(cls) -> str:
        """The key of the indices of the partitions in the catalogue."""
        return cls._partitionindex_prop[0]
        
    @classproperty
    def _partid(cls) -> str:
        """The key of the flags of the stars that are not at the center of their parent particle in the catalogue."""
        return cls._particleflag_prop[0]
    
    @classproperty
    def _pi(cls) -> str:
        """The key of the parallaxes in the catalogue."""
        return cls._parallax_prop[0]
    
    @classproperty
    def _mu(cls) -> Tuple[str, str]:
        """The keys of the two components of the equatorial proper motions in the catalogue."""
        return cls._propermotion_prop[0]
    
    @classproperty
    def _mugal(cls) -> Tuple[str, str]:
        """The keys of the two components of the Galactic proper motions in the catalogue."""
        return cls._galacticpropermotion_prop[0]
    
    @classproperty
    def _vr(cls) -> str:
        """The key of the radial velocities in the catalogue."""
        return cls._radialvelocity_prop[0]
    
    def __dir__(self) -> List[str]:
        """
        Give the names of the attributes of the output: its own, and the public ones of the table that it wraps.

        Returns
        -------
        names : list of str
            The names, sorted.
        """
        return sorted({i for i in self.__vaex.__dir__() if not i.startswith('_')}.union(
            super(Output, self).__dir__()).union(
            self._vaex_under_list if self.__vaex is not None else []))

    def __repr__(self) -> str:
        """
        Give the representation of the output: the one of the table that it wraps.

        Returns
        -------
        representation : str
            The representation.
        """
        return repr(self._vaex)
    
    def __getitem__(self, item: str) -> Any:
        """
        Get a column, or a selection, of the table that the output wraps.

        Parameters
        ----------
        item : str
            What vaex takes to get an item: the key of a column, an expression, a mask.

        Returns
        -------
        item : object
            What vaex gives.
        """
        return self._vaex[item]
    
    def __setitem__(self, item: str, value: Any) -> None:
        """
        Set a column of the table that the output wraps, which is a new one if the key is new.

        Parameters
        ----------
        item : str
            The key of the column.
        value : object
            What vaex takes for a column: an array, an expression.
        """
        self._vaex[item] = value

    def __len__(self) -> int:
        """
        Give the number of stars of the catalogue.

        Returns
        -------
        length : int
            The number of rows of the table.
        """
        return len(self._vaex)

    def __getattr__(self, item: str) -> Any:
        """
        Give an attribute of the table that the output wraps, if the output does not have it.

        Only the public ones are given, and ``_repr_html_`` once the table exists: the others are looked for in the output.

        Parameters
        ----------
        item : str
            The name of the attribute.

        Returns
        -------
        attribute : object
            The attribute of the table, or of the output.

        Raises
        ------
        AttributeError
            If neither has it.
        """
        if (item in self.__vaex.__dir__() and not item.startswith('_')) or (item in self._vaex_under_list and self.__vaex is not None):
            return getattr(self.__vaex, item)
        else:
            return self.__getattribute__(item)

    @classmethod
    def _compile_export_mag_names(cls, photosystems: list[PhotoSystem]) -> Tuple[str, ...]:
        """
        Give the keys of the magnitudes of photometric systems in the catalogue.

        Parameters
        ----------
        photosystems : list of PhotoSystem
            The systems.

        Returns
        -------
        keys : tuple of str
            The keys of the magnitudes of each system, one system after the other.
        """
        return tuple(itertools.chain.from_iterable([photosystem.to_export_keys for photosystem in photosystems]))
    
    @classmethod
    def _make_export_keys(cls, photosystems: list[PhotoSystem], extra_keys: Iterable[str] = ()) -> Tuple[str, ...]:
        """
        Give the keys of the columns that come from the files of the backend, for photometric systems.

        Parameters
        ----------
        photosystems : list of PhotoSystem
            The systems, whose magnitudes are among the columns.
        extra_keys : iterable of str, optional
            Other keys. Default to none.

        Returns
        -------
        keys : tuple of str
            The keys of the properties of the backend, of the extra ones and of the magnitudes, each once, in no particular
            order.
        """
        return tuple(set(cls._export_keys).union(extra_keys).union(cls._compile_export_mag_names(photosystems)))

    @classmethod
    def _make_catalogue_keys(cls, photosystems: list[PhotoSystem], extra_keys: Tuple[str, ...] = ()) -> Tuple[str, ...]:
        """
        Give the keys of all the columns of the catalogue, for photometric systems: those of the files of the backend, and
        those that the post-processing adds.

        Parameters
        ----------
        photosystems : list of PhotoSystem
            The systems, whose magnitudes are among the columns.
        extra_keys : tuple of str, optional
            Other keys. Default to none.

        Returns
        -------
        keys : tuple of str
            The keys, each once, in no particular order.
        """
        return cls._make_export_keys(photosystems, extra_keys=cls._postprocess_keys+extra_keys)

    def _make_input_optional_keys(self) -> Tuple[str, ...]:
        """
        Give the optional keys that the input has, as the catalogue has them: the index of the population is ``satid``.

        Returns
        -------
        keys : tuple of str
            The keys.
        """
        return tuple(k if k != self._satindex_prop[0] else self._satindex for k in self.survey.input.optional_keys())

    def _redefine_partitions_in_ebfs(self, partitioning_rule: Optional[CallableDFtoInt]) -> None:
        """
        Give the stars the partition that a rule gives them, in each file of the backend, with a process for each worker.

        Nothing is done if there is no rule.

        Parameters
        ----------
        partitioning_rule : callable or None
            The rule: see ``Survey.make_survey``.
        """
        if partitioning_rule is not None:  # TODO test instead partitioning_rule?
            ebfs: List[pathlib.Path] = self._ebfs
            export_keys: Tuple[str, ...] = self.export_keys
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
    def _singlethread_redefine_partitions(cls, ebf_file: pathlib.Path, partitioning_rule: CallableDFtoInt, export_keys: Tuple[str, ...]) -> None:
        """
        Give the stars of a file of the backend the partition that a rule gives them, and sort them by partition.

        Parameters
        ----------
        ebf_file : pathlib.Path
            The file of the backend.
        partitioning_rule : callable
            The rule: see ``Survey.make_survey``.
        export_keys : tuple of str
            The keys of the columns of the catalogue, which the empty table has.

        Notes
        -----
        The rule is first applied to an empty table to know the columns that it uses, and then to a table of those columns,
        read in the file. The partitions are written in the file, the arrays of stars are put in the order of the partitions
        (the other items of the file are not changed), and the log of the file says that Python changed it.
        """
        dummy_df = RecordingDataFrame([], columns=export_keys, dtype=float)
        _ = partitioning_rule(dummy_df)
        ebf_df = pd.DataFrame({key: ebf.read(str(ebf_file), f"/{key}") for key in dummy_df.record_of_all_used_keys})
        new_partition_id = partitioning_rule(ebf_df)
        del ebf_df
        ebf.update_ind(str(ebf_file), f'/{cls._partitionid}', new_partition_id)
        partition_id_sorter = np.argsort(new_partition_id)
        del new_partition_id
        gc.collect()
        n_stars = len(partition_id_sorter)
        for key in ebf._EbfMap.keys(str(ebf_file)):
            if key not in [b'/log', b'/center'] and not key.startswith(b'/.'):
                values = ebf.read(str(ebf_file), key)
                # only the arrays of stars are reordered: the other items (the descriptions of the isochrones) stay as they are
                if np.ndim(values) > 0 and len(values) == n_stars:
                    ebf.update_ind(str(ebf_file), key, values[partition_id_sorter])
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
        """
        Make the HDF5 file of each partition from the files of the backend, with a process for each worker.

        When the pipelines are not run by the workers on the files (``pp_auto_flush`` is not set), the catalogue is loaded
        again afterwards.
        """
        ebfs: List[pathlib.Path] = self._ebfs
        export_keys: Tuple[str, ...] = self.export_keys
        partition_ids: List[int] = list(self._hdf5s.keys())
        photosystem_names: Tuple[str, ...] = tuple(photosystem.name for photosystem in self.photosystems)
        # with concurrent.futures.ThreadPoolExecutor(max_workers=self._max_pp_workers) as executor:  # credit to https://www.squash.io/how-to-parallelize-a-simple-python-loop/
        #     # Submit tasks to the executor
        #     futures = [executor.submit(self._singlethread_ebf_to_hdf5, i, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs, ebfs, export_keys)
        #                for i, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs in common_entries(self._hdf5s, self.__ebfs_part_slices, self.__ebfs_part_lengths)]
        #     # Collect the results
        #     _ = [future.result() for future in concurrent.futures.as_completed(futures)]
        with pathos.pools.ProcessPool(self._max_pp_workers) as executor:  # credit to https://github.com/uqfoundation/pathos/issues/158#issuecomment-449636971
            # Submit tasks to the executor
            futures = [executor.apipe(self._singlethread_ebf_to_hdf5, i, partition_id, partition_ids, hdf5_file, part_slices_in_ebfs, part_lengths_in_ebfs, ebfs, export_keys, self.all_metadata, self.verbose, photosystem_names)
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
                                  ebfs: List[pathlib.Path], export_keys: Tuple[str, ...],
                                  metadata: Dict[str, Union[str,float,int]],
                                  verbose: bool = True, photosystem_names: Tuple[str, ...] = ()) -> None:
        """
        Make the HDF5 file of a partition from the files of the backend.

        The file has the three groups of the catalogue, and the stars of the partition in each of the files of the backend,
        which are taken in turn, the file that comes first depending on the place of the partition. The attributes give the
        consolidated logs of the backend, a header, the number of the partition, where the numbers of all of them are, and
        the metadata of the survey. The group of the parent particles is left to be filled by the post-processing.

        Parameters
        ----------
        i : int
            The place of the partition among all of them, which sets the file of the backend that comes first.
        partition_id : int
            The number of the partition.
        partition_ids : list of int
            The numbers of all the partitions.
        hdf5_file : pathlib.Path
            The file to write.
        part_slices_in_ebfs : dict
            For each file of the backend, by name, the slices of its stars that are in the partition.
        part_lengths_in_ebfs : dict
            For each file of the backend, by name, the number of its stars that are in the partition.
        ebfs : list of pathlib.Path
            The files of the backend.
        export_keys : tuple of str
            The columns to write.
        metadata : dict
            The attributes to write in the file, by name.
        verbose : bool, optional
            Whether to print what is exported. Default to True.
        photosystem_names : tuple of str, optional
            The photometric systems, whose information on their isochrones the logs of the backend may hold. Default to
            none.
        """
        header_ebf: str          = str(ebfs[0].resolve())
        ebfs = [ebf_path for ebf_path in ebfs if ebf_path.name in part_lengths_in_ebfs]
        n_ebfs: int              = len(ebfs)
        if ebfs:
            i_ebf: int               = i % n_ebfs
            ebfs = ebfs[i_ebf:]+ebfs[:i_ebf]
            header_ebf = str(ebfs[0].resolve())
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
            f5.create_group(STARCATALOG_GROUP)
            f5.create_group(PARENTPARTICLE_GROUP)
            f5.create_group(PARTITIONING_GROUP)
            f5datasets = {name: f5[STARCATALOG_GROUP].create_dataset(
                                                  name=name,
                                                  shape=(data_length,),
                                                  dtype=ebf.getHeader(header_ebf, f"/{name}").get_dtype()
                                                  )
                            for name in export_keys}
            ebf_logs = []
            for ebf_path in ebfs:
                ebf_name: str            = ebf_path.name
                ebf_str: str             = str(ebf_path.resolve())
                ebf_log: Dict[str, str]  = cls.__parse_galaxia_ebf_log(ebf_str, photosystem_names)
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
                    print(list(f5[STARCATALOG_GROUP].keys()))
                ebf_logs.append(ebf_log)
            f5.attrs["Galaxia_ebf_output_logs"] = metadata_dicts_to_consolidated_json(ebf_logs, metadata_name = "Galaxia_ebf_output_logs")
            f5.attrs[METADATA_HEADER_KEY] = (
                f"File generated by Python module {NAME} v{__version__}, "
                f"software available at <{__url__}>."
            )
            f5.attrs[f"Cat_partition_id"] = partition_id
            f5.attrs[f"All_partition_ids"] = f"Stored under dataset {PARTITIONING_GROUP}/partition_id"
            f5[PARTITIONING_GROUP].create_dataset(
                name='partition_id',
                data=partition_ids
                )
            for k, v in metadata.items():
                f5.attrs[k] = v

    @classmethod
    def __parse_galaxia_ebf_log(cls, ebf_str: str, photosystem_names: Tuple[str, ...] = ()) -> Dict[str, str]:
        """
        Read the log that the backend wrote in a file.

        Parameters
        ----------
        ebf_str : str
            The path of the file.
        photosystem_names : tuple of str, optional
            The photometric systems, whose information on their isochrones the backend may have written in an item of its
            own. Default to none.

        Returns
        -------
        log : dict
            The text before and after the parameters, each parameter by its name, and the information on the isochrones of
            each system that the backend gave it.
        """
        log_string: str = str(ebf.read(ebf_str, '/log'))
        log_header, log_body = log_string.split('\n# <parameterfile>\n')
        log_body, log_footer = log_body.split('\n# </parameterfile>\n')
        log_body = log_body.replace('\n# ','\n')
        log_dict: Dict[str, str] = {parts[0]: ' '.join(parts[1:]) for parts in (line.split() for line in log_body.splitlines()) if len(parts)>0}
        # the backend documents each photometric system in an item of its own, `<system>_isochrone_info`: the raw descriptor
        # and the raw header of its isochrones (a backend that does not do it leaves them out, and so does the log)
        isochrone_info: Dict[str, str] = {}
        for name in photosystem_names:
            item = f"/{name.lower()}_isochrone_info"
            if ebf.containsKey(ebf_str, item):
                isochrone_info[f"isochrone_info_{name}"] = str(ebf.read(ebf_str, item))
        return {"log_header": log_header, **log_dict, "log_footer": log_footer, **isochrone_info}

    def read_galaxia_output(self, partitioning_rule: Optional[CallableDFtoInt], max_pp_workers: int, pp_auto_flush: bool) -> None:
        """
        Make the catalogue from the files of the backend: give the stars their partitions, and write the HDF5 files.

        Each step is skipped when the output caches its results and the state file says that it was done.

        Parameters
        ----------
        partitioning_rule : callable or None
            The rule that gives the partitions: see ``Survey.make_survey``.
        max_pp_workers : int
            The number of workers.
        pp_auto_flush : bool
            Whether the pipelines are run by the workers on the files: see ``Survey.make_survey``.
        """
        self._max_pp_workers = max_pp_workers
        self._pp_auto_flush = pp_auto_flush
        self.check_state_before_running(description="redefine_partitions_in_ebfs")(self._redefine_partitions_in_ebfs)(partitioning_rule)
        self.check_state_before_running(description="convert_ebf_to_hdf5", level=1)(self._ebf_to_hdf5)()

    ### DEFINING POST PROCESSING PIPELINES BELOW # TODO consider a PostProcess class that runs postprocess pipeline at __call__ and holds flush_with_columns

    @classmethod
    def __pp_convert_cartesian_to_galactic(cls, df: pd.DataFrame) -> None:
        """
        Convert the positions and the velocities to Galactic coordinates, assuming the Sun is on the -x axis.

        The positions are Cartesian and relative to the solar position. The table is given the longitude and the latitude,
        in degrees, the distance, the proper motions in the directions of longitude and latitude, in milliarcseconds per
        year, and the radial velocity.

        Parameters
        ----------
        df : pandas.DataFrame
            The table, which is changed in place.
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
        Convert Galactic coordinates and proper motions to equatorial ones.

        The table is given the right ascension and the declination, in degrees, and the proper motions in milliarcseconds
        per year.

        Parameters
        ----------
        df : pandas.DataFrame
            The table, which is changed in place.
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
        Convert equatorial coordinates and proper motions to Galactic ones.

        The table is given the longitude and the latitude, in degrees, and the proper motions in milliarcseconds per year.

        Parameters
        ----------
        df : pandas.DataFrame
            The table, which is changed in place.
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
        """
        Make the last conversions: the parallax, in milliarcseconds, from the distance in kiloparsecs, and the temperature
        and the luminosity from their decimal logarithms, which are what the backend gives.

        Parameters
        ----------
        df : pandas.DataFrame
            The table, which is changed in place.
        """
        df[cls._pi]   = 1.0/df[cls._rad]  # parallax in mas (from distance in kpc)
        df[cls._teff] = 10**df[cls._teff]  #Galaxia returns log10(teff/K)
        df[cls._lum]  = 10**df[cls._lum]  #Galaxia returns log10(lum/lsun)

    @property
    def _max_pp_workers(self) -> int:
        """The number of workers that the steps in processes, and the pipelines, use."""
        return self.__max_pp_workers
    
    @_max_pp_workers.setter
    def _max_pp_workers(self, value: int) -> None:
        self.__max_pp_workers: int = value

    @property
    def _pp_auto_flush(self) -> bool:
        """
        Whether the pipelines are run by the workers on the HDF5 files, which they update,
        and not on the tables in memory.
        """
        return self.__pp_auto_flush
    
    @_pp_auto_flush.setter
    def _pp_auto_flush(self, value: bool) -> None:
        self.__pp_auto_flush: bool = value

    @staticmethod
    def __consolidate_partitions_per_process(partitions: NDArray, lengths: NDArray, max_workers: int) -> Dict[Hashable, pd.Index]:
        """
        Share the partitions among the workers, in groups that have about the same number of stars.

        There are at most as many groups as workers, and fewer if there are fewer partitions.

        Parameters
        ----------
        partitions : numpy.ndarray
            The numbers of the partitions.
        lengths : numpy.ndarray
            The number of stars of each partition.
        max_workers : int
            The number of workers.

        Returns
        -------
        groups : dict
            The partitions of each group, as an index, by the number of the group.
        """
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

    def apply_post_process_pipeline_and_flush(self, post_process: CallableDFtoNone, *args: Any,
                                              flush_with_columns: Iterable[str] = (), hold_flush: bool = False,
                                              hold_reload: bool = False,
                                              consolidate_partitions_per_process: bool = False,
                                              update_metadata: Optional[Dict[str, Union[str,float,int]]] = None) -> None:
        r"""
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

            consolidate_partitions_per_process : bool
                Whether to give each worker a group of partitions, of about
                the same total number of stars, instead of a task for each
                partition: there are then as many tasks as workers. This is for
                catalogues of many partitions. Default to False.

            update_metadata : dict
                Attributes to set in the HDF5 files once the pipeline is
                applied, by name. A value that is not a boolean, a number or a
                string is a function: it is given the attributes of the file,
                and gives the value to set. Default to None, which sets
                nothing.
        """
        # post_process(self._vaex, *args)
        vaex_df_or_hdf5_or_list_s: Dict[Any, Any] = self._hdf5s if self._pp_auto_flush else self._vaex_per_partition
        if consolidate_partitions_per_process:
            partitions_per_process: Dict[Hashable, pd.Index] = self.__consolidate_partitions_per_process(
                np.array(list(self.__partitions_lengths.keys())),
                np.array(list(self.__partitions_lengths.values())),
                self._max_pp_workers
                )
            vaex_df_or_hdf5_or_list_s = {process: [vaex_df_or_hdf5_or_list_s[i] for i in partitions]
                                         for process, partitions in partitions_per_process.items()}
        if update_metadata is not None:
            pass # TODO check format of update_metadata ?
        if self._pp_auto_flush:
            with pathos.pools.ProcessPool(self._max_pp_workers) as executor:  # credit to https://github.com/uqfoundation/pathos/issues/158#issuecomment-449636971
                # Submit tasks to the executor
                futures = [executor.apipe(_decorate_post_processing(post_process,
                                                                    self._pp_auto_flush,
                                                                    flush_with_columns=flush_with_columns,
                                                                    max_thread_workers=int(np.ceil((os.cpu_count() or 1)/self._max_pp_workers)),
                                                                    update_metadata=update_metadata,
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
                                                                    update_metadata=update_metadata,
                                                                    verbose=self.verbose),
                                        vaex_df_or_hdf5, *args)
                        for vaex_df_or_hdf5 in vaex_df_or_hdf5_or_list_s.values()]
                # Collect the results
                _ = [future.result() for future in concurrent.futures.as_completed(futures)]
        if not(hold_flush) and not(self._pp_auto_flush):
            self.flush_extra_columns_to_hdf5(with_columns=flush_with_columns, update_metadata=update_metadata)
        if not(hold_reload):
            self.__reload_vaex()

    def post_process_output(self) -> None:
        """
        Make the post-processing that every catalogue gets, and the data of the parent particles.

        The positions and the velocities are converted to Galactic coordinates, these to equatorial ones, and the last
        conversions are made (see ``apply_post_process_pipeline_and_flush``). The catalogue is loaded again, a warning says
        so if it has no star, and the data of the parent particles is written in the files. Each step is skipped when the
        output caches its results and the state file says that it was done.
        """
        self.check_state_before_running(description="pp_cartesian_to_galactic")(self._pp_convert_cartesian_to_galactic)(hold_reload=True)
        self.check_state_before_running(description="pp_galactic_to_icrs", level=1)(self._pp_convert_galactic_to_icrs)(hold_reload=True)
        self.check_state_before_running(description="pp_last_conversions", level=1)(self._pp_last_conversions)(hold_reload=True)
        self.__reload_vaex()
        if len(self._vaex) == 0:
            warn("The survey has no star in any partition, its catalogue is empty: check the limits of the survey (magnitudes, colors), "
                 "its input and its partitioning rule", UserWarning)
        self.check_state_before_running(description="export_parent_data", level=1)(self.__export_parent_particles_data)()

    def _pp_convert_cartesian_to_galactic(self, **kwargs: Any) -> None:
        """
        Apply the conversion of the positions and the velocities to Galactic coordinates to the catalogue, and flush what it
        made.

        Parameters
        ----------
        **kwargs
            Given to ``apply_post_process_pipeline_and_flush``.
        """
        pipeline_name = "convert_cartesian_to_galactic"
        if self.verbose:
            print(f"Running {pipeline_name} post-processing pipeline")
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_cartesian_to_galactic, flush_with_columns=self._gal+(self._rad,), **kwargs)

    def _pp_convert_galactic_to_icrs(self, **kwargs: Any) -> None:
        """
        Apply the conversion of Galactic coordinates to equatorial ones to the catalogue, and flush what it made.

        Parameters
        ----------
        **kwargs
            Given to ``apply_post_process_pipeline_and_flush``.
        """
        pipeline_name = "convert_galactic_to_icrs"
        if self.verbose:
            print(f"Running {pipeline_name} post-processing pipeline")
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_galactic_to_icrs, flush_with_columns=self._cel, **kwargs)
    
    def _pp_convert_icrs_to_galactic(self, **kwargs: Any) -> None:
        """
        Apply the conversion of equatorial coordinates to Galactic ones to the catalogue, and flush what it made.

        Parameters
        ----------
        **kwargs
            Given to ``apply_post_process_pipeline_and_flush``.
        """
        pipeline_name = "convert_icrs_to_galactic"
        if self.verbose:
            print(f"Running {pipeline_name} post-processing pipeline")
        self.apply_post_process_pipeline_and_flush(self.__pp_convert_icrs_to_galactic, flush_with_columns=self._gal+self._mugal, **kwargs)

    def _pp_last_conversions(self, **kwargs: Any) -> None:
        """
        Apply the last conversions: the parallax, the temperature and the luminosity to the catalogue, and flush what it
        made.

        Parameters
        ----------
        **kwargs
            Given to ``apply_post_process_pipeline_and_flush``.
        """
        pipeline_name = "last_conversions"
        if self.verbose:
            print(f"Running {pipeline_name} post-processing pipeline")
        self.apply_post_process_pipeline_and_flush(self.__pp_last_conversions, flush_with_columns=(self._teff, self._lum), **kwargs)

    def __name_with_ext(self, ext: str) -> pathlib.Path:
        """
        Give the path of a file of the output, from its suffix.

        Parameters
        ----------
        ext : str
            The suffix, which starts with a dot.

        Returns
        -------
        path : pathlib.Path
            The path: the base of the files of the output, followed by the suffix.
        """
        name_base = self._file_base
        return name_base.parent / f"{name_base.name}{ext}"
    
    def save(self, path: Union[str, pathlib.Path]) -> None:  # TODO Gotta update this
        """
        Save the output to a new path.

        This is not implemented.

        Parameters
        ----------
        path : str or pathlib.Path
            The new path.

        Raises
        ------
        NotImplementedError
            Always.
        """
        raise NotImplementedError
        old_path = self._path
        self.__path = pathlib.Path(path)
        self._vaex.close()
        old_path.rename(self._path)
        self.__vaex = vaex.open(self._path)

    def __make_state(self) -> None:
        """Make the state of the output, which is kept in a file next to its files."""
        self.__state: Output._State = self._State(self.__name_with_ext('.dummy'), self)

    @property
    def _state(self) -> Output._State:
        """The state of the output, which tells what was done when the results are cached."""
        return self.__state

    @property
    def check_state_before_running(self) -> Callable[..., Callable[[Callable[..., None]], Callable[..., None]]]:
        """The function of the state that makes a decorator which runs a function only if its work is not done yet."""
        return self._state.check_state_file_before_running

    @property
    def caching(self) -> bool:
        """Whether the results are cached: it is the choice of the survey."""
        return self.survey.caching

    @property
    def verbose(self) -> bool:
        """Whether the output prints what it does: it is the choice of the survey."""
        return self.survey.verbose

    @property
    def survey(self) -> Survey:
        """The survey that made the output."""
        return self.__survey
    
    @property
    def photosystems(self) -> List[PhotoSystem]:
        """The photometric systems of the survey."""
        return self.survey.photosystems

    @property
    def magnitude_systems(self) -> Dict[str, Optional[str]]:
        """
        The magnitude system of each magnitude of the catalogue, by its key: 'Vega', 'AB' or 'ST', or None where the descriptor
        of the isochrones of its photometric system does not say it.
        """
        return {key: photosystem.magnitude_system for photosystem in self.photosystems for key in photosystem.to_export_keys}

    def offsets_to_magnitude_system(self, target: str) -> Dict[str, float]:
        """
        Give what to add to each magnitude of the catalogue to have it in a magnitude system.

        It is zero for those that are in it already, whichever system they are in (``magnitude_systems``). The conversion
        between Vega and AB uses the zero points of the SVO filter profile service, which needs the network: see
        ``PhotoSystem.vega_to_ab_offsets``.

        Parameters
        ----------
        target : str
            The magnitude system to convert to: 'Vega' or 'AB', whatever the case.

        Returns
        -------
        offsets : dict
            What to add to each magnitude, by its key in the catalogue.

        Raises
        ------
        ValueError
            If the target is not given or is not a magnitude system, or if the descriptor of the isochrones of a system does
            not say which system its magnitudes are in.
        NotImplementedError
            If the conversion is not one between Vega and AB.
        """
        return {key: offset for photosystem in self.photosystems
                for key, offset in photosystem.offsets_to_magnitude_system(target).items()}

    @property
    def export_keys(self) -> Tuple[str, ...]:
        """
        The keys of the columns that come from the files of the backend, for the survey:
        with the magnitudes and the optional keys of the input.
        """
        return self._make_export_keys(self.photosystems, extra_keys=self._make_input_optional_keys())
    
    @property
    def catalogue_keys(self) -> Tuple[str, ...]:
        """The keys of all the columns of the catalogue, for the survey, once it is post-processed."""
        return self._make_catalogue_keys(self.photosystems, extra_keys=self._make_input_optional_keys())
    
    @property
    def output_dir(self) -> pathlib.Path:
        """The directory of the files of the output."""
        return pathlib.Path(cast("Union[str, pathlib.Path]", self.parameters[FTTAGS.output_dir]))

    @property
    def output_name(self) -> str:
        """The name of the files of the output: the name of the survey, then the one of the input, each with its hash."""
        return f"{self.survey.surveyname_hash}.{self.survey.inputname_hash}"

    @property
    def rsun_skycoord(self) -> coordinates.SkyCoord:
        """The position of the observer, in Galactic Cartesian coordinates, as a sky coordinate."""
        _temp = [self.parameters[k] for k in FTTAGS.rSun]
        return coordinates.SkyCoord(u=_temp[0], v=_temp[1], w=_temp[2], unit='kpc', representation_type='cartesian', frame='galactic')

    @property
    def parameters(self) -> Dict[str, Union[str,float,int]]:
        """The parameters of the survey."""
        return self.survey.parameters
    
    @property
    def parameter_mag_color_names(self) -> str:
        """The magnitude and the colors of the selection in the color-magnitude diagram, as the survey was made with them."""
        return cast(str, self.parameters[FTTAGS.mag_color_names])
    
    @property
    def parameter_magnitude_name(self) -> str:
        """The name of the magnitude of the selection in the color-magnitude diagram."""
        return self.parameter_mag_color_names.split(',')[0]

    @property
    def parameter_abs_mag_hi(self) -> float:
        """The upper limit of the absolute magnitude of the selection."""
        return cast(float, self.parameters[FTTAGS.abs_mag_lim_hi])

    @property
    def parameter_app_mag_hi(self) -> float:
        """The upper limit of the apparent magnitude of the selection."""
        return cast(float, self.parameters[FTTAGS.app_mag_lim_hi])

    @property
    def all_metadata(self) -> Dict[str, Any]:
        """
        The metadata that the HDF5 files have as attributes, by name.

        It is made of the parameters of the survey, and the metadata of the input and of the survey, each with a prefix. A
        value that is not a boolean, a number or a string is given as a string.
        """
        return {
            k: v if isinstance(v, (bool, int, float, str)) else str(v)  # TODO consider Iterable case?
            for k, v in {
                **{f"{PARAMETER_METADATA_PREFIX}{k}": v for k, v in self.parameters.items()},
                **{f"{INPUT_METADATA_PREFIX}{k}": v for k, v in self.survey.input.metadata.items()},
                **{f"{SURVEY_METADATA_PREFIX}{k}": v for k, v in self.survey.metadata.items()},
            }.items()
        }

    @cached_property
    def __ebfs_partitions(self) -> Dict[int, Dict[str, NDArray]]:
        """
        The stars of each partition in each file of the backend: for a partition, the indices of its stars in each file, by
        the name of the file.

        Only the partitions that have a number of zero or more are there, and if there is none the partition -1 is, without
        stars. Reading it raises RuntimeError if a file of the backend does not exist.
        """
        return_dict: Dict[int, Dict[str, NDArray]] = {}
        for ebf_path in self._ebfs:
            if ebf_path.exists():
                ebf_name: str = ebf_path.name
                ebf_str: str = str(ebf_path.resolve())
                indices_generator = (item
                                     for item in pd.DataFrame(
                                         ebf.read(ebf_str, f"/{self._partitionid}")
                                         ).groupby([0]).indices.items() # groupby ignores ids that don't appear
                                     if cast(int, item[0])>=0)  # hard-coded condition to only take partition id >= 0
                for key, ind in indices_generator:
                    i, indices = cast(int, key), cast(NDArray, ind)  # pandas types them widely
                    if i in return_dict:
                        return_dict[i][ebf_name] = indices
                    else:
                        return_dict[i] = {ebf_name: indices}
            else:
                raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        if not return_dict:
            return_dict[-1] = {}
        return return_dict

    @cached_property
    def __ebfs_part_lengths(self) -> Dict[int, Dict[str, int]]:
        """
        The number of stars of each partition in each file of the backend,
        by partition and then by the name of the file.
        """
        return {i: {ebf_name: len(indices)
                    for ebf_name,indices in indices_per_ebf.items()}
                for i,indices_per_ebf in self.__ebfs_partitions.items()}

    @cached_property
    def __partitions_lengths(self) -> Dict[int, int]:
        """The number of stars of each partition, by partition."""
        return {i: sum(part_lengths_in_ebfs.values()) for i, part_lengths_in_ebfs in self.__ebfs_part_lengths.items()}

    @cached_property
    def __ebfs_part_slices(self) -> Dict[int, Dict[str, List[slice]]]:
        """
        The slices of the stars of each partition in each file of the backend, by partition and then by the name of the
        file: one for each run of consecutive stars.
        """
        return {i: {ebf_name: [slice(start, stop)
                    for start, stop in zip(
                        [indices[0]]+indices[where_slice_change+1].tolist(),
                        (indices[where_slice_change]+1).tolist() + [indices[-1]+1]
                        )]
                    for ebf_name,indices in indices_per_ebf.items()
                    if (where_slice_change:=np.where(np.diff(indices)>1)[0]) is not None}
                for i,indices_per_ebf in self.__ebfs_partitions.items()}

    # (nothing uses the two properties below: commented out, with the mark UNUSED, to be removed after v1)
    # UNUSED: @cached_property
    # UNUSED: def __vaex_partitions(self) -> Dict[int, NDArray]:
    # UNUSED:     return self._vaex[self._partitionid].to_pandas_series().to_frame().groupby([0]).indices

    # UNUSED: @cached_property
    # UNUSED: def __vaex_partition_slices(self) -> Dict[int, slice]:
    # UNUSED:     return {i: slice(indices[0], indices[-1]+1) for i,indices in self.__vaex_partitions.items()}

    @property
    def _vaex(self) -> vaex.dataframe.DataFrameLocal:
        """
        The table that the output wraps: the catalogue, with all its partitions.

        Reading it raises RuntimeError if there is none, which is the case of an output that was not made by a survey.
        """
        if self.__vaex is None:
            raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        else:
            return self.__vaex

    @cached_property
    def _vaex_per_partition(self) -> Dict[int, vaex.dataframe.DataFrameLocal]:
        # opened when asked for: with the auto flush the workers open the files by path, and each frame holds
        # 2 descriptors and a mapping per partition in the main process; __reload_vaex forgets them
        """
        The table of each partition, by partition, opened when it is asked for.

        Reading it raises RuntimeError if the output was not made by a survey.
        """
        if self.__vaex is None:
            raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        return {i: vaex.open(str(hdf5_file), group=STARCATALOG_GROUP) for i, hdf5_file in self._hdf5s.items()}

    @property
    def _path(self) -> pathlib.Path:
        """The path of the output. This is not implemented."""
        raise NotImplementedError
        if self.__path is None:
            return self.__hdf5
        else:
            return self.__path

    @property
    def _file_base(self) -> pathlib.Path:
        """The path of the files of the output without their suffix: in the output directory, with the name of the output."""
        return self.output_dir / self.output_name

    @property
    def _ebf_glob_pattern(self) -> pathlib.Path:
        """The pattern of the paths of the files that the backend wrote."""
        return self.__name_with_ext('.*.ebf')
    
    @property
    def _ebf_glob(self) -> Iterator[pathlib.Path]:
        """The files that the backend wrote, as the file system lists them."""
        _temp = self._ebf_glob_pattern
        return _temp.parent.glob(_temp.name)
    
    @staticmethod
    def _ebf_generation_key(path: pathlib.Path) -> Tuple[int, Union[int, str]]:
        """
        Give what puts the files of the backend in the order of their generations.

        The files are named after the output, then the number of the generation, then ``.ebf``.

        Parameters
        ----------
        path : pathlib.Path
            The file.

        Returns
        -------
        key : tuple
            A key to sort by: the numbered generations come first, in numerical order, then the others by name.
        """
        generation = path.stem.rsplit('.', 1)[-1]
        return (0, int(generation)) if generation.isdigit() else (1, generation)

    @cached_property
    def _ebfs(self) -> List[pathlib.Path]:
        """
        The files that the backend wrote, in the order of their generations.

        Reading it raises RuntimeError if there is none.
        """
        ebfs = sorted(self._ebf_glob, key=self._ebf_generation_key)  # not in the order that the file system lists them in, which can change
        if not ebfs:
            raise RuntimeError("Galaxia didn't return any ebf output files")
        return ebfs
    
    def __clear_ebfs(self, force: bool = False) -> None:
        """
        Remove the files that the backend wrote.

        Parameters
        ----------
        force : bool, optional
            Whether to remove them without asking. Default to False, which asks for each file.
        """
        for ebf in self._ebf_glob:
            if True if force else input(f"You are about to remove {ebf}, do I proceed? [y/N] ") == 'y':
                ebf.unlink()
    
    @property
    def _hdf5_glob_pattern(self) -> pathlib.Path:
        """The pattern of the paths of the HDF5 files of the partitions."""
        return self.__name_with_ext('.*.h5')
    
    @cached_property
    def _hdf5s(self) -> Dict[int, pathlib.Path]:
        """
        The HDF5 file of each partition, by partition: the number of the partition is in the name,
        with leading zeros so that they are of one length.
        """
        pattern = self._hdf5_glob_pattern
        partitions = self.__ebfs_partitions
        length_tags = len(str(max(tuple(partitions.keys())+(0,))))
        return {i: pattern.parent / pattern.name.replace('*',f"{i:0{length_tags}d}") for i in partitions.keys()}

    def flush_extra_columns_to_hdf5(self, with_columns: Iterable[str] = (), update_metadata: Optional[Dict[str, Union[str,float,int]]] = None) -> None:  # temporary until vaex supports it
        """
            Flush the dataframe new columns to its backend memory-mapped file

            Parameters
            ----------
            with_columns : iterable
                If given an iterable structure of existing column keys, the
                flushing will also overwrite those in the backend file with
                their current in-memory values. Default to an empty tuple.

            update_metadata : dict
                Attributes to set in the HDF5 files, by name. A value that is
                not a boolean, a number or a string is a function: it is given
                the attributes of the file, and gives the value to set.
                Default to None, which sets nothing.
        """
        # with pathos.pools.ProcessPool(self.__max_pp_workers) as executor:  # credit to https://github.com/uqfoundation/pathos/issues/158#issuecomment-449636971
        #     # Submit tasks to the executor
        #     futures = [executor.apipe(_flush_extra_columns_to_hdf5, vaex_df, hdf5_file, with_columns, update_metadata, self.verbose)
        #                for _, hdf5_file, vaex_df in common_entries(self._hdf5s, self._vaex_per_partition)]
        #     # Collect the results
        #     _ = [future.get() for future in futures]
        with concurrent.futures.ThreadPoolExecutor(max_workers=self._max_pp_workers) as executor:  # credit to https://www.squash.io/how-to-parallelize-a-simple-python-loop/
            # Submit tasks to the executor
            futures = [executor.submit(_flush_extra_columns_to_hdf5, vaex_df, hdf5_file, with_columns, update_metadata, self.verbose)
                       for _, hdf5_file, vaex_df in common_entries(self._hdf5s, self._vaex_per_partition)]
            # Collect the results
            _ = [future.result() for future in concurrent.futures.as_completed(futures)]
        # self.__reload_vaex()
        # gc.collect()

    def __reload_vaex(self) -> None:
        """
        Load the catalogue again from the HDF5 files, and let go of the tables that were opened on them before.

        Raises
        ------
        RuntimeError
            If there is no file.
        """
        if self.__vaex is not None:
            self.__vaex.close()
        if self._hdf5s.values():
            self.__vaex = vaex.open_many(map(str,self._hdf5s.values()), group=STARCATALOG_GROUP)
        else:
            raise RuntimeError("Corrupted HDF5 internal dictionary")
        if '_vaex_per_partition' in self.__dict__ and not self._pp_auto_flush:
            for vaex_df in self._vaex_per_partition.values():
                vaex_df.close()
        self.__dict__.pop('_vaex_per_partition', None)
        gc.collect()

    @property
    def _parent_stats_per_partitionid(self) -> Dict[int, pd.DataFrame]:
        """
        The statistics of the parent particles in each partition, by partition, as tables indexed by the index of the
        parent.

        A table has, for each parent, the number of its stars in the partition and the sum of their initial masses, and
        then, over all the partitions, the number of those that have stars of the parent, and the number and the sum of the
        initial masses of these stars. There are no rows when the catalogue has no star.
        """
        if len(self._vaex) == 0:  # vaex cannot group a table without a row: no parent has a star, and each statistic is an empty column
            no_stars = pd.DataFrame({'star_count': np.empty(0, dtype='int64'), f'{self._minit}_sum': np.empty(0, dtype='float64'),
                                     'parent_partition_count': np.empty(0, dtype='int64'), 'parent_star_count': np.empty(0, dtype='int64'),
                                     f'parent_{self._minit}_sum': np.empty(0, dtype='float64')},
                                    index=pd.Index(np.empty(0, dtype='int64'), name=self._parentid))
            return {partitionid: no_stars for partitionid in self._hdf5s}
        groupby: vaex.dataframe.DataFrameLocal = self._vaex.groupby(
            [self._partitionid, self._parentid],
            agg={'star_count': vaex.agg.count(), f'{self._minit}_sum': vaex.agg.sum(self._minit)},
            sort=True, progress=True)
        parent_groupby: vaex.dataframe.DataFrameLocal = groupby.groupby(
            self._parentid,
            agg={'parent_partition_count': vaex.agg.count(),
                 'parent_star_count': vaex.agg.sum('star_count'),
                 f'parent_{self._minit}_sum': vaex.agg.sum(f'{self._minit}_sum')},
            sort=True, progress=True)
        final_groupby: vaex.dataframe.DataFrameLocal = groupby.join(parent_groupby, on=self._parentid)
        parent_stats: Dict[int, pd.DataFrame] = {
            partitionid: stats.drop(self._partitionid, inplace=True).to_pandas_df().set_index(self._parentid)
            for (partitionid,), stats in final_groupby.groupby(self._partitionid)}
        # vaex keeps the progress tree of its last computation in a module global, which refers to the frames of the
        # groupbys above and, through them, to the open file of every partition: let it go, or a survey that is dropped
        # keeps its files open until another computation of vaex replaces the tree
        vaex.progress._last_progress_tree = None
        return parent_stats

    def __export_parent_particles_data(self) -> None:
        """
        Write, in the file of each partition, the data of the parent particles that have stars in it.

        These are the index of the parent, its properties that the input kept, and the statistics of its stars: they go in
        the group of the parent particles, in the order of the particles of the input.
        """
        for _, hdf5_file, parent_stats in common_entries(self._hdf5s, self._parent_stats_per_partitionid):
            particles_mask: NDArray = np.isin(self.survey.input.particle_parentids, parent_stats.index)
            particles_parentids: NDArray = self.survey.input.particle_parentids[particles_mask]
            reordered_stats: pd.DataFrame = parent_stats.loc[particles_parentids]
            with h5.File(hdf5_file, 'r+') as f5:
                f5[PARENTPARTICLE_GROUP].create_dataset(name=self._parentid, data=particles_parentids)
                for key, array in self.survey.input.parent_particles_data_generator:
                    f5[PARENTPARTICLE_GROUP].create_dataset(name=key, data=array[particles_mask])
                for column_name, column in reordered_stats.items():
                    f5[PARENTPARTICLE_GROUP].create_dataset(name=column_name, data=column.to_numpy())


format_docstring(Output.__init__,
                 _output_properties=''.join(
                     [f"\n            * {desc} via key{'' if isinstance(key, str) else 's'} ``{str(key).replace(chr(39),'')}``"
                      for key, desc in in_description_order(
                          Output._export_properties.union(Output._postprocess_properties))]),
                 _optional_properties=''.join(
                     [f"\n            * {desc} via key{'' if isinstance(key, str) else 's'} ``{str(key).replace(chr(39),'')}``"
                      for key, desc in in_description_order(Output._all_optional_properties)]),
                 STARCATALOG_GROUP=STARCATALOG_GROUP,
                 PARENTPARTICLE_GROUP=PARENTPARTICLE_GROUP,
                 PARTITIONING_GROUP=PARTITIONING_GROUP)


if __name__ == '__main__':
    raise NotImplementedError()
