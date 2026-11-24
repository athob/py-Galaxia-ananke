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
Contains the Survey class definition

Please note that this module is private. The Survey class is
available in the main ``galaxia_ananke`` namespace - use that instead.
"""
from __future__ import annotations
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Callable, Optional, Union, Tuple, List, Set, Dict, Iterable, Iterator, Sequence, cast
from numpy.typing import NDArray
from warnings import warn
import re
import shlex
import pathlib
from pprint import PrettyPrinter

import numpy as np

from ._constants import *
from ._templates import *
from ._defaults import *
from .utils import CallableDFtoInt, execute, format_docstring, lexicalorder_dict, mark_metadata_prop, collect_metadata_marked_properties, hash_iterable
from . import photometry
from .photometry.PhotoSystem import PhotoSystem
from .Output import Output

if TYPE_CHECKING:
    from .Input import Input

__all__ = ['Survey']


@collect_metadata_marked_properties
class Survey:
    """
    Driver that runs the backend on an input, in one or several photometric systems, and gives the output that it makes.

    The survey is made by ``make_survey``, which gives its parameters their values: those that are asked for are read from
    the survey afterwards (``fileparam`` for the parameters of the backend).
    """
    _metadata: Dict[str, Any]  # a property that the decorator above adds
    def __init__(self, input: Input, photo_sys: Union[str,List[str]] = DEFAULT_PSYS, surveyname: str = DEFAULT_SURVEYNAME, verbose: bool = True) -> None:
        """
            Driver to exploit the input object and run Galaxia with it.

            Call signature::

                survey = Survey(input,
                                photo_sys={DEFAULT_PSYS},
                                surveyname='{DEFAULT_SURVEYNAME}')

            Parameters
            ----------
            input : :obj:`Input`
                Input object storing the particle data.
            
            photo_sys : string or list
                Name(s) of the photometric system(s) Galaxia should use to
                generate the survey. Default to {DEFAULT_PSYS}.
                Available photometric systems can be found with the photometry
                submodule - please refer to its documentation for further
                details.

            surveyname : string
                Optional name Galaxia should use for the output files. Default
                to '{DEFAULT_SURVEYNAME}'.

            verbose : bool
                Whether the survey prints what it does. Default to True.
        """
        self.__surveyname: str = surveyname
        self.__input: Input = input
        self.__photosystems: List[PhotoSystem] = self.prepare_photosystems(photo_sys)
        self.__verbose: bool = verbose
        self.__fileparam: Optional[MappingProxyType[str, Union[str,float,int]]] = None
        self.__extraparam: Optional[MappingProxyType[str, Union[str,float,int]]] = None
        self.__output: Optional[Output] = None

    format_docstring(__init__, DEFAULT_SURVEYNAME=DEFAULT_SURVEYNAME, DEFAULT_PSYS=DEFAULT_PSYS)
    
    def __repr__(self) -> str:
        """
        Give the representation of the survey.

        Returns
        -------
        representation : str
            The class, with the name of the survey and its photometric systems.
        """
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['surveyname', 'photo_sys']])
        return f'{cls}({description})'

    @classmethod
    def prepare_photosystems(cls, photo_sys: Union[str,List[str]]) -> list[PhotoSystem]:
        """
        Give the photometric systems that are asked for.

        Parameters
        ----------
        photo_sys : str or list of str
            The key of a photometric system, ``category/name``, or the keys of several.

        Returns
        -------
        photosystems : list of PhotoSystem
            The systems, in the order that was given.

        Raises
        ------
        KeyError
            If a key is not the one of an available photometric system.
        """
        if isinstance(photo_sys, str):
            photo_sys = [photo_sys]
        return [photometry.available_photo_systems[psys] for psys in photo_sys]

    @classmethod
    def set_isochrones_from_photosys(cls, photo_sys: Union[str,List[str]]) -> list[PhotoSystem]:
        """
        Give the photometric systems that are asked for, as ``prepare_photosystems`` does: this is the deprecated name of
        it.

        Parameters
        ----------
        photo_sys : str or list of str
            The key of a photometric system, ``category/name``, or the keys of several.

        Returns
        -------
        photosystems : list of PhotoSystem
            The systems, in the order that was given.
        """
        warn('This class method will be deprecated, please use instead class method prepare_photosystems', DeprecationWarning, stacklevel=2)
        return cls.prepare_photosystems(photo_sys)

    def _prepare_survey_parameters_and_output(self, cmd_magnames: Union[str,Dict[str,str]], n_gens: Sequence[int],
                                              **kwargs: Any) -> None:
        """
        Check the parameters of the survey, set them, and make the output that the survey will fill.

        The parameters of the survey are those of the backend (the defaults, the parameters of the survey itself, and those
        that were given), and the generations and the photometric systems after the first, which identify the survey without
        being parameters of the backend.

        Parameters
        ----------
        cmd_magnames : str or dict
            The magnitude and the colors of the selection in the color-magnitude diagram, as
            ``PhotoSystem.check_cmd_magnames`` takes them.
        n_gens : sequence of int
            The numbers of the generations.
        **kwargs
            The parameters of the survey that are not the defaults.

        Raises
        ------
        ValueError
            If there is no generation, or if the names of ``cmd_magnames`` are not magnitudes of the first photometric
            system.
        TypeError
            If the number of neighbours is given, which is the one of the input, or a parameter is not one of the survey.
        """
        if len(n_gens) == 0:
            raise ValueError("n_gens must hold at least one generation (or be an integer above 0), but none is given")
        if FTTAGS.nres in kwargs:
            raise TypeError(f"The number of neighbours of a survey is that of its input ({self.ngb}), it is not a survey parameter: "
                            f"give ngb to the input and not {FTTAGS.nres} to the survey")
        unknown = sorted(set(kwargs) - set(DEFAULTS_FOR_PARFILE))
        if unknown:
            raise TypeError(f"Unknown survey parameter(s) {unknown}: the survey parameters are {sorted(set(DEFAULTS_FOR_PARFILE) - {FTTAGS.nres})}")
        photosys = self.photosystems[0]
        checked_cmd_magnames: str = photosys.check_cmd_magnames(cmd_magnames)
        parameters: Dict[str, Union[str,float,int]] = DEFAULTS_FOR_PARFILE.copy()
        own_parameters: Dict[str, Union[str,float,int]] = {
            FTTAGS.photo_categ: photosys.category, FTTAGS.photo_sys: photosys.name,
            FTTAGS.mag_color_names: checked_cmd_magnames, FTTAGS.nres: self.ngb}
        parameters.update(**own_parameters, **kwargs)
        n_gen_tag_length = len(str(max(n_gens)))
        self.__fileparam = MappingProxyType(parameters)
        self.__extraparam = MappingProxyType({
            **{f"n_gen_{i:0{n_gen_tag_length}d}": n for i,n in enumerate(n_gens)},
            **{k: v for n,PS in enumerate(self.photosystems[1:], start=1)
                    for k,v in zip(FTTAGS.append_photo(n), PS.categ_and_name)}
                    })
        self.__output = Output(self)

    def _write_parameter_file(self) -> Tuple[pathlib.Path, Dict[str, Union[str,float,int]]]:
        """
        Write the parameter file of the backend, in the directory of the input.

        When the survey caches its results, a file that already holds this text is not written again.

        Returns
        -------
        parfile : pathlib.Path
            The file.
        parameters : dict
            The parameters of the survey that are in the file.
        """
        parameters: Dict[str, Union[str,float,int]] = self.fileparam
        surveyname_hash: str = self.surveyname_hash
        parfile: pathlib.Path = self.inputdir / PARFILENAME_TEMPLATE.substitute({FTTAGS.name: surveyname_hash})  # TODO make temporary? create a global record of temporary files?
        parfile_text: str = PARFILE_TEMPLATE.substitute({FTTAGS.output_file: surveyname_hash, **parameters})
        if ((parfile.read_text() != parfile_text # proceed if parfile_text is not in parfile,
            if parfile.exists()                  # only if parfile exist,
            else True)                           # otherwise proceed if doesn't exist
            if self.caching else True):          # -> proceed anyway if self.caching is False
            parfile.write_text(parfile_text)
        return parfile, parameters

    def _run_survey(self, parfile: pathlib.Path, n_gens: Iterable[int], max_gen_workers: int) -> None:
        """
        Run the backend, a command for each generation, all at the same time.

        Parameters
        ----------
        parfile : pathlib.Path
            The parameter file.
        n_gens : iterable of int
            The numbers of the generations.
        max_gen_workers : int
            The number of workers for the generations, which the commands do not use: all of them are started at once.
        """
        cmds = [RUN_TEMPLATE.substitute(**{
            CTTAGS.hdim_block : '' if self.hdim is None
                                else HDIMBLOCK_TEMPLATE.substitute(**{CTTAGS.hdim: self.hdim}),
            CTTAGS.nfile      : shlex.quote(self.inputname_hash),
            CTTAGS.ngen       : ngen,
            CTTAGS.parfile    : shlex.quote(str(parfile))
        }) for ngen in n_gens]
        execute(cmds, max_workers=max_gen_workers, verbose=self.verbose)

    def _append_survey(self, photosystem: PhotoSystem, max_gen_workers: Optional[int]) -> None:
        """
        Append a photometric system to the files that the backend made, a command for each file, all at the same time.

        Parameters
        ----------
        photosystem : PhotoSystem
            The system, which is not the first one.
        max_gen_workers : int or None
            The number of workers, which the commands do not use. None for the number of the files.
        """
        if max_gen_workers is None:
            max_gen_workers = len(list(self.__ebf_output_files_glob))
        cmds = [APPEND_TEMPLATE.substitute(**{
            CTTAGS.pcat     : shlex.quote(photosystem.category),
            CTTAGS.psys     : shlex.quote(photosystem.name),
            CTTAGS.filename : shlex.quote(str(filename))
        }) for filename in self.__ebf_output_files_glob]
        execute(cmds, max_workers=max_gen_workers, verbose=self.verbose)

    @staticmethod
    def _fsample_of_one_generation(sampling_factor: float, n_gens: Sequence[int]) -> float:
        """
        Give the sampling factor of one generation, for a sampling factor of the whole survey.

        The factor is shared among the generations: a single generation keeps it as it was given.

        Parameters
        ----------
        sampling_factor : float
            The sampling factor of the whole survey.
        n_gens : sequence of int
            The numbers of the generations.

        Returns
        -------
        fsample : float
            The factor that the backend is given, which is the sampling factor of each generation.
        """
        n = len(n_gens)
        # a single generation (or none, which is refused later) keeps the value as it was given
        return sampling_factor / n if n > 1 else sampling_factor

    def _vanilla_survey(self, cmd_magnames: Union[str,Dict[str,str]] = DEFAULT_CMD,
                              sampling_factor: Optional[float] = None, input_sorter: Optional[NDArray[np.int_]] = None,
                              n_jobs: Optional[int] = None, n_gens: Union[int, Sequence[int]] = (0,),
                              max_gen_workers: Optional[int] = None, fsample: Optional[float] = None, **kwargs: Any) -> None:
        """
        Run the backend on the input: make the parameters and the files, run the generations, and append the other
        photometric systems.

        Each step that is done once is skipped when the survey caches its results and the state file says that it was done.
        The keywords that ``make_survey`` does not take are given to this method: ``make_survey`` documents them.

        Parameters
        ----------
        cmd_magnames : str or dict, optional
            The magnitude and the colors of the selection in the color-magnitude diagram. Default to ``DEFAULT_CMD``.
        sampling_factor : float, optional
            The sampling factor of the whole survey. Default to 1.
        input_sorter : numpy.ndarray, optional
            The indices that order the particles of the input. Default to the order of their partition ids.
        n_jobs : int, optional
            The deprecated name of ``n_gens``: an integer here replaces it.
        n_gens : int or sequence of int, optional
            The number of generations, or their numbers. Default to a single generation, the number 0.
        max_gen_workers : int, optional
            Not implemented: a warning says so when it is given, and the generations all run at once.
        fsample : float, optional
            The deprecated sampling factor of one generation: it cannot be given with ``sampling_factor``.
        **kwargs
            The other parameters of the survey, which are those of the backend: see ``make_survey``.

        Raises
        ------
        TypeError
            If ``sampling_factor`` and ``fsample`` are both given, and for what ``_prepare_survey_parameters_and_output``
            refuses.
        ValueError
            If there is no generation.
        """
        if isinstance(n_jobs, int):
            n_gens = n_jobs
            warn('Keyword argument n_jobs will be deprecated, please use instead keyword argument n_gens. Consider also reading doc regarding keyword argument max_pp_workers.', DeprecationWarning, stacklevel=2)
        if isinstance(n_gens, int):
            n_gens = range(n_gens)
        if max_gen_workers is None:
            max_gen_workers = len(n_gens)
        else:
            warn('The keyword argument max_gen_workers is currently not implemented.', stacklevel=2)
        if fsample is not None:
            if sampling_factor is not None:
                raise TypeError("Give either sampling_factor, the sampling factor of the whole survey, "
                                "or fsample, that of one generation (deprecated), and not both")
            warn('Keyword argument fsample is deprecated and will be removed: it is the sampling factor of one single '
                 'generation, so that a survey of n_gens generations samples n_gens times fsample. Please use instead '
                 'keyword argument sampling_factor, the sampling factor of the whole survey, which is divided among '
                 'the generations.', DeprecationWarning, stacklevel=3)
        else:
            fsample = self._fsample_of_one_generation(1 if sampling_factor is None else sampling_factor, n_gens)
        self.input.input_sorter = input_sorter
        self._prepare_survey_parameters_and_output(cmd_magnames, n_gens, fsample=fsample, **kwargs)
        inputname, parfile, for_parfile = self.input.prepare_input(self)
        #
        self.check_state_before_running(description='run_survey_complete')(self._run_survey)(parfile, n_gens=n_gens, max_gen_workers=max_gen_workers)
        for photosystem in self.photosystems[1:]:
            self.check_state_before_running(description=f'append_{photosystem.name}_complete', level=1)(self._append_survey)(photosystem, max_gen_workers=max_gen_workers)

    def make_survey(self, *, verbose: bool = True, partitioning_rule: Optional[CallableDFtoInt] = None,
                             max_pp_workers: int = 1, pp_auto_flush: bool = True, **kwargs: Any) -> Output:
        """
            Driver to exploit the input object and run Galaxia with it.
            
            Call signature::
            
                output = self.make_survey(cmd_magnames= '{DEFAULT_CMD}' ,
                                          sampling_factor=1, verbose=True, **kwargs)
            
            Parameters
            ----------
            cmd_magnames : string
                Names of the filters Galaxia should use for the color-
                magnitude diagram box selection. The given string must meet
                the following format::

                    "band1,band2-band3"
                
                where ``band1`` is the magnitude filter and ``(band2, band3)``
                are the filters that define the ``band2-band3`` color index.
                The filter names must correspond to filters that are part of
                the first chosen photometric system in photo_sys. Default to
                ``'{DEFAULT_CMD}'``
            
            sampling_factor : float
                Sampling factor of the whole synthetic star survey, over all
                its independent generations (see n_gens): 1 returns a full
                sample, 0.5 a 50% undersampling and 2 a 200% oversampling.
                Galaxia is given this factor divided by the number of
                generations, which is its own parameter fsample: that of one
                generation. Default to 1.

            fsample : float
                Deprecated, use sampling_factor instead. The sampling factor
                of one single generation, as Galaxia understands it: a survey
                of n_gens generations samples n_gens times fsample.

            input_sorter : array_like of int
                The indices that put the particles of the input, and their
                kernels, in the order in which they are given to Galaxia and
                hashed. Default to the order of the partition ids of the
                particles.
            
            n_gens, n_jobs : int or iterable of int
                Number of independent catalog generations ran in parallel. Can
                also receive an iterable containing each generation number to
                run in parallel. Default to 1. Usage of n_jobs is deprecated
                and will be removed.

            max_gen_workers : int
                CURRENTLY NOT PROPERLY IMPLEMENTED
                Maximum number of workers to parallelize the initial catalog
                generations. Default to the number of independent generations
                in n_gens.
            
            max_pp_workers : int
                Maximum number of workers to parallelize the post-processing
                pipelines after the initial catalog generation. Default to 1.
                The workers stay alive once the survey is made, and the next
                survey that asks for the same number reuses them: see the
                notes below.
            
            pp_auto_flush : bool
                Whether the post-processing pipelines are run by the workers
                on the HDF5 file of each partition, which they update as they
                go (True), or on the data of the partitions in the memory of
                this process, in threads, which is written in the files
                afterwards (False). Default to True.
            
            verbose : bool
                Verbose boolean flag to allow pipeline to print what it's doing
                to stdout. Default to True.
            
            partitioning_rule : callable
                A function that gives the partition of each star of the
                catalogue. It takes the table of the stars, with the columns
                of the catalogue that it uses, and gives an integer for each
                row: the number of the partition of the star. A negative
                number puts the star in no partition, and leaves it out of the
                catalogue. The rule is first applied to an empty table, to
                know the columns that it uses. Without one, a star is in the
                partition of its parent particle, which the key
                ``partitionid`` of the particles gives (0 for all of them if
                there is none: a single partition). Default to None.
            
            parfile : string
                Name of file where Input should save the parameters for
                Galaxia. Default to '{DEFAULT_PARFILE}'
            
            output_dir : string or pathlib.Path
                Path to directory where to save the input/output files of
                Galaxia. Default to '{TTAGS_output_dir}'
            
            app_mag_lim_lo, app_mag_lim_hi, abs_mag_lim_lo, abs_mag_lim_hi, color_lim_lo, color_lim_hi : float
                These allow to specify the limits of the chosen color-magnitude
                diagram box selection (``lo`` for lower and ``hi`` for upper).
                ``app_mag``, ``abs_mag`` and ``color`` represent respectively
                limits in apparent magnitudes, absolute magnitudes and color
                index. Default values follow those set in the dictionary::
                {DEFAULT_CMD_BOX} 
            
            rSun0, rSun1, rSun2 : float
                Coordinates for the observer position in kpc. Respectively
                default to::

                    {TTAGS_rSun0}, {TTAGS_rSun1} & {TTAGS_rSun2}
            
            vSun0, vSun1, vSun2 : float
                Coordinates for the observer velocity in km/s. Respectively
                default to::

                    {TTAGS_vSun0}, {TTAGS_vSun1} & {TTAGS_vSun2}
            
            r_max, r_min : float
                Extent of the shell of radii from observer location within
                which particles should be considered by Galaxia. Respectively
                default to::

                    {TTAGS_r_max} & {TTAGS_r_min}
            
            rand_seed : int
                Seed to be used by Galaxia's pseudorandom number generator.
                Default to {TTAGS_rand_seed}

            nstart : int
                Index at which to start indexing synthetic stars. Default
                to {TTAGS_nstart}

            longitude, latitude : float
                Currently not implemented. Respectively default to::
                
                    {TTAGS_longitude} & {TTAGS_latitude}
            
            star_type : int
                Currently not implemented. Default to {TTAGS_star_type}

            geometry_opt : int
                Currently not implemented. Default to {TTAGS_geometry_opt}

            survey_area : float
                Currently not implemented. Default to {TTAGS_survey_area}
                
            pop_id : int
                Currently not implemented. Default to {TTAGS_pop_id}

            warp_flare_on : int
                Currently not implemented. Default to {TTAGS_warp_flare_on}

            photo_error : int
                Currently not implemented. Default to {TTAGS_photo_error}

            Returns
            -------
            output : :obj:`Output`
                Handler with utilities to utilize the output survey and its
                data.

            Notes
            -----
            The workers of ``max_pp_workers`` are processes that ``pathos``
            starts, by ``fork`` on Linux, and keeps in a cache: they stay
            alive and idle after ``make_survey`` returns, until the Python
            process ends, and keep their memory meanwhile (about 150 MB each
            after a small survey in our measurements).

            It is recommended to ask for the same ``max_pp_workers`` in all
            the surveys of a process: a different number starts another pool
            of workers, which stays alive next to the first one and keeps its
            memory too.

            Workers started by ``fork`` are copies of a process that has
            threads (those of ``vaex``, among others), which Python
            discourages. ``multiprocess.set_start_method('forkserver')``
            starts them without copying it, with the guard below.

            A script that calls ``make_survey`` should guard its entry point
            with ``if __name__ == '__main__':`` whenever the workers are not
            started by ``fork``, as on macOS or after
            ``multiprocess.set_start_method('forkserver')``: those start
            methods run the script again in each worker.
        """  # TODO Move documentation around to the subroutines where they should go
        self.verbose = verbose
        self._vanilla_survey(**kwargs)
        self.output.read_galaxia_output(partitioning_rule, max_pp_workers, pp_auto_flush)
        self.output.post_process_output()
        return self.output

    format_docstring(make_survey,
                     DEFAULT_CMD=DEFAULT_CMD,
                     DEFAULT_PARFILE=DEFAULT_PARFILE,
                     DEFAULT_CMD_BOX=('\n'+PrettyPrinter(width=60).pformat(DEFAULT_CMD_BOX)
                                      ).replace('\n','\n                  '),
                     **{f"TTAGS_{key}": val for key,val in DEFAULTS_FOR_PARFILE.items()})

    @property
    def has_no_fileparam(self) -> bool:
        """Whether the parameters of the survey are not set yet, which they are when the survey is made."""
        return self.__fileparam is None

    @property
    def fileparam(self) -> Dict[str, Union[str,float,int]]:
        """
        The parameters of the survey that are in the parameter file of the backend, by name.

        Reading it raises RuntimeError if the survey has not been made yet.
        """
        if self.has_no_fileparam:
            raise RuntimeError("Survey hasn't been made yet, run method `make_survey` first")
        else:
            return dict(cast("MappingProxyType[str, Union[str,float,int]]", self.__fileparam))

    @property
    def _extraparam(self) -> MappingProxyType[str, Union[str,float,int]]:
        """
        The parameters that identify the survey without being in the parameter file: the numbers of its generations, and the
        category and the name of its photometric systems after the first.

        Reading it raises RuntimeError if the survey has not been made yet.
        """
        if self.has_no_fileparam:
            raise RuntimeError("Survey hasn't been made yet, run method `make_survey` first")
        return cast("MappingProxyType[str, Union[str,float,int]]", self.__extraparam)

    @property
    def parameters(self) -> Dict[str, Union[str,float,int]]:
        """All the parameters of the survey: those of the parameter file, and those that only identify the survey."""
        return {**self.fileparam, **self._extraparam}

    @property
    def n_gens(self) -> List[int]:
        """The numbers of the generations of the survey, in order."""
        return list(lexicalorder_dict({
            int(n_gen[0]): value
            for key, value in self._extraparam.items()
            if (n_gen:=re.findall(r"n_gen_(\d*)", key))
            }).values())

    @property
    @mark_metadata_prop
    def sampling_factor(self) -> float:
        """
        The sampling factor of the whole survey: the factor of each generation, that the backend is given,
        times the number of generations.
        """
        return len(self.n_gens)*cast(float, self.fileparam[FTTAGS.fsample])

    @property
    def fsample(self) -> float:
        """The sampling factor of the whole survey. This is deprecated: it is ``sampling_factor``, and a warning says so."""
        value = self.sampling_factor  # refuses first if the survey was not made
        warn('Survey.fsample is deprecated and will be removed, please use instead Survey.sampling_factor: '
             'the sampling factor of the whole survey.', DeprecationWarning, stacklevel=2)
        return value

    @property
    def _surveyhash(self) -> bytes:  # not cached: the parameters change each time the survey is made
        """
        The hash of the values of the parameters of the survey, in the alphabetical order of their names.

        It is not kept: the parameters change each time the survey is made.
        """
        return hash_iterable(map(lambda el: str(el).encode(HASH_ENCODING),
                                 lexicalorder_dict(self.parameters).values()))

    @property
    @mark_metadata_prop
    def hash(self) -> str:
        """The hash of the parameters of the survey, as a string."""
        return self._surveyhash.decode()

    @property
    @mark_metadata_prop
    def surveyname(self) -> str:
        """The name of the survey."""
        return self.__surveyname
    
    @property
    def append_hash(self) -> bool:
        """Whether the names of the files get the beginning of the hash: it is the choice of the input."""
        return self.input.append_hash

    @property
    def surveyname_hash(self) -> str:
        """
        The name of the survey, followed by an underscore and the first 7 characters of its hash
        if the names of the files get it.
        """
        return self.surveyname + (f"_{self.hash[:7]}" if self.append_hash else "")
    
    @property
    def metadata(self) -> Dict[str, Any]:
        """The metadata of the survey: the values of its properties that are marked for it."""
        return self._metadata

    @property
    def input(self) -> Input:
        """The input that the survey runs the backend on."""
        return self.__input

    @property
    def photosystems(self) -> List[PhotoSystem]:
        """The photometric systems of the survey, the first one being the one that the backend runs on."""
        return self.__photosystems

    @property
    def isochrones(self) -> List[PhotoSystem]:
        """The photometric systems of the survey. This is deprecated: it is ``photosystems``, and a warning says so."""
        warn('This property will be deprecated, please use instead property photosystems', DeprecationWarning, stacklevel=2)
        return self.photosystems

    @property
    def photo_sys(self) -> Set[str]:
        """The keys of the photometric systems of the survey, as a set."""
        return {photosystem.key for photosystem in self.photosystems}

    @property
    def verbose(self) -> bool:
        """Whether the survey prints what it does. It can be set, and ``None`` leaves it as it is."""
        return self.__verbose

    @verbose.setter
    def verbose(self, value: Optional[bool]) -> None:
        if value is not None:
            self.__verbose = value

    @property
    def has_no_output(self) -> bool:
        """Whether the survey has no output yet, which it has when it is made."""
        return self.__output is None

    @property
    def output(self) -> Output:
        """
        The output of the survey, which handles the catalogue.

        Reading it raises RuntimeError if the survey has not been made yet.
        """
        if self.has_no_output:
            raise RuntimeError("Survey hasn't been made yet, run method `make_survey` first")
        else:
            return cast(Output, self.__output)
    
    @property
    def caching(self) -> bool:
        """Whether the results are cached: it is the choice of the input."""
        return self.input.caching

    @property
    def hdim(self) -> int:
        """The dimension of the kernels, as the input says it."""
        return self.input.hdim
    
    @property
    def inputname_hash(self) -> str:
        """The name of the files of the input, with the hash: the one of the input."""
        return self.input.name_hash
    
    @property
    def inputdir(self) -> pathlib.Path:
        """The directory of the files of the input."""
        return self.input._input_dir
    
    @property
    def ngb(self) -> int:
        """The number of neighbours of the particles: the one of the input."""
        return self.input.ngb

    @property
    def check_state_before_running(self) -> Callable[..., Callable[[Callable[..., None]], Callable[..., None]]]:
        """The function of the output that makes a decorator which runs a function only if its work is not done yet."""
        return self.output.check_state_before_running

    @property
    def __ebf_output_files_glob(self) -> Iterator[pathlib.Path]:
        """The files that the backend wrote, as the output finds them."""
        return self.output._ebf_glob


if __name__ == '__main__':
    raise NotImplementedError()
