#!/usr/bin/env python
"""
Contains the Survey class definition

Please note that this module is private. The Survey class is
available in the main ``Galaxia`` namespace - use that instead.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Union, Dict
from warnings import warn
from pprint import PrettyPrinter

from ._constants import *
from ._templates import *
from ._defaults import *
from .utils import execute
from . import photometry
from .photometry.PhotoSystem import PhotoSystem
from .Output import Output

if TYPE_CHECKING:
    from . import Input

__all__ = ['Survey']


class Survey:
    def __init__(self, input: Input, photo_sys=DEFAULT_PSYS, surveyname='survey') -> None:
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
        """
        self.__surveyname = surveyname
        self.__input = input
        self.__photosystems = self.prepare_photosystems(photo_sys)
        self.__output = None

    __init__.__doc__ = __init__.__doc__.format(DEFAULT_SURVEYNAME=DEFAULT_SURVEYNAME,
                                               DEFAULT_PSYS=DEFAULT_PSYS)
    
    def __repr__(self) -> str:
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['surveyname', 'photo_sys']])
        return f'{cls}({description})'

    @classmethod
    def prepare_photosystems(cls, photo_sys: str) -> list[PhotoSystem]:
        if isinstance(photo_sys, str):
            photo_sys = [photo_sys]
        return [photometry.available_photo_systems[psys] for psys in photo_sys]

    @classmethod
    def set_isochrones_from_photosys(cls, photo_sys: str) -> list[PhotoSystem]:
        warn('This class method will be deprecated, please use instead class method prepare_photosystems', DeprecationWarning, stacklevel=2)
        return cls.prepare_photosystems(photo_sys)

    def _run_survey(self, cmd_magnames: Union[str,Dict[str,str]], fsample: float, n_jobs: int = 1, verbose: bool = True, **kwargs) -> None:
        inputname, parfile, for_parfile = self.input.prepare_input(self.photosystems[0], cmd_magnames,
                                                                   output_file=self.surveyname, fsample=fsample, **kwargs)
        self.__output = Output(self, for_parfile)
        cmds = [RUN_TEMPLATE.substitute(**{
            CTTAGS.hdim_block : '' if self.hdim is None
                                else HDIMBLOCK_TEMPLATE.substitute(**{CTTAGS.hdim: self.hdim}),
            CTTAGS.nfile      : self.inputname,
            CTTAGS.ngen       : ngen,
            CTTAGS.parfile    : parfile
        }) for ngen in range(n_jobs)]
        execute(cmds, verbose=verbose)

    def _append_survey(self, photosystem: PhotoSystem, verbose: bool = True) -> None:
        cmds = [APPEND_TEMPLATE.substitute(**{
            CTTAGS.pcat     : photosystem.category,
            CTTAGS.psys     : photosystem.name,
            CTTAGS.filename : filename
        }) for filename in self.__ebf_output_files_glob]
        execute(cmds, verbose=verbose)

    def make_survey(self, cmd_magnames: Union[str,Dict[str,str]] = DEFAULT_CMD, fsample: float = 1, n_jobs: int = 1, verbose: bool = True, **kwargs) -> Output:
        """
            Driver to exploit the input object and run Galaxia with it.
            
            Call signature::
            
                output = self.make_survey(cmd_magnames= '{DEFAULT_CMD}' ,
                                          fsample=1, verbose=True, **kwargs)
            
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
                '{DEFAULT_CMD}'
            
            fsample : float
                Sampling rate from 0 to 1 for the resulting synthetic star
                survey. 1 returns a full sample while any value under returns
                partial surveys. Default to 1.
            
            n_jobs : int
                Number of independent catalog generations ran in parallel.
                Default to 1.
            
            verbose : bool
                Verbose boolean flag to allow pipeline to print what it's doing
                to stdout. Default to True.
            
            parfile : string
                Name of file where Input should save the parameters for
                Galaxia. Default to '{DEFAULT_PARFILE}'
            
            output_dir : string
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
            
            r_maxm, r_min : float
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
        """
        self._run_survey(cmd_magnames, fsample, n_jobs=n_jobs, verbose=verbose, **kwargs)
        for photosystem in self.photosystems[1:]:
            self._append_survey(photosystem, verbose=verbose)
        self.output._ebf_to_hdf5()
        self.output._post_process()
        return self.output

    make_survey.__doc__ = make_survey.__doc__.format(DEFAULT_CMD=DEFAULT_CMD,
                                                     DEFAULT_PARFILE=DEFAULT_PARFILE,
                                                     DEFAULT_CMD_BOX=('\n'+PrettyPrinter(width=60).
                                                                      pformat(DEFAULT_CMD_BOX)).
                                                                     replace('\n','\n                  '),
                                                     **{f"TTAGS_{key}": val
                                                        for key,val in DEFAULTS_FOR_PARFILE.items()})

    @property
    def surveyname(self):
        return self.__surveyname
    
    @property
    def input(self):
        return self.__input

    @property
    def photosystems(self):
        return self.__photosystems

    @property
    def isochrones(self):
        warn('This property will be deprecated, please use instead property photosystems', DeprecationWarning, stacklevel=2)
        return self.photosystems

    @property
    def photo_sys(self):
        return {photosystem.key for photosystem in self.photosystems}

    @property
    def output(self):
        if self.__output is None:
            raise RuntimeError("Survey hasn't been made yet, run method `make_survey` first")
        else:
            return self.__output
    
    @property
    def hdim(self):
        return self.input.hdim
    
    @property
    def inputname(self):
        return self.input.name
    
    @property
    def __ebf_output_files_glob(self):
        return self.output._ebf_glob


if __name__ == '__main__':
    pass
