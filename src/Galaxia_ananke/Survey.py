#!/usr/bin/env python
"""
Contains the Survey class definition

Please note that this module is private. The Survey class is
available in the main ``Galaxia`` namespace - use that instead.
"""
import subprocess

from .constants import *
from . import photometry
from .Output import Output

__all__ = ['Survey']


class Survey:
    def __init__(self, input, photo_sys=DEFAULT_PSYS, surveyname='survey') -> None:
        self.__surveyname = surveyname
        self.__input = input
        self._set_isochrones_from_photosys(photo_sys)
        self.__output = None

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['surveyname', 'photo_sys']])
        return f'{cls}({description})'

    def _set_isochrones_from_photosys(self, photo_sys):
        if isinstance(photo_sys, str):
            photo_sys = [photo_sys]
        self.__isochrones = [photometry.available_photo_systems[psys] for psys in photo_sys]

    def _run_survey(self, cmd_magnames, fsample, **kwargs):
        inputname, parfile = self.input.prepare_input(self.isochrones[0], cmd_magnames, output_file=self.surveyname, fsample=fsample, **kwargs)
        cmd = f"{GALAXIA} -r{(' --hdim=' + str(self.hdim) if self.hdim is not None else '')} --nfile={self.inputname} {parfile}"
        print(cmd)
        subprocess.call(cmd.split(' '))
        self.__output = Output(self)

    def _append_survey(self, isochrone):
        cmd = f"{GALAXIA} -a --pcat={isochrone.category} --psys={isochrone.name} {self.__ebf_output_file}"
        print(cmd)
        subprocess.call(cmd.split(' '))

    def make_survey(self, cmd_magnames=DEFAULT_CMD, fsample=1, **kwargs):
        self._run_survey(cmd_magnames, fsample, **kwargs)
        for isochrone in self.isochrones[1:]:
            self._append_survey(isochrone)
        self.output._ebf_to_hdf5()
        return self.output

    @property
    def surveyname(self):
        return self.__surveyname
    
    @property
    def input(self):
        return self.__input

    @property
    def isochrones(self):
        return self.__isochrones

    @property
    def photo_sys(self):
        return {isochrone.key for isochrone in self.__isochrones}

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
    def __ebf_output_file(self):
        return self.output._ebf


if __name__ == '__main__':
    pass
