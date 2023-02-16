#!/usr/bin/env python
"""
Docstring
"""
import pathlib
import numpy as np
import ebf

from .constants import *
from .utils import *

__all__ = ['Input']


class Input:
    _pos = 'pos3'
    _vel = 'vel3'
    _mass = 'mass'
    _kernels = 'h_cubic'
    _density = 'density'
    _required_keys_in_particles = {_pos, _vel, _mass, 'age', 'parentid', 'dform', 'feh', 'alpha', 'helium', 'carbon', 'nitrogen', 'oxygen', 'neon', 'magnesium', 'silicon', 'sulphur', 'calcium'}
    _optional_keys_in_particles = {'id'}
    def __init__(self, particles, rho_pos, rho_vel=None, name='sim', knorm=0.596831, ngb=64) -> None:
        self.__particles = particles
        self.__verify_particles()
        self.__rho_pos = rho_pos
        self.__rho_vel = rho_vel
        self.__name = name
        self.__pname = None
        self.__kname = None
        self.__knorm = knorm
        self.__ngb = ngb

    @property
    def particles(self):
        return self.__particles
    
    @property
    def length(self):
        return len(self.particles[self._mass])
    
    @property
    def rho_pos(self):
        return self.__rho_pos
    
    @property
    def rho_vel(self):
        return self.__rho_vel

    @property
    def rho(self):
        return np.vstack([self.rho_pos, self.rho_vel]).T  # TODO what if rho_vel is None
    
    @property
    def hdim(self):
        return 3 if self.rho_vel is None else 6

    @property
    def name(self):
        return self.__name
    
    @property
    def knorm(self):
        return self.__knorm

    @property
    def ngb(self):
        return self.__ngb

    @property
    def kernels(self):
        return np.sqrt(self.ngb) * self.knorm / np.cbrt(self.rho)
    
    @property
    def kname(self):
        if self.__kname is None:
            self.__kname = GALAXIA_TMP / f"{self.name}_d{self.hdim}n{self.ngb}_den.ebf"
        return self.__kname

    @property
    def pname(self):
        if self.__pname is None:
            self.__pname = GALAXIA_TMP / f"{self.name}.ebf"
        return self.__pname

    def __verify_particles(self):
        particles_keys = set(self.particles.keys())
        if particles_keys != self._required_keys_in_particles:
            missing = self._required_keys_in_particles.difference(particles_keys)
            missing = f"misses {missing}" if missing else ""
            extra = particles_keys.difference(self._required_keys_in_particles.union(self._optional_keys_in_particles))
            extra = f"misincludes {extra}" if extra else ""
            raise ValueError(f"Given particle data covers wrong set of keys: {missing}{' & ' if missing and extra else ''}{extra}")
        # TODO check format, if dataframe-like

    def prepare_input(self, isochrone, cmd_magnames, **kwargs):
        isochrone.check_cmd_magnames(cmd_magnames)
        kname = self._write_kernels()
        pname = self._write_particles()
        temp_dir = GALAXIA_NBODY1 / self.name
        temp_dir.mkdir(parents=True, exist_ok=True)
        make_symlink(kname, temp_dir)
        make_symlink(pname, temp_dir)
        temp_filename = (GALAXIA_FILENAMES / self.name).with_suffix('.txt')
        temp_filename.write_text(FILENAME_TEMPLATE.substitute(name=self.name, pname=pname.name))
        parfile = pathlib.Path(kwargs.pop('parfile', 'survey_params'))  # TODO make temporary? create a global record of temporary files?
        if not parfile.is_absolute():
            parfile = GALAXIA_TMP / parfile
        parfile.write_text(PARFILE_TEMPLATE.substitute(DEFAULTS_FOR_PARFILE, photo_categ=isochrone.category, photo_sys=isochrone.name, mag_color_names=cmd_magnames, nres=self.ngb, **kwargs))  # output_file=surveyname, fsample=fsample))
        return self.name, parfile

    def _write_kernels(self):
        kname = self.kname
        ebf.initialize(self.kname)
        ebf.write(kname, f"/{self._density}", self.rho_pos, "a")
        ebf.write(kname, f"/{self._kernels}", self.kernels, "a")
        ebf.write(kname, f"/{self._mass}", self.particles['mass'], "a")
        return kname

    def _write_particles(self):
        pname = self.pname
        ebf.initialize(pname)
        for key in self._required_keys_in_particles:
            ebf.write(pname, f"/{key}", self.particles[key], 'a')
        for key in self._optional_keys_in_particles:
            ebf.write(pname, f"/{key}", self.particles[key] if key in self.particles.keys() else np.zeros(self.length), 'a')
        return pname


if __name__ == '__main__':
    pass
