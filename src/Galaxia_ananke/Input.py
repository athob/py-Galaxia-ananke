#!/usr/bin/env python
"""
Contains the Input class definition

Please note that this module is private. The Input class is
available in the main ``Galaxia`` namespace - use that instead.
"""
import re
import pathlib
import numpy as np
import ebf

from .constants import *
from .utils import make_symlink, compare_given_and_required

__all__ = ['Input']


class Input:
    _pos = 'pos3'
    _vel = 'vel3'
    _mass = 'mass'
    _kernels = 'h_cubic'
    _density = 'density'
    _required_keys_in_particles = {_pos, _vel, _mass, 'age', 'parentid', 'feh'}
    _optional_keys_in_particles = {'id', 'dform', 'helium', 'carbon', 'nitrogen', 'oxygen', 'neon', 'magnesium', 'silicon', 'sulphur', 'calcium', 'alpha'}
    def __init__(self, *args, **kwargs) -> None:
        if args:
            if len(args) not in [2,3]: raise  # TODO mix & match args & kwargs for particles and rho_pos
            self.__input_dir = GALAXIA_TMP
            kwargs['particles'] = args[0]
            kwargs['rho_pos'] = args[1]
            kwargs['rho_vel'] = args[2] if len(args) == 3 else kwargs.get('rho_vel', None)
            self.__input_files_exist = False
        elif {'pname', 'kname'}.issubset(kwargs.keys()):  # TODO implement case where pname and kname non-formated names
            _pname = kwargs['pname'] = pathlib.Path(kwargs['pname'])
            _kname = kwargs['kname'] = pathlib.Path(kwargs['kname'])
            if _pname.parent != _kname.parent: raise  # TODO write error message
            self.__input_dir = _pname.parent
            kwargs['name'] = re.findall("(.*).ebf", _pname.name)[0]
            _hdim, kwargs['ngb'] = map(int, re.findall(f"{kwargs['name']}_d(\d*)n(\d*)_den.ebf",
                                                       _kname.name)[0])  # TODO what if _hdim is 3 ?
            kwargs['particles'] = ebf.read(_pname)
            _k =  ebf.read(_kname)
            _mass = _k[self._mass]  # dummy line to check format
            kwargs['rho_pos'] = _k[self._density]
            _knorm = _k['h_cubic'][:,0] * np.cbrt(_k['density']) / np.sqrt(kwargs['ngb'])
            kwargs['knorm'] = np.median(_knorm) if len(np.unique(np.round(_knorm/(2*np.finfo(_knorm.dtype).eps)).astype('int')))==1 else _knorm
            kwargs['rho_vel'] = (np.sqrt(kwargs['ngb']) * kwargs['knorm'] / _k[self._kernels][:,1])**3
            self.__input_files_exist = True
        else:
            raise  # TODO write error message
        self.__particles = kwargs['particles']
        self.__verify_particles()
        self.__rho_pos = kwargs['rho_pos']
        self.__rho_vel = kwargs.get('rho_vel')
        self.__name = kwargs.get('name', 'sim')
        self.__pname = kwargs.get('pname', None)
        self.__kname = kwargs.get('kname', None)
        self.__knorm = kwargs.get('knorm', 0.596831)
        self.__ngb = kwargs.get('ngb', 64)

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
    def _input_dir(self):
        return self.__input_dir

    @property
    def kernels(self):
        return np.sqrt(self.ngb) * self.knorm / np.cbrt(self.rho)
    
    @property
    def kname(self):
        if self.__kname is None:
            self.__kname = self._input_dir / f"{self.name}_d{self.hdim}n{self.ngb}_den.ebf"
        return self.__kname

    @property
    def pname(self):
        if self.__pname is None:
            self.__pname = self._input_dir / f"{self.name}.ebf"
        return self.__pname

    def keys(self):
        return self.particles.keys()
    
    def optional_keys(self):
        return list(set(self.keys()).intersection(self._optional_keys_in_particles))

    def __verify_particles(self):
        compare_given_and_required(self.keys(), self._required_keys_in_particles, self._optional_keys_in_particles,
                                   error_message="Given particle data covers wrong set of keys")
        if 'dform' not in self.particles:
            self.particles['dform'] = 0*self.particles[self._mass]
        # TODO check format, if dataframe-like

    def prepare_input(self, isochrone, cmd_magnames, **kwargs):
        cmd_magnames = isochrone.check_cmd_magnames(cmd_magnames)
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
            parfile = self._input_dir / parfile
        for_parfile = DEFAULTS_FOR_PARFILE.copy()
        for_parfile.update(**{TTAGS.photo_categ: isochrone.category, TTAGS.photo_sys: isochrone.name, TTAGS.mag_color_names: cmd_magnames, TTAGS.nres: self.ngb}, **kwargs)
        parfile.write_text(PARFILE_TEMPLATE.substitute(for_parfile))
        return self.name, parfile, for_parfile

    def _write_kernels(self):
        kname = self.kname
        if not self.__input_files_exist:
            ebf.initialize(self.kname)
            ebf.write(kname, f"/{self._density}", self.rho_pos, "a")
            ebf.write(kname, f"/{self._kernels}", self.kernels, "a")
            ebf.write(kname, f"/{self._mass}", self.particles['mass'], "a")
        return kname

    def _write_particles(self):
        pname = self.pname
        if not self.__input_files_exist:
            ebf.initialize(pname)
            for key in self._required_keys_in_particles:
                ebf.write(pname, f"/{key}", self.particles[key], 'a')
            for key in self._optional_keys_in_particles:
                ebf.write(pname, f"/{key}", self.particles[key] if key in self.keys() else np.zeros(self.length), 'a')
        return pname


if __name__ == '__main__':
    pass
