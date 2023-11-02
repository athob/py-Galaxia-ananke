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
    _position_prop = ('pos3', "Position coordinates in kpc (Nx3)")
    _velocity_prop = ('vel3', "Velocity coordinates in km/s (Nx3)")
    _mass_prop = ('mass', "Stellar masses in solar masses")
    _age_prop = ('age', "Stellar ages in years and decimal logarithmic scale")
    _metallicity_prop = ('feh', "Stellar metallicity [Fe/H] in dex relative to solar")
    _parentindex_prop = ('parentid', "Index of parent particle")
    _required_properties = {_position_prop, _velocity_prop, _mass_prop, _age_prop, _metallicity_prop, _parentindex_prop}
    _populationindex_prop = ('id', "Index of parent particle population")
    _formationdistance_prop = ('dform', "Formation distance of parent particle in kpc")
    _heliumabundance_prop = ('helium', "Helium abundance [He/H] in dex")
    _carbonabundance_prop = ('carbon', "Carbon abundance [C/H] in dex")
    _nitrogenabundance_prop = ('nitrogen', "Nitrogen abundance [N/H] in dex")
    _oxygenabundance_prop = ('oxygen', "Oxygen abundance [O/H] in dex")
    _neonabundance_prop = ('neon', "Neon abundance [Ne/H] in dex")
    _magnesiumabundance_prop = ('magnesium', "Magnesium abundance [Mg/H] in dex")
    _siliconabundance_prop = ('silicon', "Silicon abundance [Si/H] in dex")
    _sulphurabundance_prop = ('sulphur', "Sulphur abundance [S/H] in dex")
    _calciumabundance_prop = ('calcium', "Calcium abundance [Ca/H] in dex")
    _alphaabundance_prop = ('alpha', "Alpha abundance [Mg/Fe] in dex")
    _optional_properties =  {_populationindex_prop, _formationdistance_prop, _heliumabundance_prop, _carbonabundance_prop, _nitrogenabundance_prop, _oxygenabundance_prop, _neonabundance_prop, _magnesiumabundance_prop, _siliconabundance_prop, _sulphurabundance_prop, _calciumabundance_prop, _alphaabundance_prop}
    _required_keys_in_particles = {_property[0] for _property in _required_properties}
    _optional_keys_in_particles = {_property[0] for _property in _optional_properties}
    _pos = _position_prop[0]
    _vel = _velocity_prop[0]
    _mass = _mass_prop[0]
    _age = _age_prop[0]
    _feh = _metallicity_prop[0]
    _kernels = 'h_cubic'
    _density = 'density'
    def __init__(self, *args, **kwargs) -> None:
        """
            Driver to store and prepare the input data for Galaxia.

            Call signatures::

                input = Input(particles, rho_pos, rho_vel=None, name='{DEFAULT_SIMNAME}',
                ngb={TTAGS_nres}, knorm=0.596831)

                input = Input(pname, kname)
            
            Parameters
            ----------
            particles : dict
                Dictionary where each elements represent the properties of the
                input particles, given as equal-length array_like objects. The
                dictionary must include the following properties with
                corresponding keys: {_required_properties}
                Additionally, Galaxia can optionally receive particle
                properties that will be carried over to the generated
                synthetic star, those include the following: {_optional_properties}

            rho_pos : array_like
                Contains the position-determined kernel density estimates for
                the input particles. Must have equal lengths as the elements
                in the particles dictionary.

            rho_vel : array_like
                Contains the velocity-determined kernel density estimates for
                the input particles. Must have equal lengths as the elements
                in the particles dictionary.

            simname : string
                Optional name Galaxia should use for the input files.
                Default to '{DEFAULT_SIMNAME}'.
            
            ngb : int
                Number of neighbouring particles Galaxia should consider.
                Default to {TTAGS_nres}.
            
            knorm : float
                Kernel normalization factor. Default to 0.596831.
            
            pname : string
                Path to existing pre-formatted particles EBF files to use as
                input for Galaxia. This keyword argument must be used in
                conjunction with kname. Default to None if unused.

            kname : string
                Path to existing pre-formatted kernel EBF files to use as
                input for Galaxia. This keyword argument must be used in
                conjunction with pname. Default to None if unused.
        """
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
        self.__ngb = kwargs.get('ngb', TTAGS.nres)

    __init__.__doc__ = __init__.__doc__.format(DEFAULT_SIMNAME=DEFAULT_SIMNAME,
                                               TTAGS_nres=TTAGS.nres,
                                               _required_properties=''.join(
                                                   [f"\n                 -{desc} via key `{str(key)}`"
                                                     for key, desc in _required_properties]),
                                               _optional_properties=''.join(
                                                   [f"\n                 -{desc} via key `{str(key)}`"
                                                     for key, desc in _optional_properties]))

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
        parfile = pathlib.Path(kwargs.pop('parfile', DEFAULT_PARFILE))  # TODO make temporary? create a global record of temporary files?
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
