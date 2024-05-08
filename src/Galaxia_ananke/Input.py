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
from astropy.utils import classproperty

from .constants import *
from .templates import *
from .defaults import *
from .utils import make_symlink, compare_given_and_required, confirm_equal_length_arrays_in_dict
from .photometry.PhotoSystem import PhotoSystem

__all__ = ['Input']


FOURTHIRDPI = 4*np.pi/3

class Input:
    _position_prop = ('pos3', "Position coordinates in kpc (Nx3)")
    _velocity_prop = ('vel3', "Velocity coordinates in km/s (Nx3)")
    _mass_prop = ('mass', "Stellar masses in solar masses")
    _age_prop = ('age', "Stellar ages in years and decimal logarithmic scale")
    _metallicity_prop = ('feh', "Stellar metallicity [Fe/H] in dex relative to solar")
    _parentindex_prop = ('parentid', "Index of parent particle")
    _populationindex_prop = ('id', "Index of parent particle population")
    _partitionindex_prop = ('partitionid', "Index of the data partition that contains the particle")
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
    _kernels = 'h_cubic'
    _density = 'density'
    _positiondensity_prop = ('rho_pos', 'Position space density in kpc^-3')
    _velocitydensity_prop = ('rho_vel', 'Velocity space density in [km/s]^-3')
    def __init__(self, *args, **kwargs) -> None:
        """
            Driver to store and prepare the input data for Galaxia.

            Call signatures::

                input = Input(particles, {rho_pos}, {rho_vel}=None,
                input_dir='{GALAXIA_TMP}', name='{DEFAULT_SIMNAME}',
                ngb={TTAGS_nres}, k_factor=1., former_kernel=False)

                input = Input(pname, kname, former_kernel=False)
            
            Parameters
            ----------
            particles : dict
                Dictionary where each elements represent the properties of the
                input particles, given as equal-length array_like objects.
                {particles_dictionary_description}
                Input comes with class method make_dummy_particles_input to
                produce a dummy example of that dictionary.

            {rho_pos} : array_like
                Contains the position-determined kernel density estimates for
                the input particles. Must have equal lengths as the elements
                in the particles dictionary.

            {rho_vel} : array_like
                Contains the velocity-determined kernel density estimates for
                the input particles. Must have equal lengths as the elements
                in the particles dictionary.

            input_dir : string or pathlib.Path
                Optional arguments to specify path for the directory where
                Galaxia's input data should be generated. Default to '{GALAXIA_TMP}'.

            name : string
                Optional name Galaxia should use for the input files.
                Default to '{DEFAULT_SIMNAME}'.
            
            ngb : int
                Number of neighbouring particles Galaxia should consider.
                Default to {TTAGS_nres}. ONLY SUPPORT 8, 32, 64 & 128.

            k_factor : float
                Scaling factor applied to the kernels lengths to adjust all
                the kernels sizes uniformly. Lower values reduces the kernels
                extents, while higher values increases them.
                Default to 1 (no adjustment).
            
            pname : string
                Path to existing pre-formatted particles EBF files to use as
                input for Galaxia. This keyword argument must be used in
                conjunction with kname. Default to None if unused.

            kname : string
                Path to existing pre-formatted kernel EBF files to use as
                input for Galaxia. This keyword argument must be used in
                conjunction with pname. Default to None if unused.

            former_kernel : bool or dict
                Flag that allow the utilization of the former implementation
                for the kernels lengths with the consideration of a kernel
                normalization factor. If providing a dictionary, you can
                configure the former kernel normalization factor knorm by
                including the value knorm under key 'knorm'. Default to False,
                if True, knorm defaults to 0.596831.
        """
        if args:
            if len(args) not in [2,3]: raise  # TODO mix & match args & kwargs for particles and rho_pos
            kwargs['particles'] = args[0]
            kwargs[self._rho_pos] = args[1]
            kwargs[self._rho_vel] = args[2] if len(args) == 3 else kwargs.get(self._rho_vel, None)
            self.__input_files_exist = False
        elif {'pname', 'kname'}.issubset(kwargs.keys()):  # TODO implement case where pname and kname non-formated names
            _pname = kwargs['pname'] = pathlib.Path(kwargs['pname'])
            _kname = kwargs['kname'] = pathlib.Path(kwargs['kname'])
            if _pname.parent != _kname.parent: raise ValueError(f"Given pname file {_pname} and kname file {_kname} need to share the same parent directory")
            kwargs['input_dir'] = _pname.parent
            kwargs['name'] = re.findall("(.*).ebf", _pname.name)[0]
            _hdim, kwargs['ngb'] = map(int, re.findall(f"{kwargs['name']}_d(\\d*)n(\\d*)_den.ebf",
                                                       _kname.name)[0])  # TODO what if _hdim is 3 ?
            kwargs['particles'] = ebf.read(_pname)
            _k =  ebf.read(_kname)
            _mass = _k[self._mass]  # dummy line to check format
            kwargs[self._rho_pos] = _k[self._density]
            _k_factor = _k[self._kernels][:,0] * np.cbrt(FOURTHIRDPI*_k[self._density])
            if kwargs.get('former_kernel', False):
                _knorm = _k_factor/(np.sqrt(kwargs['ngb']) * np.cbrt(FOURTHIRDPI))
                # _knorm = _k[self._kernels][:,0] * np.cbrt(_k[self._density]) / np.sqrt(kwargs['ngb'])
                _knorm = np.median(_knorm) if len(np.unique(np.round(_knorm/(2*np.finfo(_knorm.dtype).eps)).astype('int')))==1 else _knorm
                kwargs[self._rho_vel] = (np.sqrt(kwargs['ngb']) * _knorm / _k[self._kernels][:,1])**3
                kwargs['former_kernel'] = {'knorm': _knorm}
            else:
                kwargs[self._rho_vel] = (_k_factor / _k[self._kernels][:,1])**3 / FOURTHIRDPI
                kwargs['k_factor'] = _k_factor
            self.__input_files_exist = True
        else:
            raise ValueError("Wrong signature: please consult help of the Input constructor")
        self.__particles = kwargs['particles'].copy()
        self.__verify_particles(self.particles)
        self.__complete_particles(self.particles)
        self.__pos_density = kwargs[self._rho_pos]
        self.__vel_density = kwargs.get(self._rho_vel)
        self.__input_dir = pathlib.Path(kwargs.get('input_dir', GALAXIA_TMP))
        self.__name = kwargs.get('name', DEFAULT_SIMNAME)
        self.__pname = kwargs.get('pname', None)
        self.__kname = kwargs.get('kname', None)
        self.__ngb = kwargs.get('ngb', TTAGS.nres)
        __old = kwargs.get('former_kernel', False)
        if __old and not isinstance(__old, dict):  __old = {}
        __knorm = __old.get('knorm', 0.596831) if isinstance(__old, dict) else None
        self.__k_factor = kwargs.get('k_factor', 1. if __knorm is None else np.sqrt(self.ngb) * __knorm * np.cbrt(FOURTHIRDPI))

    @classproperty
    def particles_dictionary_description(cls):
        description = """
                The dictionary includes the following properties with
                corresponding keys:
                {_required_properties}
                
                Additionally, Galaxia can optionally receive particle
                properties that will be carried over to the generated
                synthetic star, those include the following: 
                {_optional_properties}
        """.format(_required_properties=''.join(
                       [f"\n                 -{desc} via key `{str(key)}`"
                        for key, desc in Input._required_properties]),
                   _optional_properties=''.join(
                       [f"\n                 -{desc} via key `{str(key)}`"
                        for key, desc in Input._optional_properties]))
        return description

    @classproperty
    def _required_properties(cls):
        return {
            cls._position_prop,
            cls._velocity_prop,
            cls._mass_prop,
            cls._age_prop,
            cls._metallicity_prop
            }

    @classproperty
    def _optional_properties(cls):
        return {
            cls._parentindex_prop,
            cls._partitionindex_prop,
            cls._populationindex_prop,
            cls._formationdistance_prop,
            cls._heliumabundance_prop,
            cls._carbonabundance_prop,
            cls._nitrogenabundance_prop,
            cls._oxygenabundance_prop,
            cls._neonabundance_prop,
            cls._magnesiumabundance_prop,
            cls._siliconabundance_prop,
            cls._sulphurabundance_prop,
            cls._calciumabundance_prop,
            cls._alphaabundance_prop
            }

    @classproperty
    def _required_keys_in_particles(cls):
        return {_property[0] for _property in cls._required_properties}

    @classproperty
    def _optional_keys_in_particles(cls):
        return {_property[0] for _property in cls._optional_properties}

    @classproperty
    def all_possible_keys_in_particles(cls):
        return cls._required_keys_in_particles.union(cls._optional_keys_in_particles)
    
    @classproperty
    def _pos(cls):
        return cls._position_prop[0]

    @classproperty
    def _vel(cls):
        return cls._velocity_prop[0]

    @classproperty
    def _mass(cls):
        return cls._mass_prop[0]

    @classproperty
    def _age(cls):
        return cls._age_prop[0]

    @classproperty
    def _feh(cls):
        return cls._metallicity_prop[0]

    @classproperty
    def _He(cls):
        return cls._heliumabundance_prop[0]

    @classproperty
    def _C(cls):
        return cls._carbonabundance_prop[0]

    @classproperty
    def _N(cls):
        return cls._nitrogenabundance_prop[0]

    @classproperty
    def _O(cls):
        return cls._oxygenabundance_prop[0]

    @classproperty
    def _Ne(cls):
        return cls._neonabundance_prop[0]

    @classproperty
    def _Mg(cls):
        return cls._magnesiumabundance_prop[0]

    @classproperty
    def _Si(cls):
        return cls._siliconabundance_prop[0]

    @classproperty
    def _S(cls):
        return cls._sulphurabundance_prop[0]

    @classproperty
    def _Ca(cls):
        return cls._calciumabundance_prop[0]

    @classproperty
    def _elem_list(cls):
        return [cls._He, cls._C, cls._N, cls._O, cls._Ne, cls._Mg, cls._Si, cls._S, cls._Ca]

    @classproperty
    def _alph(cls):
        return cls._alphaabundance_prop[0]

    @classproperty
    def _parentid(cls):
        return cls._parentindex_prop[0]

    @classproperty
    def _partitionid(cls):
        return cls._partitionindex_prop[0]

    @classproperty
    def _dform(cls):
        return cls._formationdistance_prop[0]

    @classproperty
    def _pop_id(cls):
        return cls._populationindex_prop[0]

    @classproperty
    def _rho_pos(cls):
        return cls._positiondensity_prop[0]

    @classproperty
    def _rho_vel(cls):
        return cls._velocitydensity_prop[0]

    @property
    def particles(self):
        return self.__particles
    
    @property
    def length(self):
        return len(self.particles[self._mass])
    
    @property
    def rho_pos(self):
        return self.__pos_density
    
    @property
    def rho_vel(self):
        return self.__vel_density

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
    def ngb(self):
        return self.__ngb
    
    @property
    def k_factor(self):
        return self.__k_factor

    @property
    def _input_dir(self):
        return self.__input_dir

    @property
    def kernels(self):
        return self.k_factor/np.cbrt(FOURTHIRDPI*self.rho)
    
    @property
    def kname(self):  # TODO replace with functools cached_property?
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

    def prepare_input(self, photosys: PhotoSystem, cmd_magnames, **kwargs):
        cmd_magnames = photosys.check_cmd_magnames(cmd_magnames)
        parfile, for_parfile = self._write_parameter_file(photosys, cmd_magnames, **kwargs)
        sorter = np.lexsort((self.particles[self._partitionid],))
        kname = self._write_kernels(sorter)
        pname = self._write_particles(sorter)
        temp_filename = self._prepare_nbody1(kname, pname)
        return self.name, parfile, for_parfile

    def _write_parameter_file(self, photosys: PhotoSystem, cmd_magnames, **kwargs):
        parfile = pathlib.Path(kwargs.pop('parfile', DEFAULT_PARFILE))  # TODO make temporary? create a global record of temporary files?
        if not parfile.is_absolute():
            parfile = self._input_dir / parfile
        for_parfile = DEFAULTS_FOR_PARFILE.copy()
        for_parfile.update(**{TTAGS.photo_categ: photosys.category, TTAGS.photo_sys: photosys.name, TTAGS.mag_color_names: cmd_magnames, TTAGS.nres: self.ngb}, **kwargs)
        parfile.write_text(PARFILE_TEMPLATE.substitute(for_parfile))
        return parfile, for_parfile

    def _write_kernels(self, sorter):
        kname = self.kname
        if not self.__input_files_exist:
            ebf.initialize(self.kname)
            ebf.write(kname, f"/{self._density}", self.rho_pos[sorter], "a")
            ebf.write(kname, f"/{self._kernels}", self.kernels[sorter], "a")
            ebf.write(kname, f"/{self._mass}", self.particles[self._mass][sorter], "a")
        return kname

    def _write_particles(self, sorter):
        pname = self.pname
        if not self.__input_files_exist:
            ebf.initialize(pname)
            for key in self._required_keys_in_particles:
                ebf.write(pname, f"/{key}", self.particles[key][sorter], 'a')
            for key in self._optional_keys_in_particles:
                ebf.write(pname, f"/{key}", self.particles[key][sorter] if key in self.keys() else np.zeros(self.length), 'a')
        return pname
    
    def _prepare_nbody1(self, kname: pathlib.Path, pname: pathlib.Path):
        temp_dir = GALAXIA_NBODY1 / self.name
        temp_dir.mkdir(parents=True, exist_ok=True)
        make_symlink(kname, temp_dir)
        make_symlink(pname, temp_dir)
        temp_filename = (GALAXIA_FILENAMES / self.name).with_suffix('.txt')
        temp_filename.write_text(FILENAME_TEMPLATE.substitute(name=self.name, pname=pname.name))
        return temp_filename
    
    @classmethod
    def __verify_particles(cls, particles):
        compare_given_and_required(particles.keys(), cls._required_keys_in_particles, cls._optional_keys_in_particles,
                                   error_message="Given particle data covers wrong set of keys")
        confirm_equal_length_arrays_in_dict(particles, cls._mass, error_message_dict_name='particles')
        # TODO check format, if dataframe-like
    
    @classmethod
    def __complete_particles(cls, particles):
        if cls._parentid not in particles:
            particles[cls._parentid] = np.arange(particles[cls._mass].shape[0])
        if cls._partitionid not in particles:
            particles[cls._partitionid] = np.zeros(particles[cls._mass].shape[0], dtype='int')
        if cls._dform not in particles:
            particles[cls._dform] = 0*particles[cls._mass]

    @classmethod
    def make_dummy_particles_input(cls, n_parts=10**5):
        """
            Generate an example dummy input particles dictionary for Input
            made of randomly generated arrays.

            Parameters
            ----------
            n_parts : int
                Number of particles the example include. Default to 10**5.

            Returns
            -------
            p : dict
                Dummy example input particles dictionary for Input.
                {particles_dictionary_description}
        """
        p = {}
        p[cls._pos] = 30*np.random.randn(n_parts, 3)
        p[cls._vel] = 50*np.random.randn(n_parts, 3)
        p[cls._mass] = 5500 + 700*np.random.randn(n_parts)
        p[cls._age] = 9.7 + 0.4*np.random.randn(n_parts)
        p[cls._feh] = -0.7 + 0.4*np.random.randn(n_parts)
        for el in cls._elem_list:
            p[el] = -0.6 + 0.4*np.random.randn(n_parts)
        p[cls._alph] = p[cls._Mg] - p[cls._feh]
        p[cls._parentid] = np.arange(n_parts).astype('int')
        p[cls._partitionid] = np.zeros(n_parts, dtype='int')
        p[cls._dform] = np.zeros(n_parts, dtype='float32')
        p[cls._pop_id] = np.zeros(n_parts, dtype='int')
        return p

    @classmethod
    def make_dummy_densities_input(cls, n_parts=10**5):
        """
            Generate an example dummy input densities for Input
            made of randomly generated arrays.

            Parameters
            ----------
            n_parts : int
                Number of particles the example include. Default to 10**5.

            Returns
            -------
            rho_pos : array
                Dummy example input position densities for Input.

            rho_vel : array
                Dummy example input velocity densities for Input.
        """
        rho_pos = np.exp(-2.9 + 1.1*np.random.randn(n_parts))
        rho_vel = np.exp(-4.4 + 1.1*np.random.randn(n_parts))
        return rho_pos, rho_vel


Input.__init__.__doc__ = Input.__init__.__doc__.format(GALAXIA_TMP=GALAXIA_TMP,
                                                       DEFAULT_SIMNAME=DEFAULT_SIMNAME,
                                                       TTAGS_nres=TTAGS.nres,
                                                       particles_dictionary_description=Input.particles_dictionary_description,
                                                       rho_pos=Input._rho_pos,
                                                       rho_vel=Input._rho_vel)

Input.make_dummy_particles_input.__func__.__doc__ = Input.make_dummy_particles_input.__doc__.format(
    particles_dictionary_description=Input.particles_dictionary_description)


if __name__ == '__main__':
    pass
