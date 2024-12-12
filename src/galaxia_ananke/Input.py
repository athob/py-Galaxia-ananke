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
Contains the Input class definition

Please note that this module is private. The Input class is
available in the main ``galaxia_ananke`` namespace - use that instead.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Any, Optional, Union, Tuple, List, Dict, OrderedDict
from numpy.typing import NDArray, ArrayLike
from warnings import warn
from functools import cached_property
import itertools
import re
import pathlib
import numpy as np
import ebf
# from astropy.utils import classproperty

from ._constants import *
from ._templates import *
from ._defaults import *
from .utils import classproperty, make_symlink, compare_given_and_required, confirm_equal_length_arrays_in_dict, lexicalorder_dict, mark_metadata_prop, collect_metadata_marked_properties, hash_iterable
from .photometry.PhotoSystem import PhotoSystem

if TYPE_CHECKING:
    from . import Survey

__all__ = ['Input']


@collect_metadata_marked_properties
class Input:
    _position_prop = ('pos3', "Position coordinates in $kpc$ (Nx3)")
    _velocity_prop = ('vel3', "Velocity coordinates in $km/s$ (Nx3)")
    _masscurrent_prop = ('mass', "Present-day stellar masses in solar masses")
    _massinitial_prop = ('massinit', "Initial stellar masses in solar masses")
    _age_prop = ('age', "Stellar ages in years and decimal logarithmic scale")
    _metallicity_prop = ('feh', "Stellar metallicity $[Fe/H]$ in dex relative to solar")
    _parentindex_prop = ('parentid', "Index of parent particle")
    _populationindex_prop = ('id', "Index of parent particle population")
    _partitionindex_prop = ('partitionid', "Index of the data partition that contains the particle")
    _formationdistance_prop = ('dform', "Formation distance of parent particle in kpc")
    _heliumabundance_prop = ('helium', "Helium abundance $[He/H]$ in $dex$")
    _carbonabundance_prop = ('carbon', "Carbon abundance $[C/H]$ in $dex$")
    _nitrogenabundance_prop = ('nitrogen', "Nitrogen abundance $[N/H]$ in $dex$")
    _oxygenabundance_prop = ('oxygen', "Oxygen abundance $[O/H]$ in $dex$")
    _neonabundance_prop = ('neon', "Neon abundance $[Ne/H]$ in $dex$")
    _magnesiumabundance_prop = ('magnesium', "Magnesium abundance $[Mg/H]$ in $dex$")
    _siliconabundance_prop = ('silicon', "Silicon abundance $[Si/H]$ in $dex$")
    _sulphurabundance_prop = ('sulphur', "Sulphur abundance $[S/H]$ in $dex$")
    _calciumabundance_prop = ('calcium', "Calcium abundance $[Ca/H]$ in $dex$")
    _alphaabundance_prop = ('alpha', "Alpha abundance $[Mg/Fe]$ in $dex$")
    _kernels_prop = ('h_cubic', "Phase-space kernel radii in $kpc$ and $km/s$ (Nx2)")
    def __init__(self, *args, **kwargs) -> None:
        """
            Driver to store and prepare the input data for Galaxia.

            Call signatures::

                input = Input(particles, kernels,
                              input_dir='{GALAXIA_TMP}',
                              name='{DEFAULT_SIMNAME}', caching=False,
                              ngb={TTAGS_nres},
                              k_factor=1.)

                input = Input(pname, kname, caching=False)
            
            Parameters
            ----------
            particles : dict
                Dictionary where each elements represent the properties of the
                input particles, given as equal-length array_like objects.
                {particles_dictionary_description}
                Input comes with class method make_dummy_particles_input to
                produce a dummy example of that dictionary.

            kernels : array_like
                Array containing the kernel characteristic lengths for each
                input particle. Must be of equal length N as the arrays
                provided in the particles dictionary, as either a (N) or a
                (Nx2) array if representing respectively kernels in position
                space, or in phase space (using the same position and velocity
                unit as that of the corresponding particles dictionary entry).

            input_dir : string or pathlib.Path
                Optional arguments to specify path for the directory where
                Galaxia's input data should be generated. Default to '{GALAXIA_TMP}'.

            name : string
                Optional name Galaxia should use for the input files.
                Default to '{DEFAULT_SIMNAME}'.

            caching : bool
                TODO

            append_hash : bool
                TODO

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
        """
        self.caching: bool = kwargs.get('caching', False)
        self.__append_hash: bool = kwargs.get('append_hash', self.caching)
        # check which args/kwargs signature case and populate kwargs accordingly
        if args:
            if len(args) not in [2]: raise  # TODO mix & match args & kwargs for particles and kernels
            kwargs['particles'] = args[0]
            kwargs['kernels'] = args[1]
            self.__input_files_exist: bool = False
        elif {'pname', 'kname'}.issubset(kwargs.keys()):  # TODO implement case where pname and kname non-formated names
            _pname = kwargs['pname'] = pathlib.Path(kwargs['pname'])
            _kname = kwargs['kname'] = pathlib.Path(kwargs['kname'])
            if _pname.parent != _kname.parent: raise ValueError(f"Given pname file {_pname} and kname file {_kname} need to share the same parent directory")
            kwargs['input_dir'] = _pname.parent
            kwargs['name'] = re.findall("(.*).ebf", _pname.name)[0]
            _hdim, kwargs['ngb'] = map(int, re.findall(f"{kwargs['name']}_d(\\d*)n(\\d*)_den.ebf",
                                                       _kname.name)[0])  # TODO what if _hdim is 3 ?
            kwargs['particles'] = ebf.read(_pname)
            _k: Dict[str, NDArray] =  ebf.read(_kname)
            kwargs['kernels'] = _k[self._kernels]
            kwargs['k_factor'] = 1.
            self.__input_files_exist: bool = True
        else:
            raise ValueError("Wrong signature: please consult help of the Input constructor")
        # verify and assign attributes based on kwargs
        self.__particles: OrderedDict[str, NDArray] = lexicalorder_dict(kwargs['particles'].copy())
        self.__verify_particles(self.particles)
        self.__complete_particles(self.particles)
        self.__kernels: NDArray = self.__conform_kernels(kwargs['kernels'])
        self.__input_dir: pathlib.Path = pathlib.Path(kwargs.get('input_dir', GALAXIA_TMP))
        self.__name: str = kwargs.get('name', DEFAULT_SIMNAME)
        self.__pname: Optional[pathlib.Path] = kwargs.get('pname', None)
        self.__kname: Optional[pathlib.Path] = kwargs.get('kname', None)
        self.__ngb: int = kwargs.get('ngb', FTTAGS.nres)  # TODO Galaxia makes no use of it, need to remove in the future
        self.__k_factor: Union[float, NDArray] = kwargs.get('k_factor', 1.)

    @classproperty
    def particles_dictionary_description(cls):
        description = """
            The particle dictionary includes the following properties with
            corresponding keys:
            {_required_properties}
            
            Additionally, Galaxia can optionally receive particle properties
            that will be carried over to the generated synthetic star, those
            include the following: 
            {_optional_properties}
        """.format(_required_properties=''.join(
                       [f"\n            * {desc} via key ``{str(key)}``"
                        for key, desc in Input._required_properties]),
                   _optional_properties=''.join(
                       [f"\n            * {desc} via key ``{str(key)}``"
                        for key, desc in Input._optional_properties]))
        return description

    @classproperty
    def _required_properties(cls):
        return {
            cls._position_prop,
            cls._velocity_prop,
            cls._masscurrent_prop,
            cls._massinitial_prop,
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
        return cls._masscurrent_prop[0]
    
    @classproperty
    def _massinit(cls):
        return cls._massinitial_prop[0]

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
    def _kernels(cls):
        return cls._kernels_prop[0]

    @property
    def caching(self) -> bool:
        return self.__caching

    @caching.setter
    def caching(self, value: bool) -> None:
        if value:
            warn(f"You have requested caching mode, be aware this feature is currently experimental and may result in unintended behaviour.", DeprecationWarning, stacklevel=2)
        self.__caching: bool = value
    
    @property
    def particles(self) -> Dict[str, NDArray]:
        return self.__particles
    
    @property
    def length(self) -> int:
        return len(self.particles[self._mass])
    
    @property
    def hdim(self) -> int:  # TODO how to adapt for ellipsoidal kernels?
        return 3 if self.kernels.ndim == 1 else 6

    @property
    @mark_metadata_prop
    def name(self) -> str:
        return self.__name
    
    @property
    def append_hash(self) -> bool:
        return self.__append_hash

    @property
    def name_hash(self) -> str:
        return self.name + (f"_{self.hash[:7]}" if self.append_hash else "")

    @property
    def ngb(self) -> int:
        return self.__ngb
    
    @property
    def k_factor(self) -> Union[float, NDArray]:
        return self.__k_factor

    @property
    def _input_dir(self) -> pathlib.Path:
        return self.__input_dir

    @property
    def kernels(self) -> NDArray:
        return self.__kernels
    
    @property
    def _base_inputfile(self) -> pathlib.Path:  # TODO what if pname and kname already exist (case where input args are pname and kname)?
        return self._input_dir / self.name_hash

    @property
    def _hashname(self) -> pathlib.Path:
        return self._base_inputfile.with_suffix(f".{HASH_EXT}")

    @property
    def pname(self) -> pathlib.Path:
        if self.__pname is None:
            self.__pname = self._base_inputfile.with_suffix(".ebf")
        return self.__pname

    @property
    def kname(self) -> pathlib.Path:
        if self.__kname is None:
            self.__kname = self.pname.with_name(f"{self.pname.stem}_d{self.hdim}n{self.ngb}_den.ebf")
        return self.__kname

    def keys(self):
        return self.particles.keys()
    
    def optional_keys(self) -> List[str]:
        return list(set(self.keys()).intersection(self._optional_keys_in_particles))

    def prepare_input(self, survey: Survey) -> Tuple[str, pathlib.Path, Dict[str, Union[str,float,int]]]:
        """
            TODO
        """
        parfile, for_parfile = survey._write_parameter_file()
        temp_filename = self._write_ebf_files()
        return self.name_hash, parfile, for_parfile

    def _write_ebf_files(self):
        particlefile: pathlib.Path = self.pname
        kernelfile: pathlib.Path = self.kname
        inputhashfile: pathlib.Path = self._hashname
        inputhash = self._inputhash
        if ((inputhashfile.read_bytes() != inputhash # proceed if hashes don't match,
            if (particlefile.exists() and            # only if particlefile exists,
                kernelfile.exists() and              # and kernelfile exists,
                inputhashfile.exists())              # and inputhashfile exists,
            else True)                               # otherwise proceed if both don't exist
            if self.caching else True):              # -> proceed anyway if caching is False
            self.__write_particles(particlefile)
            self.__write_kernels(kernelfile)
            inputhashfile.write_bytes(inputhash)
        temp_filename = self.__prepare_nbody1(kernelfile, particlefile)
        return temp_filename

    @property
    def input_sorter(self) -> NDArray[np.int_]:
        return self._input_sorter
    
    @input_sorter.setter
    def input_sorter(self, value: Optional[NDArray[np.int_]]) -> None:
        if value is None:
            value: NDArray[np.int_] = self.__lex_partitionid_sorter
        else:
            pass # TODO check validity of input_sorter?
        self._input_sorter: NDArray[np.int_] = value
        if '_inputhash' in self.__dict__:
            del self._inputhash

    @cached_property
    def _inputhash(self) -> bytes:
        return hash_iterable(map(lambda array: array[self.input_sorter].copy(order='C'),
                                 itertools.chain(self.particles.values(),[self.kernels])))

    @property
    @mark_metadata_prop
    def hash(self) -> str:
        return self._inputhash.decode()

    @property
    def __lex_partitionid_sorter(self) -> NDArray[np.int_]:
        return np.lexsort((self.particles[self._partitionid],))

    @property
    def metadata(self) -> Dict[str, Any]:
        return self._metadata

    def __write_particles(self, pname: pathlib.Path):
        if not self.__input_files_exist:
            ebf.initialize(pname)
            for key in self._required_keys_in_particles:
                ebf.write(pname, f"/{key}", self.particles[key][self.input_sorter], 'a')
            for key in self._optional_keys_in_particles:
                ebf.write(pname, f"/{key}", self.particles[key][self.input_sorter] if key in self.keys() else np.zeros(self.length), 'a')
   
    def __write_kernels(self, kname: pathlib.Path):
        if not self.__input_files_exist:
            ebf.initialize(self.kname)
            ebf.write(kname, f"/{self._kernels}", self.k_factor*self.kernels[self.input_sorter], "a")
 
    def __prepare_nbody1(self, kname: pathlib.Path, pname: pathlib.Path):
        temp_dir = GALAXIA_NBODY1 / self.name_hash
        temp_dir.mkdir(parents=True, exist_ok=True)
        make_symlink(kname, temp_dir)
        make_symlink(pname, temp_dir)
        temp_filename = (GALAXIA_FILENAMES / self.name_hash).with_suffix('.txt')
        temp_filename.write_text(FILENAME_TEMPLATE.substitute(name=self.name_hash, pname=pname.name))
        return temp_filename
    
    @classmethod
    def __verify_particles(cls, particles: Dict[str, NDArray]):
        if cls._mass in particles and cls._massinit not in particles:
            from __metadata__ import __email__
            raise KeyError(f"BACKWARD INCOMPATIBILITY: please note of recent changes in the py-Galaxia-ananke implementation.\nFrom now on, the input must include the initial mass as well as the present-day mass of the parent particles. Also this quantity should be stored in the input dictionary under key '{cls._massinit}', while keeping the present-day mass under key '{cls._mass}'.\nPlease contact {__email__} if you have any question regarding that change.")
        compare_given_and_required(particles.keys(), cls._required_keys_in_particles, cls._optional_keys_in_particles,
                                   error_message="Given particle data covers wrong set of keys")
        confirm_equal_length_arrays_in_dict(particles, cls._massinit, error_message_dict_name='particles')
        # TODO check format, if dataframe-like
    
    @classmethod
    def __complete_particles(cls, particles: Dict[str, NDArray]):
        if cls._parentid not in particles:
            particles[cls._parentid] = np.arange(particles[cls._massinit].shape[0])
        if cls._partitionid not in particles:
            particles[cls._partitionid] = np.zeros(particles[cls._massinit].shape[0], dtype='int')
        # if cls._dform not in particles:
        #     particles[cls._dform] = 0*particles[cls._massinit]

    @classmethod
    def __conform_kernels(cls, kernels: NDArray) -> NDArray:
        if kernels.ndim == 2:
            if kernels.shape[1] == 1:
                kernels = kernels.flatten()
            elif kernels.shape[1] not in [2]:  # TODO in future, this could consider (Nx6) and (Nx21) cases for ellipsoidal kernels
                raise ValueError(f"Given 2D-array kernels is (Nx{kernels.shape[1]}) but should be either (Nx1) (position space) or (Nx2) (phase space)")
        if kernels.ndim not in [1,2]:
            raise ValueError(f"Given kernels ndim is {kernels.ndim} but should be either 1 (position space) or 2 (phase space)")
        return kernels

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
            
            Notes
            -----
            {particles_dictionary_description}
        """
        p = {}
        p[cls._pos] = 30*np.random.randn(n_parts, 3)
        p[cls._vel] = 50*np.random.randn(n_parts, 3)
        p[cls._massinit] = 5500 + 700*np.random.randn(n_parts)
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
    def make_dummy_kernels_input(cls, n_parts=10**5):
        """
            Generate an example dummy input kernels for Input
            made of randomly generated arrays.

            Parameters
            ----------
            n_parts : int
                Number of particles the example include. Default to 10**5.

            Returns
            -------
            kernels : array
                Dummy example input phase space kernels for Input.
        """
        kernels = np.exp([0.5,1] + 0.4*np.random.randn(n_parts,2))
        return kernels


Input.__init__.__doc__ = Input.__init__.__doc__.format(GALAXIA_TMP=GALAXIA_TMP,
                                                       DEFAULT_SIMNAME=DEFAULT_SIMNAME,
                                                       TTAGS_nres=FTTAGS.nres,
                                                       particles_dictionary_description=Input.particles_dictionary_description)

Input.make_dummy_particles_input.__func__.__doc__ = Input.make_dummy_particles_input.__doc__.format(
    particles_dictionary_description=Input.particles_dictionary_description)


if __name__ == '__main__':
    raise NotImplementedError()
