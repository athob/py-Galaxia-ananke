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
Contains the Input class definition

Please note that this module is private. The Input class is
available in the main ``galaxia_ananke`` namespace - use that instead.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Any, Optional, Union, Tuple, List, Set, Dict, OrderedDict, Iterator, KeysView
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
from .utils import classproperty, format_docstring, in_description_order, make_symlink, compare_given_and_required, confirm_equal_length_arrays_in_dict, lexicalorder_dict, mark_metadata_prop, collect_metadata_marked_properties, hash_iterable
from .photometry.PhotoSystem import PhotoSystem

if TYPE_CHECKING:
    from .Survey import Survey

__all__ = ['Input']


@collect_metadata_marked_properties
class Input:
    """
    Driver that stores the particles and their kernels, and prepares them as the input of the backend.

    The input is made from the particles and the kernels, or from the files of an input that was prepared before. It writes
    its files in a directory, with a hash of its data in their names, and it hands the backend what it needs to run on them.
    """
    _metadata: Dict[str, Any]  # a property that the decorator above adds
    _input_sorter: NDArray[np.int_]  # set by the setter of input_sorter, which a Survey calls: a new Input has none yet (M14)
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
    def __init__(self, *args: Any, **kwargs: Any) -> None:
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
            *args
                The particles and the kernels, in this order, which are
                described below. They are not given when ``pname`` and
                ``kname`` are.

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
                Whether the work that is already done for the same input, in the
                same directory, is not done again: the files of the input are
                not written again if their hash is the one of the data, and the
                steps of the survey that a state file says are done are skipped.
                This is experimental, and a warning says so. Default to False.

            append_hash : bool
                Whether the names of the files of the input end with the first
                7 characters of the hash of the data, so that data that differ
                do not share their files. Default to True.

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
        self.__append_hash: bool = kwargs.get('append_hash', True)
        # check which args/kwargs signature case and populate kwargs accordingly
        if args:
            if len(args) not in [2]:  # TODO mix & match args & kwargs for particles and kernels
                raise TypeError(f"Input takes the particles and the kernels as its two arguments, and not {len(args)}: please consult help of the Input constructor")
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
            self.__input_files_exist = True
        else:
            raise ValueError("Wrong signature: please consult help of the Input constructor")
        # verify and assign attributes based on kwargs
        self.__particles: OrderedDict[str, NDArray] = lexicalorder_dict({key: np.asarray(value) for key, value in kwargs['particles'].items()})
        self.__verify_particles(self.particles)
        self.__complete_particles(self.particles)
        self.__kernels: NDArray = self.__conform_kernels(np.asarray(kwargs['kernels']), self.particles[self._massinit].shape[0])
        self.__input_dir: pathlib.Path = pathlib.Path(kwargs.get('input_dir', GALAXIA_TMP))
        self.__name: str = kwargs.get('name', DEFAULT_SIMNAME)
        self.__pname: Optional[pathlib.Path] = kwargs.get('pname', None)
        self.__kname: Optional[pathlib.Path] = kwargs.get('kname', None)
        self.__ngb: int = kwargs.get('ngb', DEFAULTS_FOR_PARFILE[FTTAGS.nres])  # TODO Galaxia makes no use of it, need to remove in the future
        self.__k_factor: Union[float, NDArray] = kwargs.get('k_factor', 1.)

    @classproperty
    def particles_dictionary_description(cls) -> str:
        """
        The description of the dictionary of the particles that the documentation of the input includes.

        It lists the properties that the particles must have, and those that they may have, each with its key.
        """
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
                        for key, desc in in_description_order(Input._required_properties)]),
                   _optional_properties=''.join(
                       [f"\n            * {desc} via key ``{str(key)}``"
                        for key, desc in in_description_order(Input._optional_properties)]))
        return description

    @classproperty
    def _required_properties(cls) -> Set[Tuple[str, str]]:
        """The properties that the particles must have, as pairs of their key and their description."""
        return {
            cls._position_prop,
            cls._velocity_prop,
            cls._massinitial_prop,
            cls._age_prop,
            cls._metallicity_prop
            }

    @classproperty
    def _optional_properties(cls) -> Set[Tuple[str, str]]:
        """The properties that the particles may have, as pairs of their key and their description."""
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
    def _required_keys_in_particles(cls) -> Set[str]:
        """The keys of the properties that the particles must have."""
        return {_property[0] for _property in cls._required_properties}

    @classproperty
    def _optional_keys_in_particles(cls) -> Set[str]:
        """The keys of the properties that the particles may have."""
        return {_property[0] for _property in cls._optional_properties}

    @classproperty
    def all_possible_keys_in_particles(cls) -> Set[str]:
        """The keys that the particles may have: the required ones and the optional ones."""
        return cls._required_keys_in_particles.union(cls._optional_keys_in_particles)
    
    @classproperty
    def _pos(cls) -> str:
        """The key of the positions in the dictionary of the particles."""
        return cls._position_prop[0]

    @classproperty
    def _vel(cls) -> str:
        """The key of the velocities in the dictionary of the particles."""
        return cls._velocity_prop[0]

    @classproperty
    def _mass(cls) -> str:
        """The key of the current masses in the dictionary of the particles."""
        return cls._masscurrent_prop[0]
    
    @classproperty
    def _massinit(cls) -> str:
        """The key of the initial masses in the dictionary of the particles."""
        return cls._massinitial_prop[0]

    @classproperty
    def _age(cls) -> str:
        """The key of the ages in the dictionary of the particles."""
        return cls._age_prop[0]

    @classproperty
    def _feh(cls) -> str:
        """The key of the metallicities in the dictionary of the particles."""
        return cls._metallicity_prop[0]

    @classproperty
    def _He(cls) -> str:
        """The key of the helium abundances in the dictionary of the particles."""
        return cls._heliumabundance_prop[0]

    @classproperty
    def _C(cls) -> str:
        """The key of the carbon abundances in the dictionary of the particles."""
        return cls._carbonabundance_prop[0]

    @classproperty
    def _N(cls) -> str:
        """The key of the nitrogen abundances in the dictionary of the particles."""
        return cls._nitrogenabundance_prop[0]

    @classproperty
    def _O(cls) -> str:
        """The key of the oxygen abundances in the dictionary of the particles."""
        return cls._oxygenabundance_prop[0]

    @classproperty
    def _Ne(cls) -> str:
        """The key of the neon abundances in the dictionary of the particles."""
        return cls._neonabundance_prop[0]

    @classproperty
    def _Mg(cls) -> str:
        """The key of the magnesium abundances in the dictionary of the particles."""
        return cls._magnesiumabundance_prop[0]

    @classproperty
    def _Si(cls) -> str:
        """The key of the silicon abundances in the dictionary of the particles."""
        return cls._siliconabundance_prop[0]

    @classproperty
    def _S(cls) -> str:
        """The key of the sulphur abundances in the dictionary of the particles."""
        return cls._sulphurabundance_prop[0]

    @classproperty
    def _Ca(cls) -> str:
        """The key of the calcium abundances in the dictionary of the particles."""
        return cls._calciumabundance_prop[0]

    @classproperty
    def _elem_list(cls) -> List[str]:
        """The keys of the abundances of the elements, other than the metallicity and the alpha abundance."""
        return [cls._He, cls._C, cls._N, cls._O, cls._Ne, cls._Mg, cls._Si, cls._S, cls._Ca]

    @classproperty
    def _alph(cls) -> str:
        """The key of the alpha abundances in the dictionary of the particles."""
        return cls._alphaabundance_prop[0]

    @classproperty
    def _parentid(cls) -> str:
        """The key of the indices of the parent particles in the dictionary of the particles."""
        return cls._parentindex_prop[0]

    @classproperty
    def _partitionid(cls) -> str:
        """The key of the indices of the partitions that contain the particles in the dictionary of the particles."""
        return cls._partitionindex_prop[0]

    @classproperty
    def _dform(cls) -> str:
        """The key of the formation distances of the parent particles in the dictionary of the particles."""
        return cls._formationdistance_prop[0]

    @classproperty
    def _pop_id(cls) -> str:
        """The key of the indices of the populations of the parent particles in the dictionary of the particles."""
        return cls._populationindex_prop[0]

    @classproperty
    def _kernels(cls) -> str:
        """The key of the kernels in the dictionary of the particles."""
        return cls._kernels_prop[0]

    @property
    def caching(self) -> bool:
        """
        Whether the work that is already done for the same input is not done again.

        It can be set: a true value warns that this is experimental.
        """
        return self.__caching

    @caching.setter
    def caching(self, value: bool) -> None:
        if value:
            warn(f"You have requested caching mode, be aware this feature is currently experimental and may result in unintended behaviour.", DeprecationWarning, stacklevel=2)
        self.__caching: bool = value
    
    @property
    def particles(self) -> Dict[str, NDArray]:
        """
        The particles: the arrays of their properties, by key, in the alphabetical order of the keys.

        The keys that the particles had to have are there, and so are the indices of the parents and of the partitions,
        which are made if they were not given.
        """
        return self.__particles
    
    @property
    def particle_parentids(self) -> NDArray:
        """The index of the parent particle of each particle."""
        return self.particles[self._parentid]
    
    @property
    def parent_particles_data_generator(self) -> Iterator[Tuple[str, NDArray]]:
        """
        The arrays of the particles that the catalogue keeps for the parents, as pairs of their key and their array.

        They are the required properties and then the kernels, which are multiplied by the scaling factor.
        """
        return itertools.chain(
            ((key, self.particles[key]) for key in self._required_keys_in_particles),
            ((self._kernels, self.__actual_kernels) for _ in range(1))
            )
    
    @property
    def length(self) -> int:
        """The number of particles."""
        return len(self.particles[self._massinit])
    
    @property
    def hdim(self) -> int:  # TODO how to adapt for ellipsoidal kernels?
        """The dimension of the space of the kernels: 3 for kernels in position, and 6 for kernels in phase space."""
        return 3 if self.kernels.ndim == 1 else 6

    @property
    @mark_metadata_prop
    def name(self) -> str:
        """The name that the files of the input are given."""
        return self.__name
    
    @property
    def append_hash(self) -> bool:
        """Whether the names of the files of the input end with the beginning of the hash of the data."""
        return self.__append_hash

    @property
    def name_hash(self) -> str:
        """
        The name of the files of the input: the name, followed by an underscore and
        the first 7 characters of the hash if it is appended.
        """
        return self.name + (f"_{self.hash[:7]}" if self.append_hash else "")

    @property
    def ngb(self) -> int:
        """The number of neighbouring particles that the backend is asked to consider."""
        return self.__ngb
    
    @property
    def k_factor(self) -> Union[float, NDArray]:
        """The scaling factor of the sizes of the kernels."""
        return self.__k_factor

    @property
    def _input_dir(self) -> pathlib.Path:
        """The directory of the files of the input."""
        return self.__input_dir

    @property
    def kernels(self) -> NDArray:
        """The kernels, as they were given: one length for each particle, or two in phase space."""
        return self.__kernels
    
    @cached_property
    def __actual_kernels(self) -> NDArray:
        """The kernels that the backend is given: the given ones, times the scaling factor."""
        return self.k_factor * self.kernels
    
    @property
    def _base_inputfile(self) -> pathlib.Path:  # TODO what if pname and kname already exist (case where input args are pname and kname)?
        """The path of the files of the input without their suffix: in the directory, with the name and its hash."""
        return self._input_dir / self.name_hash

    @property
    def _hashname(self) -> pathlib.Path:
        """The path of the file that holds the hash of the input, next to the others."""
        return self._base_inputfile.with_name(f"{self._base_inputfile.name}.{HASH_EXT}")      # added to the name, which can hold dots

    @property
    def pname(self) -> pathlib.Path:
        """The path of the file of the particles: the one that was given, or one in the directory, named after the input."""
        if self.__pname is None:
            self.__pname = self._base_inputfile.with_name(f"{self._base_inputfile.name}.ebf")
        return self.__pname

    @property
    def kname(self) -> pathlib.Path:
        """
        The path of the file of the kernels: the one that was given, or one named after the file of the particles,
        with the dimension of the kernels and the number of neighbours.
        """
        if self.__kname is None:
            self.__kname = self.pname.with_name(f"{self.pname.stem}_d{self.hdim}n{self.ngb}_den.ebf")
        return self.__kname

    def keys(self) -> KeysView[str]:
        """
        Give the keys of the particles.

        Returns
        -------
        keys : collections.abc.KeysView
            The keys.
        """
        return self.particles.keys()
    
    def optional_keys(self) -> List[str]:
        """
        Give the optional keys that the particles have.

        Returns
        -------
        keys : list of str
            The keys that the particles may have, and have.
        """
        return list(set(self.keys()).intersection(self._optional_keys_in_particles))

    def prepare_input(self, survey: Survey) -> Tuple[str, pathlib.Path, Dict[str, Union[str,float,int]]]:
        """
        Write what the backend needs to run on the input: the parameter file of the survey, and the files of the input.

        Parameters
        ----------
        survey : Survey
            The survey that the input is prepared for, which writes its parameter file.

        Returns
        -------
        name_hash : str
            The name of the files of the input, which the backend is given.
        parfile : pathlib.Path
            The parameter file.
        parameters : dict
            The parameters of the survey that are in the parameter file.
        """
        parfile, for_parfile = survey._write_parameter_file()
        temp_filename = self._write_ebf_files()
        return self.name_hash, parfile, for_parfile

    def _write_ebf_files(self) -> pathlib.Path:
        """
        Write the files of the particles and of the kernels, and what makes the backend find them.

        When the input caches its results, the files are not written again if they exist and their hash is the one of the
        data.

        Returns
        -------
        filename : pathlib.Path
            The file that names the input for the backend.
        """
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
        """
        The indices that put the particles, and their kernels, in the order in which they are written and hashed.

        It can be set: ``None`` gives the order of the partition ids of the particles, and the hash is computed again.
        """
        return self._input_sorter
    
    @input_sorter.setter
    def input_sorter(self, value: Optional[NDArray[np.int_]]) -> None:
        if value is None:
            value = self.__lex_partitionid_sorter
        else:
            pass # TODO check validity of input_sorter?
        self._input_sorter = value
        if '_inputhash' in self.__dict__:
            del self._inputhash

    @cached_property
    def _inputhash(self) -> bytes:
        """The hash of the data of the input, with the particles and the kernels in the order of the sorter."""
        return hash_iterable(map(lambda array: array[self.input_sorter].copy(order='C'),
                                 itertools.chain(self.particles.values(),[self.__actual_kernels])))

    @property
    @mark_metadata_prop
    def hash(self) -> str:
        """The hash of the data of the input, as a string."""
        return self._inputhash.decode()

    @property
    def __lex_partitionid_sorter(self) -> NDArray[np.int_]:
        """The indices that put the particles in the order of their partition ids."""
        return np.lexsort((self.particles[self._partitionid],))

    @property
    def metadata(self) -> Dict[str, Any]:
        """The metadata of the input: the values of its properties that are marked for it."""
        return self._metadata

    def __write_particles(self, pname: pathlib.Path) -> None:
        """
        Write the particles in a file, in the order of the sorter.

        The required properties are written, and the optional ones too, as zeros for those that the particles do not have.
        Nothing is written when the files were given and exist.

        Parameters
        ----------
        pname : pathlib.Path
            The file of the particles.
        """
        if not self.__input_files_exist:
            ebf.initialize(pname)
            for key in self._required_keys_in_particles:
                ebf.write(pname, f"/{key}", self.particles[key][self.input_sorter], 'a')
            for key in self._optional_keys_in_particles:
                ebf.write(pname, f"/{key}", self.particles[key][self.input_sorter] if key in self.keys() else np.zeros(self.length), 'a')
   
    def __write_kernels(self, kname: pathlib.Path) -> None:
        """
        Write the kernels, times the scaling factor, in the order of the sorter.

        Nothing is written when the files were given and exist.

        Parameters
        ----------
        kname : pathlib.Path
            The file of the kernels.
        """
        if not self.__input_files_exist:
            ebf.initialize(self.kname)
            ebf.write(kname, f"/{self._kernels}", self.__actual_kernels[self.input_sorter], "a")
 
    def __prepare_nbody1(self, kname: pathlib.Path, pname: pathlib.Path) -> pathlib.Path:
        """
        Make the backend find the files: link them in a folder of its own input, and write the file that names the input.

        Parameters
        ----------
        kname : pathlib.Path
            The file of the kernels.
        pname : pathlib.Path
            The file of the particles.

        Returns
        -------
        filename : pathlib.Path
            The file that names the input for the backend.
        """
        temp_dir = GALAXIA_NBODY1 / self.name_hash
        temp_dir.mkdir(parents=True, exist_ok=True)
        make_symlink(kname, temp_dir)
        make_symlink(pname, temp_dir)
        temp_filename = GALAXIA_FILENAMES / f"{self.name_hash}.txt"
        temp_filename.write_text(FILENAME_TEMPLATE.substitute(name=self.name_hash, pname=pname.name))
        return temp_filename
    
    @classmethod
    def __verify_particles(cls, particles: Dict[str, NDArray]) -> None:
        """
        Check the particles.

        Parameters
        ----------
        particles : dict
            The arrays of the particles, by key.

        Raises
        ------
        KeyError
            If the current mass is given without the initial mass, which an earlier version took instead.
        ValueError
            If the keys are not the required ones with optional ones only besides, or the arrays do not have the same
            length.
        """
        if cls._mass in particles and cls._massinit not in particles:
            from .__metadata__ import __email__
            raise KeyError(f"BACKWARD INCOMPATIBILITY: please note of recent changes in the py-Galaxia-ananke implementation.\nFrom now on, the input must include the initial mass as well as the present-day mass of the parent particles. Also this quantity should be stored in the input dictionary under key '{cls._massinit}', while keeping the present-day mass under key '{cls._mass}'.\nPlease contact {__email__} if you have any question regarding that change.")
        compare_given_and_required(particles.keys(), cls._required_keys_in_particles, cls._optional_keys_in_particles,
                                   error_message="Given particle data covers wrong set of keys")
        confirm_equal_length_arrays_in_dict(particles, cls._massinit, error_message_dict_name='particles')
        # TODO check format, if dataframe-like
    
    @classmethod
    def __complete_particles(cls, particles: Dict[str, NDArray]) -> None:
        """
        Give the particles what they may not have: the indices of the parents, which are their own, and the partitions,
        which are all the same one.

        Parameters
        ----------
        particles : dict
            The arrays of the particles, by key, which is changed in place.
        """
        if cls._parentid not in particles:
            particles[cls._parentid] = np.arange(particles[cls._massinit].shape[0])
        if cls._partitionid not in particles:
            particles[cls._partitionid] = np.zeros(particles[cls._massinit].shape[0], dtype='int')
        # if cls._dform not in particles:
        #     particles[cls._dform] = 0*particles[cls._massinit]

    @classmethod
    def __conform_kernels(cls, kernels: NDArray, n_particles: int) -> NDArray:
        """
        Check the kernels, and give them as one array of lengths for the positions or two for the phase space.

        Parameters
        ----------
        kernels : numpy.ndarray
            The kernels.
        n_particles : int
            The number of particles, which the kernels have to match.

        Returns
        -------
        kernels : numpy.ndarray
            The kernels, flattened if they had a single column.

        Raises
        ------
        ValueError
            If they have another shape than N, Nx1 or Nx2, or are not as many as the particles.
        """
        if kernels.ndim == 2:
            if kernels.shape[1] == 1:
                kernels = kernels.flatten()
            elif kernels.shape[1] not in [2]:  # TODO in future, this could consider (Nx6) and (Nx21) cases for ellipsoidal kernels
                raise ValueError(f"Given 2D-array kernels is (Nx{kernels.shape[1]}) but should be either (Nx1) (position space) or (Nx2) (phase space)")
        if kernels.ndim not in [1,2]:
            raise ValueError(f"Given kernels ndim is {kernels.ndim} but should be either 1 (position space) or 2 (phase space)")
        if kernels.shape[0] != n_particles:
            raise ValueError(f"Given kernels are {kernels.shape[0]} but should be as many as the particles, {n_particles}")
        return kernels

    @classmethod
    def make_dummy_particles_input(cls, n_parts: int = 10**5) -> Dict[str, NDArray]:
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
        p: Dict[str, NDArray] = {}
        p[cls._pos] = 30*np.random.randn(n_parts, 3)
        p[cls._vel] = 50*np.random.randn(n_parts, 3)
        p[cls._massinit] = 5500 + 700*np.random.randn(n_parts)
        p[cls._age] = 9.7 + 0.4*np.random.randn(n_parts)
        p[cls._feh] = -0.7 + 0.4*np.random.randn(n_parts)
        # p[cls._mass] = np.clip(0.7 + 0.05*np.random.randn(n_parts), 0, 1)*p[cls._massinit]
        for el in cls._elem_list:
            p[el] = -0.6 + 0.4*np.random.randn(n_parts)
        p[cls._alph] = p[cls._Mg] - p[cls._feh]
        p[cls._parentid] = np.arange(n_parts).astype('int')
        p[cls._partitionid] = np.zeros(n_parts, dtype='int')
        p[cls._dform] = np.zeros(n_parts, dtype='float32')
        p[cls._pop_id] = np.zeros(n_parts, dtype='int')
        return p

    @classmethod
    def make_dummy_kernels_input(cls, n_parts: int = 10**5) -> NDArray:
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


format_docstring(Input.__init__,
                 GALAXIA_TMP=GALAXIA_TMP,
                 DEFAULT_SIMNAME=DEFAULT_SIMNAME,
                 TTAGS_nres=DEFAULTS_FOR_PARFILE[FTTAGS.nres],
                 particles_dictionary_description=Input.particles_dictionary_description)

format_docstring(Input.make_dummy_particles_input,
                 particles_dictionary_description=Input.particles_dictionary_description)


if __name__ == '__main__':
    raise NotImplementedError()
