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
The isochrones of a photometric system, made from a folder of files that exists or from data.
"""
import re
import sys
from typing import Any, Dict, KeysView, Tuple, List, Callable, Sequence, Optional
from numpy.typing import NDArray
from functools import cached_property
import pathlib
import shutil

import numpy as np
from astropy.table import QTable

from .._constants import *
from ..utils import compare_given_and_required, tuple_itemgetter
from .IsochroneFile import IsochroneFile
from .Formatting import *
from .Formatting import Formatting

__all__ = ['Isochrone']


class Isochrone:
    """
    The isochrones of a photometric system: a folder of files, one for each metallicity, and the descriptor of those files
    that the backend reads.

    The folder is in a category, which is the folder that holds it. An isochrone is made from a folder that exists, or from
    its data, which is then written in the folder of the custom photometric systems of the installation. Ages need to be in
    log scale.
    """
    _zini = python_formatting._zini
    _age = python_formatting._age
    _mini = python_formatting._mini
    _mact = python_formatting._mact
    _lum = python_formatting._lum
    _teff = python_formatting._teff
    _grav = python_formatting._grav
    _required_keys = [_age, _mini, _mact, _lum, _teff, _grav]
    _required_metallicities = {0.0001, 0.0002, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.0012, 0.0014, 0.0016, 0.0018, 0.002, 0.0022, 0.0024, 0.0026, 0.003, 0.0034, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.024, 0.028, 0.03}
    _file_descriptor = "IsoFileDescriptor.txt"
    _magnitude_systems = ('Vega', 'AB', 'ST')
    _magnitude_system_line = re.compile(r'^#\s*Magnitude system\s*:\s*(.*?)\s*$', re.IGNORECASE)
    def __init__(self, *args: Any, overwrite: bool = False, magnitude_system: Optional[str] = None) -> None:
        """
        Make an isochrone from a folder that exists, or from its data.

        Parameters
        ----------
        *args
            The path of the folder of an isochrone that exists. Or the name of a new isochrone and its data, which is a
            dictionary of tables, one for each metallicity (the mass fraction Z): the folder is made in the custom
            photometric systems of the installation.
        overwrite : bool, optional
            Whether a new isochrone replaces one of the same name. Default to False.
        magnitude_system : str, optional
            The magnitude system of the magnitudes of a new isochrone, 'Vega', 'AB' or 'ST', whatever the case. That of an
            isochrone that exists is in its descriptor. Default to none, which says nothing of it in the descriptor.

        Raises
        ------
        TypeError
            If no argument is given, or more than two, or a magnitude system with the path of an isochrone that exists.
        FileExistsError
            If a new isochrone has the name of one that exists, and ``overwrite`` is not set.
        ValueError
            If the magnitude system is not one of the three, or if the data is refused (see ``_check_isochrone_data``):
            nothing is written then, and the isochrone that it was to replace is left as it was.
        """
        if not args:
            raise TypeError("Isochrone requires at least one argument")
        elif len(args) == 1:
            if magnitude_system is not None:
                raise TypeError("magnitude_system is for an isochrone that is made from its data: that of an isochrone that "
                                "exists is the line '# Magnitude system: ...' of its descriptor")
            self._path = pathlib.Path(args[0])
        elif len(args) == 2:
            magnitude_system = self._check_magnitude_system(magnitude_system)  # before anything is written or removed
            self._path = ISOCHRONES_PATH / CUSTOM_PHOTOCAT / args[0]
            if self.path.exists() and not overwrite:
                raise FileExistsError(f"Isochrone '{args[0]}' already exists at '{self.path}' (use `overwrite` kwarg to ignore)")
            iso_column_order = self._check_isochrone_data(args[1])  # before the isochrone that is to be replaced is removed
            if self.path.exists():
                shutil.rmtree(self.path)
            self._write_isochrone_files(args[1], iso_column_order, magnitude_system)
        else:
            raise TypeError(f"Too many arguments ({len(args)} given)")
    
    def __repr__(self) -> str:
        """
        Give the representation of the isochrone.

        Returns
        -------
        representation : str
            The class, with the category and the name of the isochrone.
        """
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['category', 'name']])
        return f'{cls}({description})'
    
    @classmethod
    def _check_magnitude_system(cls, value: Optional[str]) -> Optional[str]:
        """
        Give the magnitude system as the descriptors spell it, whatever the case it was written in.

        Parameters
        ----------
        value : str or None
            The magnitude system, which may have spaces around it.

        Returns
        -------
        magnitude_system : str or None
            'Vega', 'AB' or 'ST', or None if None was given.

        Raises
        ------
        ValueError
            If it is not one of those.
        """
        if value is None:
            return None
        known = {system.casefold(): system for system in cls._magnitude_systems}
        if not isinstance(value, str) or value.strip().casefold() not in known:
            allowed = ', '.join(map(repr, cls._magnitude_systems))
            raise ValueError(f"The magnitude system has to be one of {allowed}, not {value!r}")
        return known[value.strip().casefold()]

    def _write_isochrone_files(self, isochrone_data: dict, iso_column_order: list,
                               magnitude_system: Optional[str] = None) -> None:
        """
        Write the folder of the isochrone: the descriptor, and a file for each metallicity.

        If the data is refused or a file cannot be written, nothing is left of what was started.

        Parameters
        ----------
        isochrone_data : dict
            The tables, by metallicity.
        iso_column_order : list of str
            The columns to write, in order.
        magnitude_system : str, optional
            The magnitude system of the magnitudes, for the descriptor. Default to none.
        """
        self.path.mkdir(parents=True, exist_ok=True)
        try:
            self._write_file_descriptor(iso_column_order, magnitude_system)
            for feh, iso in isochrone_data.items():
                path = self._path / IsochroneFile._file_format.format(format(feh, '.6f'))
                _temp = IsochroneFile(path, iso[iso_column_order], isochrone=self)
        except BaseException:  # data that is refused, or files that cannot be written: nothing is left of what was started
            shutil.rmtree(self.path, ignore_errors=True)
            raise

    def _check_isochrone_data(self, isochrone_data: dict) -> list:
        """
        Check that the data of an isochrone is as it has to be, and give the order of its columns. Nothing is written.

        Parameters
        ----------
        isochrone_data : dict
            The tables, by metallicity.

        Returns
        -------
        iso_column_order : list of str
            The required columns, then the magnitudes in alphabetical order.

        Raises
        ------
        ValueError
            If the metallicities are not the required ones, the tables do not have the same columns, a required column is
            missing, there is no magnitude column, or the tables do not have the same ages as the backend reads them: it
            takes those of the first metallicity for all.
        """
        metallicities, headers = zip(*[(feh, list(iso.keys())) for feh, iso in isochrone_data.items()])
        compare_given_and_required(metallicities, self._required_metallicities, error_message="Given isochrone data covers wrong set of metallicities")
        check = []
        for header in headers:
            if header not in check: check.append(header)
        if len(check) > 1:
            raise ValueError("Given isochrone dataframes have unequal headers")
        header = check[0]
        remain = set(self._required_keys).difference(header)
        if remain:
            raise ValueError(f"Given isochrone dataframes have incomplete headers: missing {remain}")
        magnames = sorted(list(set(header).difference(self._required_keys)))
        if not magnames:
            raise ValueError(f"Given isochrone dataframes have no magnitude columns")
        # the backend takes the ages of the first metallicity for all of them, and only checks how many there are and the first and last ones
        ages = [self._ages_as_the_backend_reads_them(iso) for iso in isochrone_data.values()]
        for feh, others in zip(metallicities[1:], ages[1:]):
            if not np.array_equal(ages[0], others):
                raise ValueError(f"Given isochrone dataframes do not have the same ages: the metallicity {metallicities[0]} has {len(ages[0])} and "
                                 f"the metallicity {feh} has {len(others)}" + ("" if len(ages[0]) != len(others) else ", at other values"))
        return self._required_keys + magnames

    def _ages_as_the_backend_reads_them(self, table: Any) -> NDArray:
        """
        Give the distinct ages of a table as the backend counts them: as 32-bit floats, and leaving out the rows that have
        no current mass.

        Parameters
        ----------
        table : table-like
            The table of one metallicity.

        Returns
        -------
        ages : numpy.ndarray
            The ages, sorted.
        """
        with_mass = ~np.isnan(np.asarray(table[self._mact], dtype=float))
        return np.unique(np.asarray(table[self._age], dtype=np.float32)[with_mass])

    def _write_file_descriptor(self, iso_column_order: list, magnitude_system: Optional[str] = None) -> None:
        """
        Write the descriptor of the files of the isochrone.

        Its first line is the only one that the backend reads: the name, the number of columns, the number of the required
        ones (which is the index of the first magnitude) and the number of magnitudes, then the names of the magnitudes. A
        comment line gives the magnitude system, if there is one.

        Parameters
        ----------
        iso_column_order : list of str
            The columns of the files, in order.
        magnitude_system : str, optional
            The magnitude system, for the comment line. Default to none.
        """
        magnames = iso_column_order[len(self._required_keys):]
        lines = [f"Python_{self.name} {len(iso_column_order)} {len(self._required_keys)} {len(magnames)} {' '.join(magnames)}"]
        if magnitude_system is not None:
            lines.append(f"# Magnitude system: {magnitude_system}")  # a comment: the backend reads the first line only
        with self.file_descriptor_path.open('w') as f: f.write('\n'.join(lines) + "\n\n")

    @property
    def path(self) -> pathlib.Path:
        """The path of the folder of the isochrone."""
        return self._path
    
    @property
    def category(self) -> str:
        """The category of the isochrone: the name of the folder that holds its folder."""
        return self.path.parent.name
    
    @property
    def name(self) -> str:
        """The name of the isochrone: the name of its folder."""
        return self.path.name

    @property
    def file_descriptor_path(self) -> pathlib.Path:
        """The path of the descriptor of the files of the isochrone."""
        return self._path / self._file_descriptor

    @property
    def has_file_descriptor(self) -> bool:
        """Whether the folder of the isochrone has its descriptor."""
        return self.file_descriptor_path.exists()

    @cached_property
    def isochrone_files(self) -> List[IsochroneFile]:
        """The files of the isochrone, one for each metallicity, in the order of their names."""
        return list(map(lambda path: IsochroneFile(path, isochrone=self),
                        sorted(self._path.glob(IsochroneFile._file_format.format('*')))))
    
    @cached_property
    def file_descriptor_content(self) -> List[str]:
        """
        The words of the first line of the descriptor.

        Reading it raises NotImplementedError if the isochrone has no descriptor.
        """
        if self.has_file_descriptor:
            with open(self.file_descriptor_path,'r') as f: return f.readline().strip('\n').split()
        else:
            raise NotImplementedError("Photometric system doesn't have an IsoFileDescriptor.txt file")
    
    @cached_property
    def magnitude_system(self) -> Optional[str]:
        """
        The magnitude system of the magnitudes of the isochrone, 'Vega', 'AB' or 'ST', as the line '# Magnitude system: ...'
        that follows the first line of its descriptor says it, or None where the descriptor does not say it.
        """
        if not self.has_file_descriptor:
            return None
        for line in self.file_descriptor_path.read_text().splitlines()[1:]:  # after the line that the backend reads
            match = self._magnitude_system_line.match(line)
            if match:
                try:
                    return self._check_magnitude_system(match.group(1))
                except ValueError as error:
                    raise ValueError(f"{error} (in {self.file_descriptor_path})") from None
        return None

    # (nothing uses the two properties below: commented out, with the mark UNUSED, to be removed after v1)
    # UNUSED: @cached_property
    # UNUSED: def __fd_photosysname(self) -> str:
    # UNUSED:     return self.file_descriptor_content[0]
    
    # UNUSED: @cached_property
    # UNUSED: def __fd_fields(self) -> int:
    # UNUSED:     return int(self.file_descriptor_content[1])
    
    @cached_property
    def _fd_startid(self) -> int:
        """The index of the first magnitude among the columns of the files: the third word of the descriptor."""
        return int(self.file_descriptor_content[2])
    
    @cached_property
    def _fd_nmags(self) -> int:
        """The number of magnitudes: the fourth word of the descriptor."""
        return int(self.file_descriptor_content[3])
    
    @cached_property
    def magnitude_itemgetter(self) -> Callable[[Sequence], tuple]:
        """The getter that takes the magnitudes from the names of the columns of a file, as a tuple."""
        return tuple_itemgetter(*range(self._fd_startid, self._fd_startid+self._fd_nmags))

    @cached_property
    def mag_names(self) -> List[str]:
        """The names of the magnitudes: the words of the descriptor after the fourth."""
        return self.file_descriptor_content[4:]

    @cached_property
    def formatting(self) -> Formatting:
        """
        The layout of the files of the isochrone: where their columns are, as the version of CMD that wrote them put them.

        The systems of the padova category that CMD made have the layout of their version. What is made in Python, and a
        system that is not known, have the layout of Python.
        """
        if False: # sys.version_info >= (3,10):  # TODO uncomment below when upgrading officially to 3.10
            # match self.name:
            #     case ('APO__SDSS' | 'DCMC' | 'Generic__Johnson_UBVRIJHK' | 'Generic__Stroemgren') if self.category == 'padova':
            #         return cmd_2_1_formatting
            #     case 'GAIA__0' if self.category == 'padova':
            #         return cmd_3_0_formatting
            #     case ('CTIO__DECam' | 'LSST_DP0' | 'Roman' | 'Euclid' | 'JWST' | 'Pan_STARRS1' | 'LSST' | 'GAIA__DR2' | 'HST__WFC3' | '2MASS+Spitzer+WISE') if self.category == 'padova':
            #         return cmd_3_6_formatting
            #     case _:
            #         return python_formatting
            pass
        else:
            if self.name == 'GAIA__0' and self.category == 'padova':
                return cmd_3_0_formatting
            if self.name in ['APO__SDSS','DCMC','Generic__Johnson_UBVRIJHK','Generic__Stroemgren'] and self.category == 'padova':
                return cmd_2_1_formatting
            if self.name in ['CTIO__DECam','LSST_DP0','Roman','Euclid','JWST__NIRCam','JWST__MIRI','Pan_STARRS1','LSST','GAIA__DR2','HST__WFC3','2MASS+Spitzer+WISE'] and self.category == 'padova':
                return cmd_3_6_formatting
            else:
                return python_formatting  # what is made in Python, and a system that is not known

    @cached_property
    def qtables_dictionary(self) -> Dict[float, Dict[float, QTable]]:  # by log10(Z/Zsun) (Zsun of the layout), then by age
        """
        The tables of the isochrone, with the names and the units of the package, by metallicity and then by age.

        The metallicity is the decimal logarithm of its ratio to the solar one, which is the one of the layout.
        """
        return {self.formatting.metallicity_converter(iso_file.metallicity): self.formatting.qtable_per_age_from_isochronefile(iso_file)
                for iso_file in self.isochrone_files}
    
    @cached_property
    def qtables_unique_ages_dictionary(self) -> Dict[float, KeysView[float]]:
        """The ages of the tables of each metallicity, by metallicity."""
        return {metallicity: qtables_per_age_dict.keys() for metallicity, qtables_per_age_dict in self.qtables_dictionary.items()}

    @cached_property
    def consolidated_points_and_values(self) -> Tuple[NDArray,NDArray]:
        """
        The points and the values that the magnitudes are interpolated from: all the tables of the isochrone, completed.

        The points are the metallicity, the age and the initial mass of each row, and the values are its magnitudes. The
        tables only hold the stars that are alive at their age, so the space of mass and age is completed in a simplistic
        way, which is experimental: a star that is a white dwarf (a luminosity of -9.999) is given at the oldest age too,
        and for an age that does not end with a white dwarf the stars that are more massive than the ones of the next age
        are added at that age, with the magnitudes set to 999, before they are given at the oldest.

        It is a tuple of the array of the points, of shape (n, 3), and the array of the values, of shape (n, number of
        magnitudes). Reading it raises ValueError if the oldest age of a metallicity does not end with a white dwarf.
        """
        list_of_points = []
        list_of_values = []
        for metallicity, unique_ages_view in self.qtables_unique_ages_dictionary.items():
            unique_ages = list(unique_ages_view)
            qtables = self.qtables_dictionary[metallicity]
            if qtables[unique_ages[-1]][-1][self._lum].value != -9.999:  # the same marker as below for the white dwarfs
                raise ValueError(f"The oldest age of the isochrone (log age {unique_ages[-1]}) at metallicity {metallicity} does not end with a white dwarf "
                                 f"(luminosity -9.999): completing the (mass, age) space needs an oldest age old enough for its most massive star to have "
                                 f"died as a white dwarf, please add older ages to the isochrone")
            for i_age, age in enumerate(unique_ages):
                qtable = qtables[age]
                list_of_points.append(qtable[[self._zini, self._age, self._mini]].to_pandas().to_numpy())
                list_of_values.append(qtable[self.mag_names].to_pandas().to_numpy())
                if qtable[-1][self._lum].value == -9.999:  # bit of a stretch for white dwarf case, but that will do for now for what I need
                    if age != unique_ages[-1]:
                        list_of_points.append(list_of_points[-1][[-1]].copy())
                        list_of_points[-1][:,1] = unique_ages[-1]
                        list_of_values.append(list_of_values[-1][[-1]].copy())
                else:  # also a stretch, no neutron star, just black holes
                    next_age = unique_ages[i_age+1]
                    max_mini_next_age = qtables[next_age][self._mini].max()
                    dead_next_age_qtable = qtable[qtable[self._mini] > max_mini_next_age]
                    list_of_points.append(dead_next_age_qtable[[self._zini, self._age, self._mini]].to_pandas().to_numpy())
                    list_of_points[-1][:,1] = next_age
                    list_of_values.append(dead_next_age_qtable[self.mag_names].to_pandas().to_numpy())
                    list_of_values[-1][:,:] = 999
                    list_of_points.append(list_of_points[-1][[-1]].copy())
                    list_of_points[-1][:,1] = unique_ages[-1]
                    list_of_values.append(list_of_values[-1][[-1]].copy())
        list_of_lengths = np.array(list(map(len, list_of_points)))
        successive_indices = np.hstack([0,np.cumsum(list_of_lengths)])
        points = np.empty((successive_indices[-1], 3))
        values = np.empty((successive_indices[-1], len(self.mag_names)))
        for first, last, sub_points, sub_values in zip(successive_indices[:-1], successive_indices[1:], list_of_points, list_of_values):
            points[first:last] = sub_points
            values[first:last] = sub_values
        return points, values


if __name__ == '__main__':
    raise NotImplementedError()
