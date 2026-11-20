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
Docstring
"""
import re
import sys
from typing import Tuple, List, Callable, Sequence, Optional
from numpy.typing import NDArray
from functools import cached_property
import pathlib
import shutil

import numpy as np

from .._constants import *
from ..utils import compare_given_and_required, tuple_itemgetter
from .IsochroneFile import IsochroneFile
from .Formatting import *

__all__ = ['Isochrone']


class Isochrone:
    """
    Ages need to be in log scale
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
    def __init__(self, *args, overwrite=False, magnitude_system: Optional[str] = None) -> None:
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
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['category', 'name']])
        return f'{cls}({description})'
    
    @classmethod
    def _check_magnitude_system(cls, value: Optional[str]) -> Optional[str]:
        """The magnitude system as it is spelled in the descriptors ('Vega', 'AB' or 'ST'), whatever the case, or None."""
        if value is None:
            return None
        known = {system.casefold(): system for system in cls._magnitude_systems}
        if not isinstance(value, str) or value.strip().casefold() not in known:
            allowed = ', '.join(map(repr, cls._magnitude_systems))
            raise ValueError(f"The magnitude system has to be one of {allowed}, not {value!r}")
        return known[value.strip().casefold()]

    def _write_isochrone_files(self, isochrone_data: dict, iso_column_order: list, magnitude_system: Optional[str] = None):
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
        """The order of the columns of the files of the isochrone, if its data is as it must be (it is refused if not): nothing is written."""
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

    def _ages_as_the_backend_reads_them(self, table) -> NDArray:
        """The distinct ages of a table as the backend counts them: as floats, and leaving out the rows that have no current mass."""
        with_mass = ~np.isnan(np.asarray(table[self._mact], dtype=float))
        return np.unique(np.asarray(table[self._age], dtype=np.float32)[with_mass])

    def _write_file_descriptor(self, iso_column_order: list, magnitude_system: Optional[str] = None):
        magnames = iso_column_order[len(self._required_keys):]
        lines = [f"Python_{self.name} {len(iso_column_order)} {len(self._required_keys)} {len(magnames)} {' '.join(magnames)}"]
        if magnitude_system is not None:
            lines.append(f"# Magnitude system: {magnitude_system}")  # a comment: the backend reads the first line only
        with self.file_descriptor_path.open('w') as f: f.write('\n'.join(lines) + "\n\n")

    @property
    def path(self):
        return self._path
    
    @property
    def category(self):
        return self.path.parent.name
    
    @property
    def name(self):
        return self.path.name

    @property
    def file_descriptor_path(self):
        return self._path / self._file_descriptor

    @property
    def has_file_descriptor(self):
        return self.file_descriptor_path.exists()

    @cached_property
    def isochrone_files(self):
        return list(map(lambda path: IsochroneFile(path, isochrone=self),
                        sorted(self._path.glob(IsochroneFile._file_format.format('*')))))
    
    @cached_property
    def file_descriptor_content(self) -> List[str]:
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
        return int(self.file_descriptor_content[2])
    
    @cached_property
    def _fd_nmags(self) -> int:
        return int(self.file_descriptor_content[3])
    
    @cached_property
    def magnitude_itemgetter(self) -> Callable[[Sequence], tuple]:
        return tuple_itemgetter(*range(self._fd_startid, self._fd_startid+self._fd_nmags))

    @cached_property
    def mag_names(self) -> List[str]:
        return self.file_descriptor_content[4:]

    @cached_property
    def formatting(self):
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
    def qtables_dictionary(self):  # by log10(Z/Zsun), with the Zsun of the layout
        return {self.formatting.metallicity_converter(iso_file.metallicity): self.formatting.qtable_per_age_from_isochronefile(iso_file)
                for iso_file in self.isochrone_files}
    
    @cached_property
    def qtables_unique_ages_dictionary(self):
        return {metallicity: qtables_per_age_dict.keys() for metallicity, qtables_per_age_dict in self.qtables_dictionary.items()}

    @cached_property
    def consolidated_points_and_values(self) -> Tuple[NDArray,NDArray]:
        list_of_points = []
        list_of_values = []
        for metallicity, unique_ages in self.qtables_unique_ages_dictionary.items():
            unique_ages = list(unique_ages)
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
                    # copies, as these two are changed below and pandas 3 gives read-only views of the data of a frame
                    dead_points = dead_next_age_qtable[[self._zini, self._age, self._mini]]
                    list_of_points.append(dead_points.to_pandas().to_numpy(copy=True))
                    list_of_points[-1][:,1] = next_age
                    list_of_values.append(dead_next_age_qtable[self.mag_names].to_pandas().to_numpy(copy=True))
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
