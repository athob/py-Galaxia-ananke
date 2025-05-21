#!/usr/bin/env python
"""
Docstring
"""
import sys
from typing import Tuple, List
from numpy.typing import NDArray
from functools import cached_property
from operator import itemgetter
import pathlib
import shutil

import numpy as np

from .._constants import *
from ..utils import compare_given_and_required
from .IsochroneFile import IsochroneFile
from .Formatting import *

__all__ = ['Isochrone']


class Isochrone:
    """
    Ages need to be in log scale
    """
    _zini = default_formatting._zini
    _age = default_formatting._age
    _mini = default_formatting._mini
    _mact = default_formatting._mact
    _lum = default_formatting._lum
    _teff = default_formatting._teff
    _grav = default_formatting._grav
    _required_keys = [_age, _mini, _mact, _lum, _teff, _grav]
    _required_metallicities = {0.0001, 0.0002, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.0012, 0.0014, 0.0016, 0.0018, 0.002, 0.0022, 0.0024, 0.0026, 0.003, 0.0034, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.024, 0.028, 0.03}
    _file_descriptor = "IsoFileDescriptor.txt"
    def __init__(self, *args, overwrite=False) -> None:
        if not args:
            raise TypeError("Isochrone requires at least one argument")
        elif len(args) == 1:
            self._path = pathlib.Path(args[0])
        elif len(args) == 2:
            self._path = ISOCHRONES_PATH / CUSTOM_PHOTOCAT / args[0]
            if self.path.exists():
                if overwrite:
                    shutil.rmtree(self.path)
                else:
                    raise FileExistsError(f"Isochrone '{args[0]}' already exists at '{self.path}' (use `overwrite` kwarg to ignore)")
            self._write_isochrone_files(args[1])
        else:
            raise TypeError(f"Too many arguments ({len(args)} given)")
    
    def __repr__(self) -> str:
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['category', 'name']])
        return f'{cls}({description})'
    
    def _write_isochrone_files(self, isochrone_data: dict):
        self.path.mkdir(parents=True, exist_ok=True)
        iso_column_order = self._write_file_descriptor(isochrone_data)
        for feh, iso in isochrone_data.items():
            path = self._path / IsochroneFile._file_format.format(format(feh, '.6f'))
            _temp = IsochroneFile(path, iso[iso_column_order], isochrone=self)

    def _write_file_descriptor(self, isochrone_data):
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
        iso_column_order = self._required_keys + magnames
        with self.file_descriptor_path.open('w') as f: f.write(
            f"Python_{self.name} {len(iso_column_order)} {len(self._required_keys)} {len(magnames)} {' '.join(magnames)}\n\n")
        return iso_column_order

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
    def __fd_photosysname(self) -> str:
        return self.file_descriptor_content[0]
    
    @cached_property
    def __fd_fields(self) -> int:
        return int(self.file_descriptor_content[1])
    
    @cached_property
    def _fd_startid(self) -> int:
        return int(self.file_descriptor_content[2])
    
    @cached_property
    def _fd_nmags(self) -> int:
        return int(self.file_descriptor_content[3])
    
    @cached_property
    def magnitude_itemgetter(self) -> itemgetter:
        return itemgetter(*range(self._fd_startid, self._fd_startid+self._fd_nmags))

    @cached_property
    def mag_names(self) -> List[str]:
        return self.file_descriptor_content[4:]

    @cached_property
    def formatting(self):
        if False: # sys.version_info >= (3,10):  # TODO uncomment below when upgrading officially to 3.10
            # match self.name:
            #     case 'WFIRST+HST__WFC3' if self.category == 'padova':
            #         return oldpadova_fomatting_withlogage
            #     case ('WFIRST' | 'LSST' | 'GAIA__0' | 'GAIA__DR2') if self.category == 'padova':
            #         return oldpadova_fomatting
            #     case ('CTIO__DECam' | 'LSST_DP0' | 'Roman' | 'Euclid' | 'JWST') if self.category == 'padova':
            #         return padova_formatting
            #     case _:
            #         return default_formatting
            pass
        else:
            if self.name == 'WFIRST+HST__WFC3' and self.category == 'padova':
                return oldpadova_fomatting_withlogage
            if self.name in ['WFIRST','LSST','GAIA__0','GAIA__DR2','HST__WFC3'] and self.category == 'padova':
                return oldpadova_fomatting
            if self.name in ['CTIO__DECam','LSST_DP0','Roman','Euclid','JWST__NIRCam','JWST__MIRI'] and self.category == 'padova':
                return padova_formatting
            else:
                return default_formatting

    @cached_property
    def qtables_dictionary(self):  # log10(fe[i]/0.0152
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
            print(first)
            points[first:last] = sub_points
            values[first:last] = sub_values
        return points, values


if __name__ == '__main__':
    raise NotImplementedError()
