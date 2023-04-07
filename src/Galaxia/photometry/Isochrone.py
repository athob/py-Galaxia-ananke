#!/usr/bin/env python
"""
Docstring
"""
import re
import pathlib
import shutil
from glob import glob

from ..constants import *
from . import Photometry
from .IsochroneFile import IsochroneFile

__all__ = ['Isochrone']


class Isochrone:
    """
    Ages need to be in log scale
    """
    _Age = 'Age'
    _M_ini = 'M_ini'
    _M_act = 'M_act'
    _Lum = 'Lum'
    _T_eff = 'T_eff'
    _Grav = 'Grav'
    _required_keys = [_Age, _M_ini, _M_act, _Lum, _T_eff, _Grav]
    _required_metallicities = {0.0001, 0.0002, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.0012, 0.0014, 0.0016, 0.0018, 0.002, 0.0022, 0.0024, 0.0026, 0.003, 0.0034, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.024, 0.028, 0.03}
    _required_cmd_magnames_dictkeys = {'magnitude', 'color_minuend', 'color_subtrahend'}
    _file_descriptor = "IsoFileDescriptor.txt"
    def __init__(self, *args, overwrite=False) -> None:
        if not args:
            raise TypeError("Isochrone requires at least one argument")
        elif len(args) == 1:
            self._path = pathlib.Path(args[0])
        elif len(args) == 2:
            self._path = Photometry.ISOCHRONES_PATH / CUSTOM_PHOTOCAT / args[0]
            if self.path.exists:
                if overwrite:
                    shutil.rmtree(self.path)
                else:
                    raise FileExistsError(f"Isochrone '{args[0]}' already exists at '{self.path}' (use `overwrite` kwarg to ignore)")
            self._write_isochrone_files(args[1])
        else:
            raise TypeError(f"Too many arguments ({len(args)} given)")
        self._isochrone_files = None
    
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
        metallicities = set(metallicities)
        if metallicities != self._required_metallicities:
            missing = self._required_metallicities.difference(metallicities)
            missing = f"misses {missing}" if missing else ""
            extra = metallicities.difference(self._required_metallicities)
            extra = f"misincludes {extra}" if extra else ""
            raise ValueError(f"Given isochrone data covers wrong set of metallicities: {missing}{' & ' if missing and extra else ''}{extra}")
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

    def _prepare_isochrone_files(self):
        self._isochrone_files = list(map(lambda path: IsochroneFile(path, isochrone=self),
                                         glob(str(self._path / IsochroneFile._file_format.format('*')))))
    
    def check_cmd_magnames(self, cmd_magnames):
        if isinstance(cmd_magnames, str):
            check = set(re.split('[ ,-]', cmd_magnames))
        elif isinstance(cmd_magnames, dict):
            cmd_magnames_keys = set(cmd_magnames.keys())
            if cmd_magnames_keys != self._required_cmd_magnames_dictkeys:
                missing = self._required_cmd_magnames_dictkeys.difference(cmd_magnames_keys)
                missing = f"misses {missing}" if missing else ""
                extra = cmd_magnames_keys.difference(self._required_cmd_magnames_dictkeys)
                extra = f"misincludes {extra}" if extra else ""
                raise ValueError(f"Given cmd_magnames dict covers wrong set of keys: {missing}{' & ' if missing and extra else ''}{extra}")
            check = set(cmd_magnames.values())
            cmd_magnames = f"{cmd_magnames['magnitude']},{cmd_magnames['color_minuend']}-{cmd_magnames['color_subtrahend']}"
        else:
            raise ValueError(f"cmd_magnames should be either a string or a dict: given {type(cmd_magnames)} instead")
        if not check.issubset(set(self.mag_names)):
            raise ValueError(f"CMD magnitude names '{cmd_magnames}' don't match isochrone names: choose among {self.mag_names}")
        else:
            return cmd_magnames

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
    def key(self):
        return f"{self.category}/{self.name}"

    @property
    def file_descriptor_path(self):
        return self._path / self._file_descriptor

    @property
    def has_file_descriptor(self):
        return self.file_descriptor_path.exists()

    @property
    def isochrone_files(self):
        if self._isochrone_files is None:
            self._prepare_isochrone_files()
        return self._isochrone_files
    
    @property
    def mag_names(self):
        if self.has_file_descriptor:
            with open(self.file_descriptor_path,'r') as f: return f.readline().strip('\n').split()[4:]
        else:
            raise NotImplementedError("Photometric system doesn't have an IsoFileDescriptor.txt file")

    @property
    def to_export_keys(self):
        return [f"{self.name.lower()}_{mag_name.lower()}" for mag_name in self.mag_names]


if __name__ == '__main__':
    pass
