#!/usr/bin/env python
"""
Docstring
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import pathlib
import itertools
import h5py as h5
import ebf
import vaex

from .constants import *

if TYPE_CHECKING:
    from . import Survey

__all__ = ['Output']


class Output:
    _pos = ['px', 'py', 'pz']
    _vel = ['vx', 'vy', 'vz']
    _cel = ['ra', 'dec']
    _gal = ['glon', 'glat']
    _rad = 'rad'
    _dmod = 'dmod'
    _export_keys = _pos + _vel + _cel + _gal + [_rad, _dmod, 'teff', 'mtip', 'feh', 'lum', 'mact', 'parentid', 'partid', 'age', 'grav', 'smass']
    _vaex_under_list = ['_repr_html_']
    def __init__(self, survey: Survey) -> None:  # TODO kwargs given to parameter file to run survey should be accessible from here: bonus TODO SkyCoord for center point: SkyCoord(u=-rSun[0], v=-rSun[1], w=-rSun[2], unit='kpc', representation_type='cartesian', frame='galactic')
        self.__survey = survey
        self.__vaex = None
        self.__path = None

    def __dir__(self):
        return sorted({i for i in self.__vaex.__dir__() if not i.startswith('_')}.union(
            super(Output, self).__dir__()).union(
            self._vaex_under_list if self.__vaex is not None else []))

    def __repr__(self) -> str:
        return repr(self._vaex)
    
    def __getitem__(self, item):
        return self._vaex[item]
    
    def __setitem__(self, item, value):
        self._vaex[item] = value
    
    def __getattr__(self, item):
        if (item in self.__vaex.__dir__() and not item.startswith('_')) or (item in self._vaex_under_list and self.__vaex is not None):
            return getattr(self.__vaex, item)
        else:
            return self.__getattribute__(item)

    @classmethod
    def _compile_export_mag_names(cls, isochrones):
        return list(itertools.chain.from_iterable([isochrone.to_export_keys for isochrone in isochrones]))
    
    @classmethod
    def _make_export_keys(cls, isochrones, extra_keys=[]):
        return cls._export_keys + extra_keys + cls._compile_export_mag_names(isochrones)

    def _ebf_to_hdf5(self):
        hdf5_file = self._hdf5
        with h5.File(hdf5_file, 'w') as f5:
            for k in self.export_keys:
                f5.create_dataset(name=k, data=ebf.read(self._ebf, f"/{k}"))
            print(f"Exported the following quantities to {hdf5_file}")
            print(list(f5.keys()))
        self.__vaex = vaex.open(hdf5_file)

    def __name_with_ext(self, ext):
        name_base = self._file_base
        return name_base.parent / f"{name_base.name}{ext}"
    
    def save(self, path):
        """
            Save output to new path
        """
        old_path = self._path
        self.__path = pathlib.Path(path)
        self._vaex.close()
        old_path.rename(self._path)
        self.__vaex = vaex.open(self._path)

    @property
    def survey(self):
        return self.__survey
    
    @property
    def isochrones(self):
        return self.survey.isochrones

    @property
    def export_keys(self):
        return self._make_export_keys(self.isochrones, extra_keys=[k if k != 'id' else 'satid' for k in self.survey.input.optional_keys()])

    @property
    def output_name(self):
        return f"{self.survey.surveyname}.{self.survey.inputname}"

    @property
    def _vaex(self):
        if self.__vaex is None:
            raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        else:
            return self.__vaex

    @property
    def _path(self):
        if self.__path is None:
            return self._hdf5
        else:
            return self.__path

    @property
    def _file_base(self):
        return GALAXIA_TMP / self.output_name
    
    @property
    def _ebf(self):
        return self.__name_with_ext('.ebf')

    @property
    def _hdf5(self):
        return self.__name_with_ext('.h5')
    
    def flush_extra_columns_to_hdf5(self, with_columns=[]):  # temporary until vaex supports it
        hdf5_file = self._hdf5
        old_column_names = set(vaex.open(hdf5_file).column_names)
        with h5.File(hdf5_file, 'r+') as f5:
            extra_columns = [k for k in set(self.column_names)-old_column_names if not k.startswith('__')]
            for k in extra_columns:
                f5.create_dataset(name=k, data=self[k].to_numpy())
            if extra_columns:
                print(f"Exported the following quantities to {hdf5_file}")
                print(extra_columns)
            for k in with_columns:
                f5[k][...] = self[k].to_numpy()
            if with_columns:
                print(f"Overwritten the following quantities to {hdf5_file}")
                print(with_columns)
        self.__vaex = vaex.open(hdf5_file)


if __name__ == '__main__':
    pass
