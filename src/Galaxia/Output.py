#!/usr/bin/env python
"""
Docstring
"""
import itertools
import h5py as h5
import ebf

from .constants import *

__all__ = ['Output']

class Output:
    _export_keys = ['rad', 'teff', 'alpha', 'mtip', 'pz', 'px', 'py', 'feh', 'lum', 'mact', 'parentid', 'dmod', 'partid', 'age', 'grav', 'smass']
    def __init__(self, survey) -> None:
        self.__survey = survey

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['name']])
        return f'{cls}({description})'

    def _ebf_to_hdf5(self):
        export_keys = list(itertools.chain.from_iterable([self._export_keys]+[isochrone.to_export_keys for isochrone in self.isochrones]))
        data = ebf.read(self._ebf)
        hdf5_file = self._hdf5
        with h5.File(hdf5_file, 'w') as f5:
            for k in export_keys:
                f5.create_dataset(name=k, data=data[k])
            print(f"Exported the following quantities to {hdf5_file}")
            print(list(f5.keys()))

    def __name_with_ext(self, ext):
        name_base = self._file_base
        return name_base.parent / f"{name_base.name}{ext}"

    @property
    def survey(self):
        return self.__survey
    
    @property
    def isochrones(self):
        return self.survey.isochrones

    @property
    def name(self):
        return f"{self.survey.surveyname}.{self.survey.inputname}"

    @property
    def _file_base(self):
        return GALAXIA_TMP / self.name
    
    @property
    def _ebf(self):
        return self.__name_with_ext('.ebf')

    @property
    def _hdf5(self):
        return self.__name_with_ext('.h5')


if __name__ == '__main__':
    pass
