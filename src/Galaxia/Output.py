#!/usr/bin/env python
"""
Docstring
"""
import pathlib
import itertools
import h5py as h5
import ebf
import vaex

from .constants import *

__all__ = ['Output']

class Output:
    _export_keys = ['ra', 'dec', 'glon', 'glat', 'rad', 'teff', 'alpha', 'mtip', 'pz', 'px', 'py', 'feh', 'lum', 'mact', 'parentid', 'dmod', 'partid', 'age', 'grav', 'smass']
    def __init__(self, survey) -> None:  # TODO kwargs given to parameter file to run survey should be accessible from here: bonus TODO SkyCoord for center point: SkyCoord(u=-rSun[0], v=-rSun[1], w=-rSun[2], unit='kpc', representation_type='cartesian', frame='galactic')
        self.__survey = survey
        self.__vaex = None
        self.__path = None

    def __repr__(self) -> str:
        # cls = self.__class__.__name__
        # description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['name']])
        # return f'{cls}({description})'
        return repr(self._vaex)

    def _ebf_to_hdf5(self):
        export_keys = list(itertools.chain.from_iterable([self._export_keys]+[isochrone.to_export_keys for isochrone in self.isochrones]))
        # data = ebf.read(self._ebf)
        hdf5_file = self._hdf5
        with h5.File(hdf5_file, 'w') as f5:
            for k in export_keys:
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
    def name(self):
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
        return GALAXIA_TMP / self.name
    
    @property
    def _ebf(self):
        return self.__name_with_ext('.ebf')

    @property
    def _hdf5(self):
        return self.__name_with_ext('.h5')


if __name__ == '__main__':
    pass
