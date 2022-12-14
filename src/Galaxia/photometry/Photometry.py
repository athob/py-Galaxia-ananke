#!/usr/bin/env python
"""
Docstring
"""
import pathlib
from glob import glob

import pandas as pd

from ..constants import *
from ..utils import Singleton
from .Isochrone import Isochrone

__all__ = ['Photometry']


class nested_dict(dict):
    def __getitem__(self, __key):
        _key = __key.split('/')
        _temp = super().__getitem__(_key[0])
        return _temp if len(_key) == 1 else nested_dict(_temp)['/'.join(_key[1:])]
    
    def __setitem__(self, __key, __value):
        _key = __key.split('/')
        if len(_key) == 1:
            return super().__setitem__(_key[0], __value)
        else:
            if _key[0] not in self.keys():
                super().__setitem__(_key[0], {})
            return nested_dict(super().__getitem__(_key[0])).__setitem__('/'.join(_key[1:]), __value)
    
    # def __delitem__(self, __key):
    #     return super().__delitem__(__key)


class Photometry(nested_dict, metaclass=Singleton):
    def __new__(cls):
        cls.instance = super(Photometry, cls).__new__(cls)
        _temp = pd.DataFrame([[iso_path.parent.name, Isochrone(iso_path)]
                               for iso_path in map(pathlib.Path, glob(str(ISOCHRONES_PATH / '*' / '*')))
                               if iso_path.is_dir() and iso_path.parent.name != 'BolCorr'])
        _temp = {key: _temp[1].loc[item].to_list()
                 for key, item in _temp.groupby(0).groups.items()}
        for key, item in _temp.items():
            cls.instance[key] = {iso.name: iso
                                 for iso in item}
        if not CUSTOM_PHOTOCAT in cls.instance.keys(): cls.instance[CUSTOM_PHOTOCAT] = {}
        return cls.instance
    
    def __setitem__(self, __key, __value) -> None:
        __key_split = __key.split('/')
        if __key_split[0] is LEGACY_PHOTOCAT and LEGACY_PHOTOCAT in self.keys():
            raise ValueError(f"The '{LEGACY_PHOTOCAT}' entry is protected")
        # elif len(__key_split) == 1:
        #         __key = f"{CUSTOM_CATEGORY}/{__key}"
        return super().__setitem__(__key, __value)

    def __delitem__(self, __key):
        raise NotImplementedError()
        # return super().__delitem__(__key)
        
    def add_isochrone(self, name, isochrone_data, **kwargs):
        if '/' in name:
            raise ValueError()
        self[CUSTOM_PHOTOCAT][name] = Isochrone(name, isochrone_data, **kwargs)


if __name__ == '__main__':
    pass
