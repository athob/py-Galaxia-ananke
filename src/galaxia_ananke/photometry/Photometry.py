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
The registry of the photometric systems that are available, by category and then by name.
"""
import shutil
from glob import glob
from typing import Any, ClassVar, cast

import pandas as pd

from .._constants import *
from ..utils import Singleton
from .Isochrone import Isochrone
from .PhotoSystem import PhotoSystem

__all__ = ['Photometry']


def _set_with_slashes(dictionary: dict, key: str, value: Any) -> None:
    """
    Set a value where a name with slashes leads to in a dictionary, making the dictionaries on the way.

    The dictionary that is there is the one that is changed, and not a copy of it.

    Parameters
    ----------
    dictionary : dict
        The dictionary, which is changed.
    key : str
        The name, whose parts between slashes are the keys to go down by.
    value : object
        What is set at the end.
    """
    head, *rest = key.split('/', 1)
    if not rest:
        dict.__setitem__(dictionary, head, value)
    else:
        if head not in dictionary:
            dict.__setitem__(dictionary, head, {})
        _set_with_slashes(dict.__getitem__(dictionary, head), rest[0], value)


class nested_dict(dict):
    """Dictionary of dictionaries whose items are got and set by a name that has slashes: ``d['a/b']`` is ``d['a']['b']``."""
    def __getitem__(self, __key: str) -> Any:
        """
        Get an item by a name that may have slashes, which go down from dictionary to dictionary.

        Parameters
        ----------
        __key : str
            The name.

        Returns
        -------
        item : object
            What the name leads to.

        Raises
        ------
        KeyError
            If a part of the name is not in the dictionary that it is looked for in.
        """
        _key = __key.split('/')
        _temp = super().__getitem__(_key[0])
        return _temp if len(_key) == 1 else nested_dict(_temp)['/'.join(_key[1:])]
    
    def __setitem__(self, __key: str, __value: Any) -> None:
        """
        Set an item by a name that may have slashes, making the dictionaries that it goes down by if they are not there.

        Parameters
        ----------
        __key : str
            The name.
        __value : object
            The item.
        """
        _key = __key.split('/', 1)
        if len(_key) == 1:
            return super().__setitem__(_key[0], __value)
        else:
            if _key[0] not in self.keys():
                super().__setitem__(_key[0], {})
            return _set_with_slashes(super().__getitem__(_key[0]), _key[1], __value)  # in the dictionary that is there, not in a copy of it
    
    # def __delitem__(self, __key):
    #     return super().__delitem__(__key)


class Photometry(nested_dict, metaclass=Singleton):
    """
    The photometric systems that are available, by category and then by name: ``available_photo_systems['category/name']``.

    There is one, made when it is first asked for, from the isochrones of the installation: each folder of a category that
    has the descriptor of the files of an isochrone (the category of the bolometric corrections, ``BolCorr``, is left out).
    The category of the custom systems always exists, and the one of the legacy systems is protected: nothing can be set in
    it.
    """
    _category = 'category'
    _isochrone = 'isochrone'
    instance: ClassVar['Photometry']  # the only one, made by __new__
    def __new__(cls) -> 'Photometry':
        """
        Make the one registry of photometric systems from the isochrones of the installation.

        Returns
        -------
        photometry : Photometry
            The registry.
        """
        cls.instance = super(Photometry, cls).__new__(cls)
        _temp = pd.DataFrame([[iso_path.parent.name, Isochrone(iso_path)]
                               for iso_path in ISOCHRONES_PATH.glob('*/*')
                               if iso_path.is_dir() and iso_path.parent.name != 'BolCorr'],
                             columns=[cls._category, cls._isochrone])
        isochrones_per_category = {key: _temp[cls._isochrone].loc[item].to_list()
                                   for key, item in _temp.groupby(cls._category).groups.items()}
        for key, item in isochrones_per_category.items():
            cls.instance[cast(str, key)] = {iso.name: PhotoSystem(iso)
                                            for iso in item if iso.has_file_descriptor}
        if not CUSTOM_PHOTOCAT in cls.instance.keys(): cls.instance[CUSTOM_PHOTOCAT] = {}
        return cls.instance
    
    def __setitem__(self, __key: str, __value: Any) -> None:
        """
        Set an item by a name that may have slashes, unless it is in the protected category of the legacy systems.

        Parameters
        ----------
        __key : str
            The name.
        __value : object
            The item.

        Raises
        ------
        ValueError
            If the first part of the name is the legacy category, and it is in the registry.
        """
        __key_split = __key.split('/')
        if __key_split[0] == LEGACY_PHOTOCAT and LEGACY_PHOTOCAT in self.keys():
            raise ValueError(f"The '{LEGACY_PHOTOCAT}' entry is protected")
        # elif len(__key_split) == 1:
        #         __key = f"{CUSTOM_CATEGORY}/{__key}"
        return super().__setitem__(__key, __value)

    def __delitem__(self, __key: str) -> None:
        """
        Refuse to delete an item.

        Parameters
        ----------
        __key : str
            The name.

        Raises
        ------
        NotImplementedError
            Always: the items of the registry are not deleted.
        """
        raise NotImplementedError()
        # return super().__delitem__(__key)
        
    def add_isochrone(self, name: str, isochrone_data: dict, **kwargs: Any) -> None:
        """
        Add a custom photometric system, made from the data of its isochrone.

        The isochrone is written in the folder of the custom systems of the installation, and its photometric system is
        added to the registry in the custom category.

        Parameters
        ----------
        name : str
            The name of the system, which has no slash.
        isochrone_data : dict
            The tables of the isochrone, by metallicity (see ``Isochrone``).
        **kwargs
            Given to ``Isochrone``, such as ``overwrite`` and ``magnitude_system``.

        Raises
        ------
        ValueError
            If the name has a slash.
        FileExistsError
            If a system has this name and ``overwrite`` is not given.
        """
        if '/' in name:
            raise ValueError()
        self[CUSTOM_PHOTOCAT][name] = PhotoSystem(Isochrone(name, isochrone_data, **kwargs))

    def _purge(self, hard: bool = False, yes: bool = False) -> None:
        """
        Remove all the isochrones that are cached in the installation.

        What is in the folder of the isochrones is removed, whatever it is: a folder, a file, or a link, which is removed
        and not what it points to. Nothing is removed unless ``hard`` is set, and the user is asked first unless ``yes`` is.

        Parameters
        ----------
        hard : bool, optional
            Whether to remove. Default to False, which does nothing.
        yes : bool, optional
            Whether to remove without asking. Default to False.
        """
        if hard and (True if yes else (input("WARNING!! This will remove all the isochrones cached in your installation. Are you certain (y/[n])? ").lower() in ['y', 'yes'])):
            for path in ISOCHRONES_PATH.glob('*'):
                if path.is_dir() and not path.is_symlink():
                    shutil.rmtree(path)
                else:  # a file, or a link, which is removed and not what it points to (rmtree refuses a link)
                    path.unlink()


if __name__ == '__main__':
    raise NotImplementedError()
