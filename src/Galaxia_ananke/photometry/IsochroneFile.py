#!/usr/bin/env python
"""
Docstring
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Union
from functools import cached_property
import pathlib
import re

import pandas as pd
from astropy.io import ascii

from .._constants import *

if TYPE_CHECKING:
    from .Isochrone import Isochrone

__all__ = ['IsochroneFile']


class IsochroneFile:
    _file_format = "output_{}.dat"
    def __init__(self, *args, isochrone: Isochrone = None) -> None:
        if isochrone is None:
            raise TypeError(f"Keyword argument 'isochrone' missing")
        self._isochrone = isochrone
        if not args:
            raise TypeError("Isochrone requires at least one argument")
        elif len(args) in [1, 2]:
            self._path = pathlib.Path(args[0])
            if len(args) == 2:
                self._write_table(args[1])
        else:
            raise TypeError(f"Too many arguments ({len(args)} given)")

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['filename']])
        return f'{cls}({description})'

    def _write_table(self, data: Union[pd.DataFrame, ascii.core.Table]) -> None:
        if isinstance(data, pd.DataFrame):
            data = ascii.core.Table.from_pandas(data)
        elif not isinstance(data, ascii.core.Table):
            raise ValueError("Given data should either be a pandas DataFrame or an astropy Table")
        data.sort([self.isochrone._age, self.isochrone._mini])
        data.write(self.path, format='ascii.commented_header')

    def _open(self, *args, **kwargs):
        return open(self._path, *args, **kwargs)

    @property
    def path(self):
        return self._path

    @property
    def filename(self):
        return self._path.name

    @property
    def metallicity(self):
        return float(re.findall("(?<={})(.*)(?={})".format(*tuple(self._file_format.split('{}'))),
                                self.filename)[0])

    @cached_property
    def column_names(self):
        with self._open('r') as f:
            _temp = '#'
            while (_temp[0] == '#') if _temp else False:
                _temp = f.readline()
                if (_temp[0] == '#') if _temp else False:
                    header = _temp
        header = header.strip('#').strip(' ').strip('\n').replace('\t',' ')
        while header.count('  '):
            header = header.replace('  ', ' ')
        return header.split(' ')

    @cached_property
    def data(self) -> ascii.core.Table:
        return ascii.read(self._path, names=self.column_names)
    
    @property
    def isochrone(self):
        return self._isochrone


if __name__ == '__main__':
    pass
