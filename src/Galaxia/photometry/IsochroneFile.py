#!/usr/bin/env python
"""
Docstring
"""
import pathlib
import re

import pandas as pd
from astropy.io import ascii

from ..constants import *

__all__ = ['IsochroneFile']


class IsochroneFile:
    def __init__(self, *args, isochrone=None) -> None:
        if not args:
            raise TypeError("Isochrone requires at least one argument")
        elif len(args) in [1, 2]:
            self._path = pathlib.Path(args[0])
            if len(args) == 2:
                self._write_table(args[1])
        else:
            raise TypeError(f"Too many arguments ({len(args)} given)")
        if isochrone is None:
            raise TypeError(f"Keyword argument 'isochrone' missing")
        self._isochrone = isochrone
        self._column_names = None
        self._data = None

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['filename']])
        return f'{cls}({description})'

    def _write_table(self, data):
        if isinstance(data, pd.DataFrame):
            data = ascii.core.Table.from_pandas(data)
        elif not isinstance(data, ascii.core.Table):
            raise ValueError("Given data should either be a pandas DataFrame or an astropy Table")
        data.sort_values(by=[self.isochrone._Age, self.isochrone._M_ini]).write(self.path, format='ascii.commented_header')

    def _open(self, *args, **kwargs):
        return open(self._path, *args, **kwargs)
    
    def _load_column_names(self):
        with self._open('r') as f:
            _temp = '#'
            while (_temp[0] is '#') if _temp else False:
                _temp = f.readline()
                if (_temp[0] is '#') if _temp else False:
                    header = _temp
        header = header.strip('#').strip(' ').strip('\n').replace('\t',' ')
        while header.count('  '):
            header = header.replace('  ', ' ')
        self._column_names = header.split(' ')
    
    def _load_data(self):
        self._data = ascii.read(self._path, names=self.column_names)


    @property
    def path(self):
        return self._path

    @property
    def filename(self):
        return self._path.name

    @property
    def metallicity(self):
        return float(re.findall("(?<={})(.*)(?={})".format(*tuple(ISO_FILE_FORMAT.split('{}'))),
                                self.filename)[0])

    @property
    def column_names(self):
        if self._column_names is None:
            self._load_column_names()
        return self._column_names

    @property
    def data(self):
        if self._data is None:
            self._load_data()
        return self._data
    
    @property
    def isochrone(self):
        return self._isochrone
