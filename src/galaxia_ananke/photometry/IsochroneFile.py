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
from __future__ import annotations
from typing import TYPE_CHECKING, Any, IO, List, Optional, Union
from functools import cached_property
import pathlib
import re

import numpy as np
import pandas as pd
from astropy.io import ascii

from .._constants import *

if TYPE_CHECKING:
    from .Isochrone import Isochrone

__all__ = ['IsochroneFile']


class IsochroneFile:
    """
    The file of the isochrone of a photometric system for one metallicity.

    The name of the file is ``output_`` followed by the metallicity and ``.dat``, and the names of its columns are in a
    commented line before the table. A file is made from the path of one that exists, or from a table that is written there.
    It is a file of an isochrone, which it has to be given.
    """
    _file_format = "output_{}.dat"
    def __init__(self, *args: Any, isochrone: Optional[Isochrone] = None) -> None:
        """
        Make the file from its path, and write it if a table is given.

        Parameters
        ----------
        *args
            The path of the file, and optionally the table to write in it (see ``_write_table``).
        isochrone : Isochrone
            The isochrone that this is a file of. It has to be given.

        Raises
        ------
        TypeError
            If the isochrone is not given, or no argument, or more than two.
        ValueError
            If the table is neither a pandas DataFrame nor an astropy Table.
        """
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
        """
        Give the representation of the file.

        Returns
        -------
        representation : str
            The class, with the name of the file.
        """
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['filename']])
        return f'{cls}({description})'

    def _write_table(self, data: Union[pd.DataFrame, ascii.core.Table]) -> None:
        """
        Write a table in the file, sorted by age and then by initial mass.

        The sort keeps the order of the rows that share both, which the sort of an astropy table does not, and it does not
        change the table that was given. The columns are written with the names in a commented line, as the backend reads
        them.

        Parameters
        ----------
        data : pandas.DataFrame or astropy.table.Table
            The table, which has the columns of the age and of the initial mass of the isochrone.

        Raises
        ------
        ValueError
            If it is neither of the two.
        """
        if isinstance(data, pd.DataFrame):
            data = ascii.core.Table.from_pandas(data)
        elif not isinstance(data, ascii.core.Table):
            raise ValueError("Given data should either be a pandas DataFrame or an astropy Table")
        data = data[np.lexsort((data[self.isochrone._mini], data[self.isochrone._age]))]  # age first, then mass; stable, unlike Table.sort, and not in place
        data.write(self.path, format='ascii.commented_header')

    def _open(self, *args: Any, **kwargs: Any) -> IO[Any]:
        """
        Open the file as ``open`` does.

        Parameters
        ----------
        *args, **kwargs
            Given to ``open``, after the path of the file.

        Returns
        -------
        file : file object
            The file.
        """
        return open(self._path, *args, **kwargs)

    @property
    def path(self) -> pathlib.Path:
        """The path of the file."""
        return self._path

    @property
    def filename(self) -> str:
        """The name of the file."""
        return self._path.name

    @property
    def metallicity(self) -> float:
        """The metallicity of the file: the number in its name, between ``output_`` and ``.dat``."""
        return float(re.findall("(?<={})(.*)(?={})".format(*tuple(self._file_format.split('{}'))),
                                self.filename)[0])

    @cached_property
    def column_names(self) -> List[str]:
        """
        The names of the columns of the file: the words of the last commented line before the table.

        Whatever white space separates them, and wherever it is, is taken for a separator.
        """
        with self._open('r') as f:
            _temp = '#'
            while (_temp[0] == '#') if _temp else False:
                _temp = f.readline()
                if (_temp[0] == '#') if _temp else False:
                    header = _temp
        return header.strip('#').split()  # whatever the white space is, and wherever: tabs, spaces, the newline that ends the line

    @cached_property
    def data(self) -> ascii.core.Table:
        """The table of the file, as astropy reads it, with the names of the columns that its header gives."""
        return ascii.read(self._path, names=self.column_names)
    
    @property
    def isochrone(self) -> Isochrone:
        """The isochrone that this is a file of."""
        return self._isochrone


if __name__ == '__main__':
    raise NotImplementedError()
