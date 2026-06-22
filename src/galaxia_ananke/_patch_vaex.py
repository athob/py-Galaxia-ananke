#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 (AGPL-3.0).
# 
# The full copyright notice, including terms governing use, modification,
# and redistribution, is contained in the files LICENSE and COPYRIGHT,
# which can be found at the root of the source code distribution tree:
# - LICENSE <https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE>
# - COPYRIGHT <https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT>
#
"""
"""
import h5py
import vaex.hdf5.dataset

from vaex import open, concat

__all__ = ['apply_vaex_patch']


_original_hdf5_load = vaex.hdf5.dataset.Hdf5MemoryMapped._load

def new_hdf5_load(self: vaex.hdf5.dataset.Hdf5MemoryMapped):  # https://github.com/vaexio/vaex/blob/65ab46281939e2fd2fc291266bb08a328ff59882/packages/vaex-hdf5/vaex/hdf5/dataset.py#L186-L221
    self.ucds = {}
    self.descriptions = {}
    self.units = {}

    if self.group is None:
        if "data" in self.h5file:
            self._load_columns(self.h5file["/data"])
            self.group = "/data"
        if "table" in self.h5file:
            self._version = 2
            self._load_columns(self.h5file["/table"])
            self.group = "/table"
        root_datasets = [dataset for name, dataset in self.h5file.items() if isinstance(dataset, h5py.Dataset)]
        if len(root_datasets):
            # if we have datasets at the root, we assume 'version 1'
            self._load_columns(self.h5file)
            self.group = "/"

        # TODO: shall we rename it vaex... ?
        # if "vaex" in self.h5file:
        # self.load_columns(self.h5file["/vaex"])
        # h5table_root = "/vaex"
        if "columns" in self.h5file:
            self._load_columns(self.h5file["/columns"])
            self.group = "/columns"
    else:
        self._version = 1
        self._load_columns(self.h5file[self.group])

    if "properties" in self.h5file:
        self._load_variables(self.h5file["/properties"])  # old name, kept for portability
    if "variables" in self.h5file:
        self._load_variables(self.h5file["/variables"])
    # self.update_meta()
    # self.update_virtual_meta()


_original_open_many = vaex.open_many

def new_open_many(filenames, **kwargs):  # github.com/vaexio/vaex/tree/65ab46281939e2fd2fc291266bb08a328ff59882/packages/vaex-core/vaex/__init__.py#L273-L286
    """Open a list of filenames, and return a DataFrame with all DataFrames concatenated.

    The filenames can be of any format that is supported by :py:func:`vaex.open`, namely hdf5, arrow, parquet, csv, etc.

    :param list[str] filenames: list of filenames/paths
    :rtype: DataFrame
    """
    dfs = []
    for filename in filenames:
        filename = filename.strip()
        if filename and filename[0] != "#":
            dfs.append(open(filename, **kwargs))
    return concat(dfs)


def apply_vaex_patch():
    vaex.hdf5.dataset.Hdf5MemoryMapped._load = new_hdf5_load
    vaex.open_many = new_open_many


if __name__ == '__main__':
    raise NotImplementedError()
