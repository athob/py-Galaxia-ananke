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
Monkey-patches and adjustments for the vaex backend.

Currently includes:
- Modified HDF5 loading logic: when a group is given, its datasets are the columns (what vaex calls version 1), as they are in
  the files of the catalogues. vaex 4.19.0 and later take a group to hold its columns in a subgroup ``columns`` (version 2),
  and fail without this.
- Replacement of ``vaex.open_many``: it gives its keyword arguments, such as ``group``, to ``vaex.open`` for each file, which
  the function of vaex, whose only argument is the list of files, does not take. (vaex already leaves the comment and the blank
  lines of the list out, as this one does.)
- Wrapping of ``vaex.utils.create_hasher``: vaex hashes with one thread in a process that was made by ``fork``. It hashes
  large data with the thread pool of blake3 (``max_threads=blake3.AUTO``), which a forked child inherits without its threads,
  so that a child that hashes after its parent did waits forever (reported to vaex, https://github.com/vaexio/vaex/issues/2503,
  and to blake3-py, https://github.com/oconnor663/blake3-py/issues/127). The parent keeps its threads, and a digest does not
  depend on their number.
- Wrapping of ``Hdf5MemoryMapped._open``: the file object that vaex gives to h5py is closed when h5py closes the HDF5 file.
  vaex opens the file twice, once for h5py and once for the memory mapping, which it stores under the same key of ``file_map``:
  the file object of h5py is no longer known to the dataset, and is finalized unclosed when the end of ``__init__`` closes the
  HDF5 file, which makes a ``ResourceWarning`` at every ``vaex.open``. Nothing leaks (its descriptor is released then).
- Suppression of all vaex internal logging (which is excessively verbose for a
  backend library and can show non-critical errors from uninstalled plugins).
"""
from __future__ import annotations
from typing import Any, Iterable
import logging
import os
from warnings import warn

import h5py
import vaex.file
import vaex.hdf5.dataset
import vaex.utils

from vaex import open, concat

from ._constants import NAME

__all__ = ['apply_vaex_patch']


_original_hdf5_load = vaex.hdf5.dataset.Hdf5MemoryMapped._load

def new_hdf5_load(self: vaex.hdf5.dataset.Hdf5MemoryMapped) -> None:  # https://github.com/vaexio/vaex/blob/65ab46281939e2fd2fc291266bb08a328ff59882/packages/vaex-hdf5/vaex/hdf5/dataset.py#L186-L221
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
        self._version = 1  # the columns are in the group itself (vaex 4.19.0 and later set 2 here, for a subgroup 'columns')
        self._load_columns(self.h5file[self.group])

    if "properties" in self.h5file:
        self._load_variables(self.h5file["/properties"])  # old name, kept for portability
    if "variables" in self.h5file:
        self._load_variables(self.h5file["/variables"])
    # self.update_meta()
    # self.update_virtual_meta()


_original_open_many = vaex.open_many

def new_open_many(filenames: Iterable[str], **kwargs: Any) -> vaex.dataframe.DataFrame:  # github.com/vaexio/vaex/tree/65ab46281939e2fd2fc291266bb08a328ff59882/packages/vaex-core/vaex/__init__.py#L273-L286
    """Open a list of filenames, and return a DataFrame with all DataFrames concatenated.

    The filenames can be of any format that is supported by :py:func:`vaex.open`, namely hdf5, arrow, parquet, csv, etc.

    :param list[str] filenames: list of filenames/paths
    :param kwargs: given to :py:func:`vaex.open` for each file (for example ``group``), which vaex's own function does not do
    :rtype: DataFrame
    """
    dfs = []
    for filename in filenames:
        filename = filename.strip()
        if filename and filename[0] != "#":
            dfs.append(open(filename, **kwargs))
    return concat(dfs)


_original_hdf5_open = vaex.hdf5.dataset.Hdf5MemoryMapped._open

def new_hdf5_open(self: vaex.hdf5.dataset.Hdf5MemoryMapped, path: Any) -> None:
    """
    Open the file as vaex does, and close the file object that h5py was given when h5py closes the HDF5 file.

    A file object of the caller, and a vaex whose ``_open`` no longer leaves a file object in ``file_map``, are left as vaex
    has them.
    """
    _original_hdf5_open(self, path)
    file = getattr(self, 'file_map', {}).get(self.path)
    if file is None or not hasattr(self, 'h5file') or vaex.file.is_file_object(path):
        return
    close_hdf5 = self.h5file.close

    def close_hdf5_and_its_file() -> None:
        close_hdf5()
        file.close()
    self.h5file.close = close_hdf5_and_its_file


_original_create_hasher = vaex.utils.create_hasher

_in_forked_child = False

def _mark_forked_child() -> None:
    global _in_forked_child
    _in_forked_child = True

def new_create_hasher(initial_data: Any = None, large_data: bool = True) -> Any:
    """Create a blake3 hasher as vaex does, with one thread (``large_data=False``) in a process that was made by ``fork``."""
    return _original_create_hasher(initial_data, large_data=large_data and not _in_forked_child)


def apply_vaex_patch() -> None:
    """
    Apply all vaex-related customisations:

    1. Suppress all vaex internal logging (WARNING, ERROR, etc.) so that
       non-critical messages (e.g., plugin loading failures, memmap cleanup
       warnings) do not pollute the user's output.
    2. Replace the HDF5 loading logic with a version that takes the datasets of a given group for the columns, as the files of
       the catalogues have them.
    3. Replace ``vaex.open_many`` with a version that gives its keyword arguments (``group``) to ``vaex.open``.
    4. Wrap ``vaex.utils.create_hasher`` so that a process made by ``fork`` hashes with one thread, and register what tells it.
    5. Wrap ``Hdf5MemoryMapped._open`` so that the file object that h5py is given is closed when h5py closes the HDF5 file.
    """
    # --- Silence vaex logging ---
    logging.getLogger("vaex").setLevel(logging.CRITICAL)

    # --- Apply HDF5 patch ---
    warn(
        f"Importing {NAME} monkey-patches vaex's HDF5 reading, "
        f"which may cause unexpected behavior in vaex's other HDF5 I/O operations.",
        RuntimeWarning,
        stacklevel=3
    )
    vaex.hdf5.dataset.Hdf5MemoryMapped._load = new_hdf5_load

    # --- Replace open_many ---
    vaex.open_many = new_open_many

    # --- Hash with one thread in a forked process ---
    vaex.utils.create_hasher = new_create_hasher
    if hasattr(os, "register_at_fork"):
        os.register_at_fork(after_in_child=_mark_forked_child)

    # --- Close the file object of h5py ---
    vaex.hdf5.dataset.Hdf5MemoryMapped._open = new_hdf5_open


if __name__ == '__main__':
    raise NotImplementedError()
