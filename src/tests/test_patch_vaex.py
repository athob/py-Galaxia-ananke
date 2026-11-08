#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2026  Adrien CR Thob
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
Tests of the patch that the package makes to vaex: opening an HDF5 file by one of its groups (how the catalogues are stored: the
group of the stars, next to those of the parent particles and of the partitioning), and opening a list of files that may hold
blank lines and comments. The files are small ones that are written here.
"""
import hashlib
import os
import pathlib

import h5py
import numpy as np
import pytest
import vaex

from ..galaxia_ananke import _patch_vaex
from ..galaxia_ananke._constants import PARENTPARTICLE_GROUP, PARTITIONING_GROUP, STARCATALOG_GROUP
from .utils import write_catalogue_hdf5


def values(df, column):
    return df[column].to_numpy().tolist()


def ananke_file(tmp_path, name='file.h5', stars=(1.0, 2.0, 3.0)):
    """A file laid out as the catalogues are: groups for the stars, the parent particles and the partitioning."""
    path = write_catalogue_hdf5(tmp_path / name, {'x': np.array(stars), 'n': np.arange(len(stars))}, attrs={'Header': 'made for a test'})
    with h5py.File(path, 'r+') as f5:
        f5.create_group(PARENTPARTICLE_GROUP).create_dataset('parentid', data=np.arange(5))
        f5.create_group(PARTITIONING_GROUP).create_dataset('partition_id', data=[0])
    return path


# -- opening a group ----------------------------------------------------------------

def test_opening_a_group_gives_the_datasets_of_that_group_and_no_other(tmp_path):
    df = vaex.open(ananke_file(tmp_path), group=STARCATALOG_GROUP)
    assert sorted(df.get_column_names()) == ['n', 'x'] and values(df, 'x') == [1.0, 2.0, 3.0] and values(df, 'n') == [0, 1, 2]
    assert len(df) == 3
    df.close()


def test_each_group_of_a_file_can_be_opened_by_itself(tmp_path):
    path = ananke_file(tmp_path)
    parents = vaex.open(path, group=PARENTPARTICLE_GROUP)
    assert parents.get_column_names() == ['parentid'] and len(parents) == 5
    parents.close()


def test_the_dtypes_of_the_datasets_are_kept(tmp_path):
    df = vaex.open(ananke_file(tmp_path), group=STARCATALOG_GROUP)
    assert df['x'].to_numpy().dtype == np.float64 and df['n'].to_numpy().dtype == np.int64
    df.close()


@pytest.mark.parametrize("layout", ["data", "columns"])
def test_the_groups_that_vaex_looks_for_by_themselves_are_still_found(layout, tmp_path):
    path = write_catalogue_hdf5(tmp_path / 'file.h5', {'a': np.array([1.0, 2.0])}, group=layout)
    df = vaex.open(path)
    assert df.get_column_names() == ['a'] and values(df, 'a') == [1.0, 2.0]
    df.close()


def test_datasets_at_the_root_are_found_too(tmp_path):
    path = tmp_path / 'file.h5'
    with h5py.File(path, 'w') as f5:
        f5.create_dataset('a', data=np.array([4.0, 5.0, 6.0]))
    df = vaex.open(path)
    assert df.get_column_names() == ['a'] and values(df, 'a') == [4.0, 5.0, 6.0]
    df.close()


def test_a_file_that_cannot_be_written_to_is_opened_and_left_as_it_is(tmp_path):
    path = ananke_file(tmp_path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    path.chmod(0o444)
    try:
        df = vaex.open(path, group=STARCATALOG_GROUP)
        assert values(df, 'x') == [1.0, 2.0, 3.0]
        df.close()
    finally:
        path.chmod(0o644)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


def test_opening_a_file_does_not_change_it(tmp_path):
    path = ananke_file(tmp_path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    vaex.open(path, group=STARCATALOG_GROUP).close()
    vaex.open(path, group=PARENTPARTICLE_GROUP).close()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


# -- opening a list of files --------------------------------------------------------

def test_the_files_of_a_list_are_one_table_in_the_order_of_the_list(tmp_path):
    first, second = ananke_file(tmp_path, 'a.h5', stars=(1.0, 2.0)), ananke_file(tmp_path, 'b.h5', stars=(3.0,))
    df = _patch_vaex.new_open_many([str(first), str(second)], group=STARCATALOG_GROUP)
    assert values(df, 'x') == [1.0, 2.0, 3.0] and len(df) == 3


def test_blank_lines_comments_and_the_spaces_around_the_names_are_left_out_of_the_list(tmp_path):
    first, second = ananke_file(tmp_path, 'a.h5', stars=(1.0, 2.0)), ananke_file(tmp_path, 'b.h5', stars=(3.0,))
    df = _patch_vaex.new_open_many([f"  {first}  \n", "", "   ", "# a comment", "  # another one", f"\t{second}\n"], group=STARCATALOG_GROUP)
    assert values(df, 'x') == [1.0, 2.0, 3.0]


def test_a_list_can_be_a_generator_of_names(tmp_path):
    paths = [ananke_file(tmp_path, f"{name}.h5", stars=(float(i),)) for i, name in enumerate('abc')]
    df = _patch_vaex.new_open_many(map(str, paths), group=STARCATALOG_GROUP)
    assert values(df, 'x') == [0.0, 1.0, 2.0]


def test_the_arguments_of_opening_are_given_to_each_file(tmp_path):
    first, second = ananke_file(tmp_path, 'a.h5'), ananke_file(tmp_path, 'b.h5')
    df = _patch_vaex.new_open_many([str(first), str(second)], group=PARENTPARTICLE_GROUP)
    assert df.get_column_names() == ['parentid'] and len(df) == 10


def test_the_patched_function_is_the_one_that_vaex_has_and_the_original_is_kept():
    assert vaex.open_many is _patch_vaex.new_open_many
    assert _patch_vaex._original_open_many is not _patch_vaex.new_open_many
    assert _patch_vaex._original_hdf5_load is not _patch_vaex.new_hdf5_load


if __name__ == '__main__':
    pass
