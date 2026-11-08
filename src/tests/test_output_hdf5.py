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
Tests of the two helpers that make a post-processing pipeline persistent: the one that writes the columns that a pipeline made
into the HDF5 file of the catalogue (and overwrites those that it changed, and updates the metadata), and the decorator that applies
a pipeline to a table, or opens the file of a catalogue, applies it and flushes the result. The catalogues are small HDF5 files.
"""
import h5py
import numpy as np
import pytest
import vaex

from ..galaxia_ananke._constants import STARCATALOG_GROUP
from ..galaxia_ananke.Output import _decorate_post_processing, _flush_extra_columns_to_hdf5
from .utils import write_catalogue_hdf5


def catalogue(tmp_path, name='catalogue.h5', **attrs):
    return write_catalogue_hdf5(tmp_path / name, {'x': np.array([1.0, 2.0, 3.0]), 'y': np.array([10, 20, 30])}, attrs=attrs)


def datasets(path):
    with h5py.File(path, 'r') as f5:
        return {name: f5[STARCATALOG_GROUP][name][...].tolist() for name in f5[STARCATALOG_GROUP]}


def attributes(path):
    with h5py.File(path, 'r') as f5:
        return dict(f5.attrs)


# -- flushing the columns that were made --------------------------------------------

def test_a_column_that_was_made_is_written_in_the_file_next_to_the_others(tmp_path):
    path = catalogue(tmp_path)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    df['z'] = df.x * 2
    _flush_extra_columns_to_hdf5(df, path, verbose=False)
    df.close()
    assert datasets(path) == {'x': [1.0, 2.0, 3.0], 'y': [10, 20, 30], 'z': [2.0, 4.0, 6.0]}


def test_the_columns_that_are_not_new_are_left_as_they_are_unless_they_are_asked_for(tmp_path):
    path = catalogue(tmp_path)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    df['x'] = df.x + 100
    _flush_extra_columns_to_hdf5(df, path, verbose=False)
    assert datasets(path)['x'] == [1.0, 2.0, 3.0]
    _flush_extra_columns_to_hdf5(df, path, with_columns=['x'], verbose=False)
    df.close()
    assert datasets(path)['x'] == [101.0, 102.0, 103.0] and datasets(path)['y'] == [10, 20, 30]


def test_a_column_that_is_new_is_not_overwritten_but_made_once_when_it_is_also_asked_for(tmp_path, capsys):
    path = catalogue(tmp_path)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    df['z'] = df.x * 2
    _flush_extra_columns_to_hdf5(df, path, with_columns=['z'], verbose=True)
    df.close()
    assert datasets(path)['z'] == [2.0, 4.0, 6.0]
    out = capsys.readouterr().out
    assert "Exported the following quantities" in out and "Overwritten" not in out


def test_the_columns_that_vaex_makes_for_itself_are_not_written(tmp_path):
    path = catalogue(tmp_path)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    df['__temporary'] = df.x * 3
    df['z'] = df.x * 2
    _flush_extra_columns_to_hdf5(df, path, verbose=False)
    df.close()
    assert set(datasets(path)) == {'x', 'y', 'z'}


def test_flushing_without_anything_new_changes_nothing(tmp_path):
    path = catalogue(tmp_path)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    _flush_extra_columns_to_hdf5(df, path, verbose=False)
    df.close()
    assert datasets(path) == {'x': [1.0, 2.0, 3.0], 'y': [10, 20, 30]}


def test_the_metadata_is_updated_with_values_or_with_functions_of_the_metadata(tmp_path):
    path = catalogue(tmp_path, count=1)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    _flush_extra_columns_to_hdf5(df, path, update_metadata={'name': 'text', 'ratio': 2.5, 'flag': True, 'count': 7,
                                                            'next': lambda attrs: attrs['count'] + 1}, verbose=False)
    df.close()
    attrs = attributes(path)
    assert (attrs['name'], attrs['ratio'], bool(attrs['flag']), attrs['count'], attrs['next']) == ('text', 2.5, True, 7, 8)


def test_without_metadata_to_update_the_metadata_is_left_alone(tmp_path):
    path = catalogue(tmp_path, count=1)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    _flush_extra_columns_to_hdf5(df, path, update_metadata=None, verbose=False)
    df.close()
    assert attributes(path) == {'count': 1}


def test_what_is_flushed_is_said_when_it_is_asked_to(tmp_path, capsys):
    path = catalogue(tmp_path)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    df['z'] = df.x * 2
    df['x'] = df.x + 1
    _flush_extra_columns_to_hdf5(df, path, with_columns=['x'], verbose=True)
    df.close()
    out = capsys.readouterr().out
    assert f"Exported the following quantities to {path}\n['z']" in out
    assert "Overwritten the following quantities to" in out and "['x']" in out


def test_nothing_is_said_when_it_is_not_asked_to(tmp_path, capsys):
    path = catalogue(tmp_path)
    df = vaex.open(path, group=STARCATALOG_GROUP)
    df['z'] = df.x * 2
    df['x'] = df.x + 1
    _flush_extra_columns_to_hdf5(df, path, with_columns=['x'], verbose=False)
    df.close()
    assert capsys.readouterr().out == ''


# -- applying a pipeline ------------------------------------------------------------

def test_a_pipeline_is_applied_to_a_table_with_the_other_arguments(tmp_path):
    path = catalogue(tmp_path)
    df = vaex.open(path, group=STARCATALOG_GROUP)

    def pipeline(table, factor):
        table['z'] = table.x * factor
    _decorate_post_processing(pipeline)(df, 3)
    assert df['z'].to_numpy().tolist() == [3.0, 6.0, 9.0]
    assert datasets(path) == {'x': [1.0, 2.0, 3.0], 'y': [10, 20, 30]}                  # a table is not flushed
    df.close()


def test_a_pipeline_is_applied_to_each_table_of_a_list(tmp_path):
    frames = [vaex.open(catalogue(tmp_path, name), group=STARCATALOG_GROUP) for name in ('a.h5', 'b.h5')]
    seen = []

    def pipeline(table):
        seen.append(len(table))
        table['z'] = table.x + 1
    assert _decorate_post_processing(pipeline)(frames) is None
    assert seen == [3, 3] and all(df['z'].to_numpy().tolist() == [2.0, 3.0, 4.0] for df in frames)
    for df in frames:
        df.close()


def test_a_pipeline_is_applied_to_the_file_that_it_is_given_and_what_it_makes_is_flushed(tmp_path):
    path = catalogue(tmp_path)

    def pipeline(table, factor):
        table['z'] = table.x * factor
    _decorate_post_processing(pipeline, hdf5_path_input=True, verbose=False)(path, 5)
    assert datasets(path) == {'x': [1.0, 2.0, 3.0], 'y': [10, 20, 30], 'z': [5.0, 10.0, 15.0]}


def test_the_files_of_a_list_are_each_opened_and_flushed(tmp_path):
    paths = [catalogue(tmp_path, name) for name in ('a.h5', 'b.h5')]

    def pipeline(table):
        table['z'] = table.x * 2
    _decorate_post_processing(pipeline, hdf5_path_input=True, verbose=False)(paths)
    assert all(datasets(path)['z'] == [2.0, 4.0, 6.0] for path in paths)


def test_the_columns_that_it_is_asked_to_overwrite_and_the_metadata_are_flushed_too(tmp_path):
    path = catalogue(tmp_path, count=1)

    def pipeline(table):
        table['x'] = table.x * 10
    _decorate_post_processing(pipeline, hdf5_path_input=True, flush_with_columns=['x'], update_metadata={'count': 2}, verbose=False)(path)
    assert datasets(path)['x'] == [10.0, 20.0, 30.0] and attributes(path)['count'] == 2


def test_the_executor_of_vaex_is_what_it_was_after_a_file_is_processed(tmp_path):
    path = catalogue(tmp_path)
    before = vaex.dataframe.main_executor
    _decorate_post_processing(lambda table: None, hdf5_path_input=True, max_thread_workers=2, verbose=False)(path)
    assert vaex.dataframe.main_executor is before


def test_without_a_file_the_executor_of_vaex_is_not_touched(tmp_path):
    df = vaex.open(catalogue(tmp_path), group=STARCATALOG_GROUP)
    before = vaex.dataframe.main_executor
    _decorate_post_processing(lambda table: None, max_thread_workers=2)(df)
    assert vaex.dataframe.main_executor is before
    df.close()


if __name__ == '__main__':
    pass
