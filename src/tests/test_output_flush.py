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
Tests of the pipelines that a user applies to the catalogue with Output.apply_post_process_pipeline_and_flush (which runs them on
the HDF5 file of each partition, in processes, and flushes what they made), of the options of it, and of the flushing of the columns
of the loaded catalogue by itself. The catalogue is made of the HDF5 files of synthetic output of the backend.
"""
import h5py
import numpy as np
import pytest

from ..galaxia_ananke._constants import STARCATALOG_GROUP
from .synthetic import read_output

PARTITIONS = [0, 0, 1, 1, 2, 2]


def add_column(df, factor):
    df['z'] = df['px'] * factor


def double_px(df):
    df['px'] = df['px'] * 2


def stars(path):
    with h5py.File(path, 'r') as f5:
        return {key: f5[STARCATALOG_GROUP][key][...] for key in f5[STARCATALOG_GROUP]}


def attrs(path):
    with h5py.File(path, 'r') as f5:
        return dict(f5.attrs)


@pytest.fixture
def made(tmp_path, monkeypatch):
    output, written = read_output(tmp_path, monkeypatch, max_pp_workers=2, pp_auto_flush=True, n_stars=6, partitions=PARTITIONS)
    return output, written[0][1]


# -- a pipeline run on the file of each partition -----------------------------------

def test_a_pipeline_is_run_on_each_partition_with_its_arguments_and_what_it_makes_is_in_the_file(made):
    output, columns = made
    output.apply_post_process_pipeline_and_flush(add_column, 3)
    for partition, path in output._hdf5s.items():
        rows = np.asarray(columns['partitionid']) == partition
        assert stars(path)['z'].tolist() == (3 * columns['px'][rows]).tolist()


def test_the_catalogue_is_loaded_with_what_the_pipeline_made_once_it_is_flushed(made):
    output, columns = made
    output.apply_post_process_pipeline_and_flush(add_column, 2)
    assert len(output._vaex) == 6 and sorted(output._vaex.get_column_names()).count('z') == 1
    assert output._vaex['z'].to_numpy().tolist() == (2 * columns['px']).tolist()


def test_the_columns_that_it_is_asked_to_flush_are_overwritten_in_the_files(made):
    output, columns = made
    output.apply_post_process_pipeline_and_flush(double_px, flush_with_columns=['px'])
    assert np.concatenate([stars(path)['px'] for path in output._hdf5s.values()]).tolist() == (2 * columns['px']).tolist()


def test_the_columns_that_it_is_not_asked_to_flush_are_left_as_they_were(made):
    output, columns = made
    output.apply_post_process_pipeline_and_flush(double_px)
    assert np.concatenate([stars(path)['px'] for path in output._hdf5s.values()]).tolist() == columns['px'].tolist()


def test_the_metadata_that_it_is_asked_to_update_is_updated_in_every_file(made):
    output, _ = made
    output.apply_post_process_pipeline_and_flush(add_column, 1, update_metadata={'checked': 7, 'next': lambda a: a['Cat_partition_id'] + 100})
    for partition, path in output._hdf5s.items():
        assert attrs(path)['checked'] == 7 and attrs(path)['next'] == partition + 100


def test_the_catalogue_is_not_loaded_when_it_is_asked_not_to_reload(made):
    output, columns = made
    output.apply_post_process_pipeline_and_flush(add_column, 3, hold_reload=True)
    assert all('z' in stars(path) for path in output._hdf5s.values())
    with pytest.raises(RuntimeError, match="Don't attempt creating an Output object on your own"):
        output._vaex


def test_the_partitions_can_be_shared_between_the_processes_by_their_number_of_stars(tmp_path, monkeypatch):
    output, written = read_output(tmp_path, monkeypatch, max_pp_workers=2, pp_auto_flush=True, n_stars=6, partitions=[0, 0, 0, 1, 2, 2])
    output.apply_post_process_pipeline_and_flush(add_column, 5, consolidate_partitions_per_process=True)
    columns = written[0][1]
    for partition, path in output._hdf5s.items():
        rows = np.asarray(columns['partitionid']) == partition
        assert stars(path)['z'].tolist() == (5 * columns['px'][rows]).tolist()


@pytest.fixture
def unflushed(tmp_path, monkeypatch):
    """An output whose partitions are not flushed by themselves: the pipelines are run on the tables, and flushed afterwards."""
    output, written = read_output(tmp_path, monkeypatch, max_pp_workers=2, pp_auto_flush=False, n_stars=6, partitions=PARTITIONS)
    return output, written[0][1]


def test_a_pipeline_is_run_on_the_tables_when_the_partitions_are_not_flushed_by_themselves(unflushed):
    output, columns = unflushed
    output.apply_post_process_pipeline_and_flush(add_column, 3)
    for partition, path in output._hdf5s.items():
        rows = np.asarray(columns['partitionid']) == partition
        assert stars(path)['z'].tolist() == (3 * columns['px'][rows]).tolist()
    assert output._vaex['z'].to_numpy().tolist() == (3 * columns['px']).tolist()


def test_the_flushing_can_be_held_until_the_tables_are_flushed_by_hand(unflushed):
    output, columns = unflushed
    output.apply_post_process_pipeline_and_flush(add_column, 3, hold_flush=True, hold_reload=True)
    assert not any('z' in stars(path) for path in output._hdf5s.values())                      # not in the files yet,
    assert all('z' in table.column_names for table in output._vaex_per_partition.values())    # but in the tables
    output.flush_extra_columns_to_hdf5()
    for partition, path in output._hdf5s.items():
        rows = np.asarray(columns['partitionid']) == partition
        assert stars(path)['z'].tolist() == (3 * columns['px'][rows]).tolist()


def test_the_columns_that_are_asked_for_are_overwritten_in_the_files_when_the_partitions_are_not_flushed_by_themselves(unflushed):
    output, columns = unflushed
    output.apply_post_process_pipeline_and_flush(double_px, flush_with_columns=['px'])
    assert np.concatenate([stars(path)['px'] for path in output._hdf5s.values()]).tolist() == (2 * columns['px']).tolist()


# -- flushing what the loaded catalogue has made ------------------------------------

@pytest.fixture
def loaded(tmp_path, monkeypatch):
    output, written = read_output(tmp_path, monkeypatch, max_pp_workers=2, pp_auto_flush=False, n_stars=6, partitions=PARTITIONS)
    return output, written[0][1]


def test_the_new_columns_of_the_tables_of_the_partitions_are_flushed_to_their_files(loaded):
    output, columns = loaded
    for df in output._vaex_per_partition.values():
        df['z'] = df.px * 4
    output.flush_extra_columns_to_hdf5()
    for partition, path in output._hdf5s.items():
        rows = np.asarray(columns['partitionid']) == partition
        assert stars(path)['z'].tolist() == (4 * columns['px'][rows]).tolist()


def test_the_columns_that_are_asked_for_are_overwritten_when_the_tables_are_flushed(loaded):
    output, columns = loaded
    for df in output._vaex_per_partition.values():
        df['px'] = df.px + 1000
    output.flush_extra_columns_to_hdf5(with_columns=['px'])
    assert np.concatenate([stars(path)['px'] for path in output._hdf5s.values()]).tolist() == (columns['px'] + 1000).tolist()


def test_the_metadata_is_updated_when_the_tables_are_flushed(loaded):
    output, _ = loaded
    output.flush_extra_columns_to_hdf5(update_metadata={'flushed': True})
    assert all(bool(attrs(path)['flushed']) for path in output._hdf5s.values())


def test_flushing_without_anything_new_changes_nothing_in_the_files(loaded):
    output, columns = loaded
    before = {partition: stars(path) for partition, path in output._hdf5s.items()}
    output.flush_extra_columns_to_hdf5()
    after = {partition: stars(path) for partition, path in output._hdf5s.items()}
    assert {p: {k: v.tolist() for k, v in cols.items()} for p, cols in before.items()} == {p: {k: v.tolist() for k, v in cols.items()} for p, cols in after.items()}


if __name__ == '__main__':
    pass
