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
Tests of the two functions that the processes of Output run, called here in this process (a process pool is not followed by the
measure of the coverage, and its failures are harder to read): the one that writes the HDF5 file of a partition from the files of the
backend, which is checked against the file that the processes made, and the one that redefines the partitions of a file of the
backend with a rule of the user. Also what happens when no star is in any partition.
"""
import json

import ebf
import h5py
import numpy as np
import pytest

from ..galaxia_ananke._constants import NAME, PARENTPARTICLE_GROUP, PARTITIONING_GROUP, STARCATALOG_GROUP
from ..galaxia_ananke.__metadata__ import __url__, __version__
from ..galaxia_ananke.Output import Output
from .synthetic import log_text, prepared_survey, read_output, write_backend_output


def read_h5(path):
    with h5py.File(path, 'r') as f5:
        return {'groups': {group: {key: f5[group][key][...] for key in f5[group]} for group in f5}, 'attrs': dict(f5.attrs)}


def write_directly(output, partition, path, verbose=False):
    """The file of a partition, written by the function that the processes run, with what the output gives it."""
    partitions = list(output._hdf5s)
    Output._singlethread_ebf_to_hdf5(partitions.index(partition), partition, partitions, path,
                                     output._Output__ebfs_part_slices[partition], output._Output__ebfs_part_lengths[partition],
                                     output._ebfs, output.export_keys, output.all_metadata, verbose)
    return path


@pytest.fixture
def made(tmp_path, monkeypatch):
    return read_output(tmp_path, monkeypatch, n_gens=(0, 1), max_pp_workers=2, pp_auto_flush=True, n_stars=5, partitions=[0, 0, 1, 1, 1])


# -- the file of a partition --------------------------------------------------------

@pytest.mark.parametrize("partition", [0, 1])
def test_the_function_of_the_processes_writes_the_file_that_the_processes_wrote(partition, made, tmp_path):
    output, _ = made
    directly = read_h5(write_directly(output, partition, tmp_path / 'directly.h5'))
    by_processes = read_h5(output._hdf5s[partition])
    assert set(directly['groups']) == set(by_processes['groups'])
    for group in by_processes['groups']:
        assert {k: v.tolist() for k, v in directly['groups'][group].items()} == {k: v.tolist() for k, v in by_processes['groups'][group].items()}
    assert {k: v.tolist() if hasattr(v, 'tolist') else v for k, v in directly['attrs'].items()} == {
        k: v.tolist() if hasattr(v, 'tolist') else v for k, v in by_processes['attrs'].items()}


def test_the_function_of_the_processes_says_what_it_exports_when_it_is_asked_to(made, tmp_path, capsys):
    output, _ = made
    write_directly(output, 1, tmp_path / 'directly.h5', verbose=True)
    out = capsys.readouterr().out
    assert out.count("Exported the following quantities from") == 2 and out.count("for partition 1") == 2          # a file of the backend each
    assert str(output._ebfs[0].resolve()) in out or str(output._ebfs[0]) in out
    assert "'px'" in out and "'partitionid'" in out


def test_the_function_of_the_processes_says_nothing_when_it_is_not_asked_to(made, tmp_path, capsys):
    output, _ = made
    capsys.readouterr()
    write_directly(output, 0, tmp_path / 'directly.h5', verbose=False)
    assert capsys.readouterr().out == ''


def test_the_function_of_the_processes_writes_a_file_without_stars_for_a_partition_without_stars(made, tmp_path):
    output, written = made
    path = tmp_path / 'without_stars.h5'
    Output._singlethread_ebf_to_hdf5(0, -1, [-1], path, {}, {}, output._ebfs, output.export_keys, output.all_metadata, False)
    h5 = read_h5(path)
    stars = h5['groups'][STARCATALOG_GROUP]
    assert set(stars) == set(output.export_keys) and all(column.shape == (0,) for column in stars.values())
    assert all(stars[key].dtype == column.dtype for key, column in written[0][1].items())
    assert h5['attrs']['Cat_partition_id'] == -1 and h5['groups'][PARTITIONING_GROUP]['partition_id'].tolist() == [-1]


def test_the_file_of_the_stars_of_no_partition_has_datasets_without_stars(tmp_path, monkeypatch):
    output, written = read_output(tmp_path, monkeypatch, n_stars=3, partitions=[-1, -1, -1], pp_auto_flush=True)
    assert list(output._hdf5s) == [-1]
    h5 = read_h5(output._hdf5s[-1])
    stars = h5['groups'][STARCATALOG_GROUP]
    assert set(stars) == set(output.export_keys) and all(column.shape == (0,) for column in stars.values())
    assert all(stars[key].dtype == column.dtype for key, column in written[0][1].items())
    assert h5['attrs']['Cat_partition_id'] == -1 and h5['groups'][PARTITIONING_GROUP]['partition_id'].tolist() == [-1]
    assert h5['groups'][PARENTPARTICLE_GROUP] == {} and isinstance(json.loads(h5['attrs']['Galaxia_ebf_output_logs']), (dict, list))


def test_a_catalogue_without_a_star_in_any_partition_is_post_processed_into_an_empty_catalogue_with_all_its_columns(tmp_path, monkeypatch):
    output, _ = read_output(tmp_path, monkeypatch, n_stars=3, partitions=[-1, -1, -1], pp_auto_flush=True)
    with pytest.warns(UserWarning, match="no star in any partition"):
        output.post_process_output()
    h5 = read_h5(output._hdf5s[-1])
    stars, parents = h5['groups'][STARCATALOG_GROUP], h5['groups'][PARENTPARTICLE_GROUP]
    assert set(stars) == set(output.catalogue_keys) and all(column.shape == (0,) for column in stars.values())
    assert set(parents) == {'parentid', 'pos3', 'vel3', 'massinit', 'age', 'feh', 'h_cubic', 'star_count', 'minit_sum',
                            'parent_partition_count', 'parent_star_count', 'parent_minit_sum'}
    assert all(len(dataset) == 0 for dataset in parents.values()) and parents['pos3'].shape == (0, 3) and parents['vel3'].shape == (0, 3)
    assert {key: str(parents[key].dtype) for key in ('parentid', 'star_count', 'parent_partition_count', 'parent_star_count')} == {
        'parentid': 'int64', 'star_count': 'int64', 'parent_partition_count': 'int64', 'parent_star_count': 'int64'}
    assert str(parents['minit_sum'].dtype) == 'float64' and str(parents['parent_minit_sum'].dtype) == 'float64'
    assert len(output._vaex) == 0 and set(output._vaex.get_column_names()) == set(output.catalogue_keys)
    assert len(output) == 0 and not output                                  # as an empty table is: an output without a star is false


# -- redefining the partitions of a file of the backend -----------------------------

def odd_or_even_px(df):
    return (df['px'].to_numpy() % 2).astype(int)


def test_the_function_of_the_processes_that_redefines_the_partitions_reorders_the_file_and_logs_it(tmp_path, monkeypatch):
    survey = prepared_survey(tmp_path, monkeypatch)
    path, columns = write_backend_output(survey, n_stars=6)[0]
    Output._singlethread_redefine_partitions(path, odd_or_even_px, survey.output.export_keys)
    assert ebf.read(str(path), '/partitionid').tolist() == [0, 0, 0, 1, 1, 1]
    pairs = set(zip(ebf.read(str(path), '/px').tolist(), ebf.read(str(path), '/parentid').tolist(), ebf.read(str(path), '/minit').tolist()))
    assert pairs == set(zip(columns['px'].tolist(), columns['parentid'].tolist(), columns['minit'].tolist()))      # the rows stay together
    assert sorted(ebf.read(str(path), '/px').tolist()[:3]) == [2.0, 4.0, 6.0]
    log = str(ebf.read(str(path), '/log'))
    assert log.startswith(log_text(0).rstrip('\n')) and f"# Modified by Python's {NAME} v{__version__}," in log and f"<{__url__}>" in log


def test_the_rule_is_given_an_empty_table_first_to_find_out_which_keys_it_uses_and_then_only_those(tmp_path, monkeypatch):
    survey = prepared_survey(tmp_path, monkeypatch)
    path, _ = write_backend_output(survey, n_stars=6)[0]
    seen = []

    def rule(df):
        seen.append((len(df), sorted(df.columns), sorted(getattr(df, 'record_of_all_used_keys', []))))
        return (df['px'].to_numpy() % 2).astype(int)
    Output._singlethread_redefine_partitions(path, rule, survey.output.export_keys)
    assert [(length, columns) for length, columns, _ in seen] == [(0, sorted(survey.output.export_keys)), (6, ['px'])]
    assert seen[0][2] == [] and 'px' in seen[0][1]


if __name__ == '__main__':
    pass
