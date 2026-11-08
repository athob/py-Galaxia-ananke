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
Tests of how Output tells which stars of the files of the backend belong to which partition: the indices of the stars of each
partition in each file, how many there are, the runs of consecutive indices that they make (which are what is copied at once), and
the names of the HDF5 files that there will be one of for each partition. The partitions are those that the synthetic files give.
"""
import pytest

from .synthetic import prepared_survey, write_backend_output


def indices(output):
    """{partition: {file name: [indices of its stars in the file]}}"""
    return {int(partition): {name: found.tolist() for name, found in per_file.items()}
            for partition, per_file in output._Output__ebfs_partitions.items()}


def make(tmp_path, monkeypatch, n_gens=(0,), **backend):
    survey = prepared_survey(tmp_path, monkeypatch, n_gens=n_gens)
    written = write_backend_output(survey, **backend)
    return survey.output, [path.name for path, _ in written]


# -- the stars of each partition ----------------------------------------------------

def test_the_stars_of_each_partition_are_those_that_have_its_number(tmp_path, monkeypatch):
    output, (name,) = make(tmp_path, monkeypatch, n_stars=5, partitions=[0, 0, 1, 1, 2])
    assert indices(output) == {0: {name: [0, 1]}, 1: {name: [2, 3]}, 2: {name: [4]}}


def test_the_stars_of_a_partition_do_not_have_to_be_next_to_each_other(tmp_path, monkeypatch):
    output, (name,) = make(tmp_path, monkeypatch, n_stars=5, partitions=[1, 0, 1, 1, 0])
    assert indices(output) == {0: {name: [1, 4]}, 1: {name: [0, 2, 3]}}


def test_the_stars_that_are_in_no_partition_are_left_out(tmp_path, monkeypatch):
    output, (name,) = make(tmp_path, monkeypatch, n_stars=5, partitions=[0, -1, 0, -2, 1])
    assert indices(output) == {0: {name: [0, 2]}, 1: {name: [4]}}


def test_a_partition_has_the_stars_of_all_the_files(tmp_path, monkeypatch):
    output, (first, second) = make(tmp_path, monkeypatch, n_gens=(0, 1), n_stars=3, partitions=[0, 1, 0])
    assert indices(output) == {0: {first: [0, 2], second: [0, 2]}, 1: {first: [1], second: [1]}}


def test_without_a_single_star_in_a_partition_there_is_the_partition_of_nothing(tmp_path, monkeypatch):
    output, _ = make(tmp_path, monkeypatch, n_stars=3, partitions=[-1, -1, -3])
    assert indices(output) == {-1: {}}


def test_the_number_of_stars_of_each_partition_in_each_file_and_in_all(tmp_path, monkeypatch):
    output, (first, second) = make(tmp_path, monkeypatch, n_gens=(0, 1), n_stars=5, partitions=[0, 0, 1, 1, 1])
    assert {int(i): lengths for i, lengths in output._Output__ebfs_part_lengths.items()} == {0: {first: 2, second: 2}, 1: {first: 3, second: 3}}
    assert {int(i): length for i, length in output._Output__partitions_lengths.items()} == {0: 4, 1: 6}


# -- the runs of consecutive stars --------------------------------------------------

def slices(output):
    return {int(partition): per_file for partition, per_file in output._Output__ebfs_part_slices.items()}


@pytest.mark.parametrize("partitions, expected", [
    ([0, 0, 1, 1], {0: [slice(0, 2)], 1: [slice(2, 4)]}),                                   # one run each
    ([0, 0, 0, 0], {0: [slice(0, 4)]}),
    ([0, 1, 0, 0, 1], {0: [slice(0, 1), slice(2, 4)], 1: [slice(1, 2), slice(4, 5)]}),       # two runs each
    ([1, 0, 1, 1, 0, 0, 1], {0: [slice(1, 2), slice(4, 6)], 1: [slice(0, 1), slice(2, 4), slice(6, 7)]}),
    ([0, 1, 2, 1, 0], {0: [slice(0, 1), slice(4, 5)], 1: [slice(1, 2), slice(3, 4)], 2: [slice(2, 3)]}),
])
def test_the_stars_of_a_partition_are_copied_in_runs_of_consecutive_stars(partitions, expected, tmp_path, monkeypatch):
    output, (name,) = make(tmp_path, monkeypatch, n_stars=len(partitions), partitions=partitions)
    assert slices(output) == {partition: {name: runs} for partition, runs in expected.items()}


# -- the files of the partitions ----------------------------------------------------

def hdf5_names(output):
    return {int(partition): path.name for partition, path in output._hdf5s.items()}


def test_there_is_a_file_for_each_partition_named_after_the_output_and_its_number(tmp_path, monkeypatch):
    output, _ = make(tmp_path, monkeypatch, n_stars=5, partitions=[0, 0, 1, 1, 2])
    assert hdf5_names(output) == {i: f"{output.output_name}.{i}.h5" for i in (0, 1, 2)}
    assert all(path.parent == output.output_dir for path in output._hdf5s.values())


def test_the_numbers_of_the_partitions_in_the_names_of_the_files_have_as_many_digits_as_the_largest(tmp_path, monkeypatch):
    output, _ = make(tmp_path, monkeypatch, n_stars=12, partitions=list(range(12)))
    assert hdf5_names(output) == {i: f"{output.output_name}.{i:02d}.h5" for i in range(12)}


def test_the_partitions_are_in_the_order_of_their_numbers_whatever_the_order_of_the_stars(tmp_path, monkeypatch):
    output, _ = make(tmp_path, monkeypatch, n_stars=4, partitions=[3, 1, 3, 2])
    assert [int(partition) for partition in output._hdf5s] == [1, 2, 3]


def test_the_partition_of_nothing_has_the_file_of_the_partition_minus_one(tmp_path, monkeypatch):
    output, _ = make(tmp_path, monkeypatch, n_stars=3, partitions=[-1, -1, -1])
    assert hdf5_names(output) == {-1: f"{output.output_name}.-1.h5"}


if __name__ == '__main__':
    pass
