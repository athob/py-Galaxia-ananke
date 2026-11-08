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
Tests of the hash of an Input, the names of the files that follow from it, the files that an Input writes for the backend and
how it prepares them (the directory where the backend finds its inputs is a temporary one here, never the cache of the
installation). The hashes are pinned: they name the cache of an input, so a change of them leaves every earlier cache unused.
"""
import os
import re
import types

import ebf
import numpy as np
import pytest

from ..galaxia_ananke import Input
from .utils import minimal_kernels, minimal_particles, patch_nbody1

# The hashes of the smallest inputs, worked out outside Input: the sha256 of the hashes of the arrays in the order that the
# particles have once the keys that are missing are added (the keys in lexical order, then the parent and the partition if they
# were not given, then the kernels), each array in the order of the partitions, as raw numbers.
HASH_OF_THE_SMALLEST_INPUT = 'eed86126b03252643420dd0a683aae8c7243da76246f73784a03abc144808628'
HASH_WITH_KERNELS_TWICE_AS_LARGE = 'db92c672d4738e0b019f296c5c2d2793b0fe3060fb81cd77b7e39888561eac41'
HASH_WITH_THE_PARTICLES_IN_THE_OTHER_PARTITION_ORDER = '2aa9bfbddafdee881e1bd79ae4c4f5c07f789f359a104bba2ce5b2a56f92f23a'
HASH_WITH_THE_PARTITION_GIVEN_AS_ZEROS = 'ed282020dcc2ff54bf00d6357af111453237125c8f76b66313ba077406842af0'


def make(tmp_path, particles=None, kernels=None, sorted_by_partition=True, **kwargs):
    kwargs = {'input_dir': tmp_path, 'name': 'sim', 'ngb': 64, **kwargs}
    inp = Input(minimal_particles() if particles is None else particles, minimal_kernels() if kernels is None else kernels, **kwargs)
    if sorted_by_partition:
        inp.input_sorter = None                             # as a survey does before it uses the input
    return inp


@pytest.fixture
def nbody1(tmp_path, monkeypatch):
    """Where the backend finds its inputs, under a temporary directory and not in the cache of the installation."""
    return patch_nbody1(monkeypatch, tmp_path)


# -- the sorter ---------------------------------------------------------------------

def test_the_default_sorter_puts_the_particles_in_the_order_of_their_partitions(tmp_path):
    particles = {**minimal_particles(), 'partitionid': np.array([1, 0])}
    assert make(tmp_path, particles).input_sorter.tolist() == [1, 0]
    assert make(tmp_path).input_sorter.tolist() == [0, 1]                  # one partition: as they are


def test_a_sorter_that_is_given_is_the_one_used(tmp_path):
    inp = make(tmp_path)
    inp.input_sorter = np.array([1, 0])
    assert inp.input_sorter.tolist() == [1, 0]


def test_a_new_sorter_makes_a_new_hash_and_the_old_one_comes_back_with_the_old_sorter(tmp_path):
    inp = make(tmp_path)
    first = inp.hash
    inp.input_sorter = np.array([1, 0])
    assert inp.hash != first
    inp.input_sorter = None
    assert inp.hash == first


@pytest.mark.xfail(raises=AttributeError, reason="the sorter only exists once it is set, which a survey does: an Input that was "
                   "just made cannot give its hash, nor the names of its files, which are built from it")
@pytest.mark.parametrize("attribute", ["hash", "name_hash", "pname", "kname"])
def test_an_input_that_was_just_made_can_say_its_hash_and_the_names_of_its_files(attribute, tmp_path):
    assert getattr(make(tmp_path, sorted_by_partition=False), attribute)


# -- the hash -----------------------------------------------------------------------

@pytest.mark.parametrize("what, expected", [
    ({}, HASH_OF_THE_SMALLEST_INPUT),
    ({'k_factor': 2.0}, HASH_WITH_KERNELS_TWICE_AS_LARGE),
    ({'particles': {**minimal_particles(), 'partitionid': np.array([1, 0])}}, HASH_WITH_THE_PARTICLES_IN_THE_OTHER_PARTITION_ORDER),
    ({'particles': {**minimal_particles(), 'partitionid': np.array([0, 0])}}, HASH_WITH_THE_PARTITION_GIVEN_AS_ZEROS),
], ids=["the smallest input", "kernels twice as large", "particles of two partitions", "the partition given as zeros"])
def test_the_hash_of_an_input_is_what_it_always_was(what, expected, tmp_path):
    inp = make(tmp_path, **what)
    assert inp.hash == expected
    assert inp.name_hash == f"sim_{expected[:7]}"


def test_spelling_out_a_default_changes_the_hash():
    # How it is today: the arrays are hashed in the order the keys have once the defaults are added, so that a partition given as
    # zeros, which is what it is by default, comes in another place and makes another hash than the one of the same input without it
    assert HASH_OF_THE_SMALLEST_INPUT != HASH_WITH_THE_PARTITION_GIVEN_AS_ZEROS


def test_the_hash_is_a_sha256_in_hexadecimal(tmp_path):
    assert re.fullmatch(r"[0-9a-f]{64}", make(tmp_path).hash)


def test_the_hash_only_depends_on_the_arrays(tmp_path):
    (tmp_path / 'elsewhere').mkdir()
    other = make(tmp_path / 'elsewhere', name='another', ngb=8, append_hash=False)
    assert other.hash == make(tmp_path).hash


@pytest.mark.parametrize("change", ["age", "kernels", "k_factor", "optional key", "precision"])
def test_the_hash_changes_with_what_is_in_the_arrays(change, tmp_path):
    particles, kernels, kwargs = minimal_particles(), minimal_kernels(), {}
    if change == "age":
        particles['age'][0] = 9.75
    elif change == "kernels":
        kernels[1, 1] = 4.5
    elif change == "k_factor":
        kwargs['k_factor'] = 1.5
    elif change == "optional key":
        particles['dform'] = np.array([0.0, 1.0])
    elif change == "precision":
        particles['age'] = particles['age'].astype('float32')
    assert make(tmp_path, particles, kernels, **kwargs).hash != make(tmp_path).hash


def test_the_hash_is_made_once_so_changing_an_array_later_does_not_change_it(tmp_path):
    # How it is today: a cache would then name files that hold something else than what the name says
    particles = minimal_particles()
    inp = make(tmp_path, particles)
    first = inp.hash
    particles['age'][0] = -999.0
    assert inp.hash == first


def test_the_metadata_of_an_input_is_its_name_and_its_hash(tmp_path):
    inp = make(tmp_path, name='mine')
    assert inp.metadata == {'name': 'mine', 'hash': inp.hash}


# -- the names of the files ---------------------------------------------------------

def test_the_files_are_named_after_the_name_the_hash_the_dimension_and_the_neighbours(tmp_path):
    inp = make(tmp_path)
    assert inp.pname == tmp_path / f"{inp.name_hash}.ebf"
    assert inp.kname == tmp_path / f"{inp.name_hash}_d6n64_den.ebf"
    assert inp._hashname == tmp_path / f"{inp.name_hash}.hash"


def test_the_files_of_kernels_in_position_space_say_so(tmp_path):
    inp = make(tmp_path, kernels=np.array([1.0, 2.0]), ngb=32)
    assert inp.kname.name.endswith("_d3n32_den.ebf")


def test_the_hash_can_be_left_out_of_the_names(tmp_path):
    inp = make(tmp_path, append_hash=False)
    assert inp.name_hash == 'sim' and inp.pname == tmp_path / 'sim.ebf'


@pytest.mark.xfail(raises=AssertionError, reason="the names of the files are made by replacing the suffix of the name with the "
                   "extension, and a name with a dot (m12i.res7100, for one) has its hash as its suffix: the files of two "
                   "inputs of that name have the same names, whatever their hash")
@pytest.mark.parametrize("name", ["run.v2", "m12i.res7100"])
def test_a_name_with_a_dot_keeps_its_hash_in_the_names_of_its_files(name, tmp_path):
    inp = make(tmp_path, name=name)
    assert inp.pname.name == f"{inp.name_hash}.ebf"


# -- the files ----------------------------------------------------------------------

def read(path, key):
    return ebf.read(str(path), f"/{key}")


def test_the_files_are_written_with_their_hash(tmp_path, nbody1):
    inp = make(tmp_path)
    inp._write_ebf_files()
    assert inp.pname.is_file() and inp.kname.is_file()
    assert inp._hashname.read_bytes() == inp.hash.encode('ascii')


def test_the_particles_file_holds_each_array_in_the_order_of_the_partitions(tmp_path, nbody1):
    inp = make(tmp_path, {**minimal_particles(), 'partitionid': np.array([1, 0])})
    inp._write_ebf_files()
    np.testing.assert_array_equal(read(inp.pname, 'massinit'), [2000.0, 1000.0])
    np.testing.assert_array_equal(read(inp.pname, 'pos3'), [[4.0, 5.0, 6.0], [1.0, 2.0, 3.0]])
    np.testing.assert_array_equal(read(inp.pname, 'partitionid'), [0, 1])
    np.testing.assert_array_equal(read(inp.pname, 'parentid'), [1, 0])      # the particles that came second are first


def test_the_particles_file_holds_every_key_and_zeros_for_those_that_were_not_given(tmp_path, nbody1):
    inp = make(tmp_path, {**minimal_particles(), 'dform': np.array([0.5, 1.5])})
    inp._write_ebf_files()
    for key in Input._required_keys_in_particles | Input._optional_keys_in_particles:
        assert read(inp.pname, key) is not None
    np.testing.assert_array_equal(read(inp.pname, 'dform'), [0.5, 1.5])
    for key in ('helium', 'calcium', 'alpha', 'id'):
        np.testing.assert_array_equal(read(inp.pname, key), [0.0, 0.0])


def test_the_kernels_file_holds_the_kernels_scaled_and_in_the_order_of_the_partitions(tmp_path, nbody1):
    inp = make(tmp_path, {**minimal_particles(), 'partitionid': np.array([1, 0])}, k_factor=2.0)
    inp._write_ebf_files()
    np.testing.assert_array_equal(read(inp.kname, 'h_cubic'), [[6.0, 8.0], [2.0, 4.0]])


# -- writing again, or not ----------------------------------------------------------

OLD = 1_000_000_000


def age_files(inp):
    paths = [inp.pname, inp.kname, inp._hashname]
    for path in paths:
        os.utime(path, (OLD, OLD))
    return paths


def modification_times(paths):
    return [path.stat().st_mtime for path in paths]


def test_without_caching_the_files_are_written_again(tmp_path, nbody1):
    inp = make(tmp_path)
    inp._write_ebf_files()
    paths = age_files(inp)
    inp._write_ebf_files()
    assert all(time > OLD + 1000 for time in modification_times(paths))


def test_with_caching_files_that_are_there_for_the_same_input_are_left_alone(tmp_path, nbody1):
    with pytest.warns(DeprecationWarning, match="caching mode"):
        inp = make(tmp_path, caching=True)
    inp._write_ebf_files()
    paths = age_files(inp)
    inp._write_ebf_files()
    assert modification_times(paths) == [OLD, OLD, OLD]


@pytest.mark.parametrize("missing", ["pname", "kname", "_hashname"])
def test_with_caching_a_file_that_is_missing_makes_all_of_them_written_again(missing, tmp_path, nbody1):
    with pytest.warns(DeprecationWarning, match="caching mode"):
        inp = make(tmp_path, caching=True)
    inp._write_ebf_files()
    paths = age_files(inp)
    getattr(inp, missing).unlink()
    inp._write_ebf_files()
    assert all(path.is_file() for path in paths)
    assert all(time > OLD + 1000 for time in modification_times(paths))


def test_with_caching_files_of_another_input_with_the_same_names_are_written_over(tmp_path, nbody1):
    with pytest.warns(DeprecationWarning, match="caching mode"):
        first = make(tmp_path, caching=True, append_hash=False)
    first._write_ebf_files()
    changed = minimal_particles()
    changed['age'][0] = 9.75
    with pytest.warns(DeprecationWarning, match="caching mode"):
        second = make(tmp_path, changed, caching=True, append_hash=False)
    assert second.pname == first.pname and second.hash != first.hash       # the same files, for other numbers
    second._write_ebf_files()
    np.testing.assert_array_equal(read(second.pname, 'age'), [9.75, 10.0])
    assert second._hashname.read_bytes() == second.hash.encode('ascii')


# -- what the backend is given ------------------------------------------------------

def test_the_backend_finds_the_files_through_a_directory_and_a_list_of_its_inputs(tmp_path, nbody1):
    inp = make(tmp_path)
    listing = inp._write_ebf_files()
    assert listing == nbody1 / 'filenames' / f"{inp.name_hash}.txt"
    assert listing.read_text() == f"nbody1/{inp.name_hash}/\n\t1\t1\n{inp.pname.name}\n"
    for path in (inp.pname, inp.kname):
        link = nbody1 / inp.name_hash / path.name
        assert link.is_symlink() and link.resolve() == path.resolve()


def test_preparing_the_input_writes_the_parameters_and_the_files_and_says_where_they_are(tmp_path, nbody1):
    inp = make(tmp_path)
    survey = types.SimpleNamespace(_write_parameter_file=lambda: (tmp_path / 'survey_params', {'fsample': 0.01}))
    assert inp.prepare_input(survey) == (inp.name_hash, tmp_path / 'survey_params', {'fsample': 0.01})
    assert inp.pname.is_file() and (nbody1 / 'filenames' / f"{inp.name_hash}.txt").is_file()


# -- an Input made from files that exist --------------------------------------------

def test_an_input_made_from_the_files_of_another_has_what_they_hold(tmp_path, nbody1):
    first = make(tmp_path, {**minimal_particles(), 'dform': np.array([0.5, 1.5])}, k_factor=2.0, ngb=32)
    first._write_ebf_files()
    again = Input(pname=first.pname, kname=first.kname)
    assert first.kname.name.endswith('_d6n32_den.ebf')
    assert again.name == first.name_hash and again.ngb == 32 and again.k_factor == 1.0      # ngb is read from the name of the kernels file
    assert again._input_dir == tmp_path and again.pname == first.pname and again.kname == first.kname
    assert again.hdim == 6 and again.length == 2
    np.testing.assert_array_equal(again.particles['massinit'], [1000.0, 2000.0])
    np.testing.assert_array_equal(again.particles['dform'], [0.5, 1.5])
    np.testing.assert_array_equal(again.kernels, [[2.0, 4.0], [6.0, 8.0]])         # as the file holds them: already scaled
    assert set(again.keys()) == Input.all_possible_keys_in_particles


def test_an_input_made_from_files_does_not_write_them_again(tmp_path, nbody1):
    first = make(tmp_path)
    first._write_ebf_files()
    again = Input(pname=first.pname, kname=first.kname)
    again.input_sorter = None
    paths = [first.pname, first.kname]
    for path in paths:
        os.utime(path, (OLD, OLD))
    again._write_ebf_files()
    assert modification_times(paths) == [OLD, OLD]
    assert again._hashname.is_file()


if __name__ == '__main__':
    pass
