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
import gc
import hashlib
import os
import pathlib
import re
import select
import signal
import types
import warnings

import blake3
import h5py
import numpy as np
import pytest
import vaex
import vaex.file
import vaex.hdf5.dataset
import vaex.utils

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
    df.close()


def test_blank_lines_comments_and_the_spaces_around_the_names_are_left_out_of_the_list(tmp_path):
    first, second = ananke_file(tmp_path, 'a.h5', stars=(1.0, 2.0)), ananke_file(tmp_path, 'b.h5', stars=(3.0,))
    df = _patch_vaex.new_open_many([f"  {first}  \n", "", "   ", "# a comment", "  # another one", f"\t{second}\n"], group=STARCATALOG_GROUP)
    assert values(df, 'x') == [1.0, 2.0, 3.0]
    df.close()


def test_a_list_can_be_a_generator_of_names(tmp_path):
    paths = [ananke_file(tmp_path, f"{name}.h5", stars=(float(i),)) for i, name in enumerate('abc')]
    df = _patch_vaex.new_open_many(map(str, paths), group=STARCATALOG_GROUP)
    assert values(df, 'x') == [0.0, 1.0, 2.0]
    df.close()


def test_the_arguments_of_opening_are_given_to_each_file(tmp_path):
    first, second = ananke_file(tmp_path, 'a.h5'), ananke_file(tmp_path, 'b.h5')
    df = _patch_vaex.new_open_many([str(first), str(second)], group=PARENTPARTICLE_GROUP)
    assert df.get_column_names() == ['parentid'] and len(df) == 10
    df.close()


def test_the_patched_function_is_the_one_that_vaex_has_and_the_original_is_kept():
    assert vaex.open_many is _patch_vaex.new_open_many
    assert _patch_vaex._original_open_many is not _patch_vaex.new_open_many
    assert _patch_vaex._original_hdf5_load is not _patch_vaex.new_hdf5_load
    assert vaex.utils.create_hasher is _patch_vaex.new_create_hasher
    assert _patch_vaex._original_create_hasher is not _patch_vaex.new_create_hasher
    assert vaex.hdf5.dataset.Hdf5MemoryMapped._open is _patch_vaex.new_hdf5_open
    assert _patch_vaex._original_hdf5_open is not _patch_vaex.new_hdf5_open


# -- hashing in a process made by fork -------------------------------------------------
# vaex hashes large data with the thread pool of blake3 (max_threads=AUTO), which a child made by fork inherits without its
# threads: a child that hashes after its parent did waits forever. The patch makes vaex hash with one thread in such a child.

def run_in_a_forked_child(function, timeout):
    """What the function gives in a child made by fork, as a string; None if the child does not answer in time (killed)."""
    read, write = os.pipe()
    pid = os.fork()
    if pid == 0:                                    # the child
        try:
            os.close(read)
            os.write(write, str(function()).encode())
            os._exit(0)
        except BaseException:
            os._exit(1)
    os.close(write)
    try:
        answer = os.read(read, 4096).decode() if select.select([read], [], [], timeout)[0] else None
    finally:
        os.close(read)
        os.kill(pid, signal.SIGKILL)
        os.waitpid(pid, 0)
    return answer


def recording_hasher(seen):
    return lambda initial_data=None, large_data=True: seen.append((initial_data, large_data))


@pytest.mark.parametrize("forked, large_data, asked",
                         [(False, True, True), (False, False, False), (True, True, False), (True, False, False)])
def test_threads_are_asked_only_where_the_process_is_not_a_forked_child(forked, large_data, asked, monkeypatch):
    seen = []
    monkeypatch.setattr(_patch_vaex, '_original_create_hasher', recording_hasher(seen))
    monkeypatch.setattr(_patch_vaex, '_in_forked_child', forked)
    _patch_vaex.new_create_hasher(b'data', large_data=large_data)
    assert seen == [(b'data', asked)]


def test_the_defaults_of_the_hasher_are_those_of_vaex(monkeypatch):
    seen = []
    monkeypatch.setattr(_patch_vaex, '_original_create_hasher', recording_hasher(seen))
    monkeypatch.setattr(_patch_vaex, '_in_forked_child', False)
    _patch_vaex.new_create_hasher()
    assert seen == [(None, True)]


def test_the_digest_does_not_depend_on_the_threads(monkeypatch):
    data = b'\xab' * (4 << 20)
    expected = _patch_vaex._original_create_hasher(data, large_data=True).hexdigest()
    monkeypatch.setattr(_patch_vaex, '_in_forked_child', True)
    assert vaex.utils.create_hasher(data, large_data=True).hexdigest() == expected


@pytest.mark.filterwarnings("ignore:This process:DeprecationWarning")
def test_a_forked_child_knows_that_it_is_one_and_its_parent_does_not():
    assert _patch_vaex._in_forked_child is False
    assert run_in_a_forked_child(lambda: _patch_vaex._in_forked_child, timeout=20) == 'True'


@pytest.mark.filterwarnings("ignore:This process:DeprecationWarning")
def test_a_forked_child_hashes_with_vaex_after_its_parent_used_the_threads_of_blake3():
    data = b'\xab' * (64 << 20)
    expected = blake3.blake3(data, max_threads=blake3.blake3.AUTO).hexdigest()   # the parent now has the pool of threads
    original = lambda: _patch_vaex._original_create_hasher(data, large_data=True).hexdigest()
    without_the_patch = run_in_a_forked_child(original, timeout=5)
    if without_the_patch == expected:
        pytest.skip("a child made by fork does not hang here without the patch: another platform, or blake3 or vaex fixed it")
    assert without_the_patch is None, "the child failed, it did not hang"
    patched = lambda: vaex.utils.create_hasher(data, large_data=True).hexdigest()
    assert run_in_a_forked_child(patched, timeout=15) == expected


def test_no_module_of_vaex_keeps_the_hasher_under_a_name_of_its_own():
    # the patch replaces the attribute of vaex.utils: a module that kept the function under another name would keep
    # the original, so every mention of create_hasher in vaex has to be a call (or its definition)
    mention = re.compile(r"\bcreate_hasher\b(?!\()")
    offenders = [str(path) for directory in map(pathlib.Path, vaex.__path__) for path in directory.rglob('*.py')
                 if not {'test', 'tests'} & set(path.relative_to(directory).parts)
                 and mention.search(path.read_text(encoding='utf-8', errors='replace'))]
    assert offenders == []


def test_applying_the_patch_again_changes_nothing():
    def functions():
        dataset = vaex.hdf5.dataset.Hdf5MemoryMapped
        return vaex.utils.create_hasher, vaex.open_many, dataset._load, dataset._open
    patched = functions()
    with pytest.warns(RuntimeWarning, match="monkey-patches"):
        _patch_vaex.apply_vaex_patch()
    assert functions() == patched


# -- the file object that h5py is given ------------------------------------------------
# vaex opens a file twice, for h5py and for the memory mapping, under the same key of file_map: the file object of h5py is
# finalized unclosed when the HDF5 file is closed, a ResourceWarning at every vaex.open. The patch closes it with that file.

def descriptors_and_mappings_of(path):
    """How many descriptors and memory mappings of this process are of the file (read in /proc)."""
    target = os.path.realpath(path)
    descriptors = 0
    for name in os.listdir('/proc/self/fd'):
        try:
            descriptors += os.path.realpath(os.readlink(f'/proc/self/fd/{name}')) == target
        except OSError:
            pass
    with open('/proc/self/maps') as maps:
        return descriptors, sum(1 for line in maps if line.rstrip().endswith(target))


def resource_warnings_of_opening(path):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        df = vaex.open(str(path), group=STARCATALOG_GROUP)
        gc.collect()
    df.close()
    return [str(w.message) for w in caught if issubclass(w.category, ResourceWarning) and path.name in str(w.message)]


def test_opening_a_file_makes_no_resource_warning(tmp_path, monkeypatch):
    dataset = vaex.hdf5.dataset.Hdf5MemoryMapped
    monkeypatch.setattr(dataset, '_open', _patch_vaex._original_hdf5_open)
    if not resource_warnings_of_opening(ananke_file(tmp_path, 'original.h5')):
        pytest.skip("vaex makes no ResourceWarning when it opens a file: fixed there, the patch is no longer needed")
    monkeypatch.setattr(dataset, '_open', _patch_vaex.new_hdf5_open)
    assert resource_warnings_of_opening(ananke_file(tmp_path, 'patched.h5')) == []


@pytest.mark.skipif(not os.path.isdir('/proc/self/fd'), reason="the descriptors and the mappings are read in /proc")
@pytest.mark.filterwarnings("ignore::ResourceWarning")
def test_the_patch_holds_no_descriptor_or_mapping_that_vaex_does_not(tmp_path, monkeypatch):
    def counts(function, name):
        monkeypatch.setattr(vaex.hdf5.dataset.Hdf5MemoryMapped, '_open', function)
        path = ananke_file(tmp_path, name)
        df = vaex.open(str(path), group=STARCATALOG_GROUP)
        gc.collect()
        while_open = descriptors_and_mappings_of(path)
        df.close()
        del df
        gc.collect()
        return while_open, descriptors_and_mappings_of(path)
    original = counts(_patch_vaex._original_hdf5_open, 'original.h5')
    patched = counts(_patch_vaex.new_hdf5_open, 'patched.h5')
    assert patched == original
    assert patched[1] == (0, 0)


class RecordingFile:
    def __init__(self, journal):
        self.journal = journal

    def close(self):
        self.journal.append('file')


def a_dataset_that_vaex_opened(journal, leaves_a_file=True):
    """What the _open of vaex leaves in a dataset: a file object in file_map, and the HDF5 file of h5py."""
    dataset = types.SimpleNamespace(path='file.h5', file_map={})

    def original_open(self, path):
        if leaves_a_file:
            self.file_map[self.path] = RecordingFile(journal)
        self.h5file = types.SimpleNamespace(close=lambda: journal.append('h5py'))
    return dataset, original_open


def test_the_file_object_of_h5py_is_closed_after_h5py_closes_the_hdf5_file(monkeypatch):
    journal = []
    dataset, original_open = a_dataset_that_vaex_opened(journal)
    monkeypatch.setattr(_patch_vaex, '_original_hdf5_open', original_open)
    monkeypatch.setattr(vaex.file, 'is_file_object', lambda path: False)
    _patch_vaex.new_hdf5_open(dataset, 'file.h5')
    dataset.h5file.close()
    dataset.h5file.close()
    assert journal == ['h5py', 'file', 'h5py', 'file']


@pytest.mark.parametrize("given_by_the_caller, leaves_a_file", [(True, True), (False, False)])
def test_a_file_of_the_caller_and_a_vaex_that_moved_are_left_alone(given_by_the_caller, leaves_a_file, monkeypatch):
    journal = []
    dataset, original_open = a_dataset_that_vaex_opened(journal, leaves_a_file)
    monkeypatch.setattr(_patch_vaex, '_original_hdf5_open', original_open)
    monkeypatch.setattr(vaex.file, 'is_file_object', lambda path: given_by_the_caller)
    _patch_vaex.new_hdf5_open(dataset, 'file.h5')
    dataset.h5file.close()
    assert journal == ['h5py']


# -- what vaex does without the patch: why it is made ---------------------------------
# These two are expected to fail: they say what vaex alone does not do, and the patch makes up for. When one of them passes
# (strict: the suite fails), vaex does it by itself and the patch, or the part of it, is no longer needed: see _patch_vaex.py.

LOADER_REASON = ("vaex 4.19.0 and later take a group to hold its columns in a subgroup 'columns': the catalogues have them in "
                 "the group itself, which the patched loader opens")
OPEN_MANY_REASON = ("the open_many of vaex takes the list of files and nothing else: the package gives it the group, "
                    "which the patched function passes on to vaex.open")


@pytest.mark.xfail(strict=True, raises=KeyError, reason=LOADER_REASON)
def test_vaex_alone_opens_a_group_that_holds_its_columns_itself(tmp_path, monkeypatch):
    monkeypatch.setattr(vaex.hdf5.dataset.Hdf5MemoryMapped, '_load', _patch_vaex._original_hdf5_load)
    df = vaex.open(ananke_file(tmp_path, 'alone.h5'), group=STARCATALOG_GROUP)
    assert values(df, 'x') == [1.0, 2.0, 3.0]
    df.close()


@pytest.mark.xfail(strict=True, raises=TypeError, reason=OPEN_MANY_REASON)
def test_vaex_alone_opens_a_list_of_files_by_their_group(tmp_path):
    first, second = ananke_file(tmp_path, 'a.h5', stars=(1.0, 2.0)), ananke_file(tmp_path, 'b.h5', stars=(3.0,))
    df = _patch_vaex._original_open_many([str(first), str(second)], group=STARCATALOG_GROUP)
    assert values(df, 'x') == [1.0, 2.0, 3.0]
    df.close()


if __name__ == '__main__':
    pass
