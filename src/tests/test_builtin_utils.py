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
Tests of the small helpers of _builtin_utils that find things and load them: get_version_of_command,
import_source_file and classproperty. The first two are the ones that py-EnBiD-ananke has too, and are tested the same.
"""
import subprocess
import sys
from operator import itemgetter
from unittest import mock

import numpy as np
import pytest

from ..galaxia_ananke._builtin_utils import classproperty, get_version_of_command, import_source_file, tuple_itemgetter
from .utils import in_tmp_wd, fake_commands


@pytest.mark.parametrize('output, version', [
    ('GNU Make 4.3\nBuilt for x86_64-pc-linux-gnu', '4.3'),
    ('g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0', '11.4.0'),          # the first of the numbers: not 22.04
    ('tar (GNU tar) 1.34', '1.34'),
    ("it's a tool's version 5.2", '5.2'),                              # quotes are not a problem
    ('version 2.7.18.1 (build 5.3)', '2.7.18.1'),
    ('release 7 of the tool 2.1', '2.1'),                              # a number without a dot is not a version
])
def test_the_version_is_the_first_dotted_number_of_the_output(output, version):
    with fake_commands({'tool': output}):
        assert get_version_of_command('tool') == version


def test_the_command_is_asked_for_its_version_with_the_version_option(tmp_path):
    tool = tmp_path / "tool"
    tool.write_text('#!/bin/sh\nif [ "$1" = "--version" ]; then echo "tool 3.4"; else exit 2; fi\n')
    tool.chmod(0o755)
    assert get_version_of_command(str(tool)) == '3.4'


def test_an_output_without_a_version_is_an_index_error():
    # Not a friendly error, but the one there is: the callers only expect FileNotFoundError.
    with fake_commands({'tool': 'no number here'}):
        with pytest.raises(IndexError):
            get_version_of_command('tool')


def test_a_command_that_is_not_there_is_a_file_not_found_error():
    with fake_commands({}):
        with pytest.raises(FileNotFoundError):
            get_version_of_command('tool')


def test_a_command_that_fails_raises():
    with fake_commands({'tool': ('1.0', 4)}):                            # it prints a version, then exits with an error
        with pytest.raises(subprocess.CalledProcessError):
            get_version_of_command('tool')


@in_tmp_wd
def test_a_source_file_is_imported_as_the_module_of_that_name():
    with open('mine.py', 'w') as f:
        f.write('VALUE = 42\ndef double(x):\n    return 2 * x\n')
    try:
        module = import_source_file('a_module_made_for_this_test', 'mine.py')
        assert module.VALUE == 42 and module.double(4) == 8
        assert sys.modules['a_module_made_for_this_test'] is module
    finally:
        sys.modules.pop('a_module_made_for_this_test', None)


@in_tmp_wd
def test_importing_a_source_file_again_reads_it_again():
    # Python trusts the bytecode it cached for a file whose size and modification time (to the second) are those
    # it was made from, and the two versions of this file are written within the same second and have the same
    # size: no bytecode is written here, so that what is tested is the function, not that cache.
    try:
        with mock.patch.object(sys, 'dont_write_bytecode', True):
            with open('mine.py', 'w') as f:
                f.write('VALUE = 1\n')
            assert import_source_file('a_module_made_for_this_test', 'mine.py').VALUE == 1
            with open('mine.py', 'w') as f:
                f.write('VALUE = 2\n')
            assert import_source_file('a_module_made_for_this_test', 'mine.py').VALUE == 2
    finally:
        sys.modules.pop('a_module_made_for_this_test', None)


@in_tmp_wd
def test_a_source_file_that_does_not_exist_is_a_file_not_found_error():
    with pytest.raises(FileNotFoundError):
        import_source_file('a_module_made_for_this_test', 'nothing here.py')


def test_a_classproperty_is_read_from_the_class_and_from_an_instance():
    class Thing:
        @classproperty
        def label(cls):
            return cls.__name__.lower()
    assert Thing.label == 'thing'
    assert Thing().label == 'thing'


def test_a_classproperty_sees_the_class_it_is_read_from():
    class Thing:
        @classproperty
        def label(cls):
            return cls.__name__.lower()
    class Special(Thing):
        pass
    assert Special.label == 'special'
    assert Special().label == 'special'


if __name__ == '__main__':
    pass


# -- the items that are always a tuple ----------------------------------------------

@pytest.mark.parametrize("indices", [(1,), (0, 2), (2, 0, 1)], ids=["one", "two", "three"])
def test_the_items_are_always_given_as_a_tuple(indices):
    items = ['a', 'b', 'c']
    assert tuple_itemgetter(*indices)(items) == tuple(items[i] for i in indices)


def test_a_single_index_gives_a_tuple_of_one_where_itemgetter_gives_the_item():
    assert itemgetter(1)(['a', 'b']) == 'b' and tuple_itemgetter(1)(['a', 'b']) == ('b',)


def test_the_indices_can_be_arrays_as_the_groups_of_particles_are():
    ages = np.array([9.1, 8.2, 7.3])
    (group,) = tuple_itemgetter(np.array([2, 0]))(ages)
    assert group.tolist() == [7.3, 9.1]
    first, second = tuple_itemgetter(np.array([2, 0]), np.array([1]))(ages)
    assert first.tolist() == [7.3, 9.1] and second.tolist() == [8.2]


def test_no_index_is_a_type_error():
    with pytest.raises(TypeError):
        tuple_itemgetter()
