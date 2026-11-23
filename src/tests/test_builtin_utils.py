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
format_docstring, which fills the docstrings of the package with its defaults, is tested here too.
"""
import subprocess
import sys
import textwrap
from operator import itemgetter
from unittest import mock

import numpy as np
import pytest

from ..galaxia_ananke._builtin_utils import (classproperty, format_docstring, get_version_of_command, import_source_file,
                                            in_description_order, tuple_itemgetter)
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


# -- the docstrings that are filled with the defaults ---------------------------------

def test_the_placeholders_of_a_docstring_are_filled_with_the_values():
    def function():
        """The default is {default}, and the name {name}."""
    format_docstring(function, default=3, name='sim')
    assert function.__doc__ == "The default is 3, and the name sim."


def test_the_docstring_of_a_method_is_filled_through_its_function():
    class Example:
        @classmethod
        def made(cls):
            """Made with {default}."""
        def called(self):
            """Called with {default}."""
    format_docstring(Example.made, default=1)
    format_docstring(Example.called, default=2)
    assert Example.made.__doc__ == "Made with 1." and Example.called.__doc__ == "Called with 2."


def test_a_value_that_the_docstring_asks_for_and_is_not_given_is_a_key_error():
    def function():
        """The default is {missing}."""
    with pytest.raises(KeyError):
        format_docstring(function)


def test_a_function_that_has_no_docstring_is_left_without_it():
    def function():
        pass
    function.__doc__ = None
    format_docstring(function, anything=1)
    assert function.__doc__ is None


def test_python_that_removes_the_docstrings_can_fill_them():
    # python -OO removes every docstring, so the package is imported where each is None: the helper does nothing there, and
    # does not raise as `function.__doc__.format(...)` does. It is loaded from its file, in a process of its own.
    helper_file = sys.modules[format_docstring.__module__].__file__
    code = textwrap.dedent("""
        import importlib.util, sys
        spec = importlib.util.spec_from_file_location("helper_under_OO", sys.argv[1])
        module = importlib.util.module_from_spec(spec)
        sys.modules["helper_under_OO"] = module
        spec.loader.exec_module(module)
        def function():
            \"\"\"The default is {default}.\"\"\"
        if function.__doc__ is not None:
            raise SystemExit("the docstrings were not removed: this is not python -OO")
        module.format_docstring(function, default=1)
        print("left without docstring" if function.__doc__ is None else "made a docstring")
    """)
    result = subprocess.run([sys.executable, '-OO', '-c', code, helper_file], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "left without docstring"


# -- the order of the properties that a docstring lists -------------------------------

def test_properties_are_put_in_the_alphabetical_order_of_their_descriptions():
    properties = {('b', 'Second'), (('x', 'y'), 'First'), ('a', 'Third')}
    assert in_description_order(properties) == [(('x', 'y'), 'First'), ('b', 'Second'), ('a', 'Third')]


def test_the_order_does_not_depend_on_the_order_in_which_the_properties_come():
    pairs = [(f'key{i}', f'Description {i % 3}') for i in range(9)]
    expected = in_description_order(pairs)
    for shuffled in (pairs[::-1], pairs[3:] + pairs[:3], set(pairs), frozenset(pairs)):
        assert in_description_order(shuffled) == expected


def test_keys_that_are_strings_and_tuples_are_put_in_order_together_by_how_they_are_written():
    # the same description: the keys are compared as the text that the docstring shows, not a string with a tuple
    assert in_description_order([('a', 'Same'), (('b', 'a'), 'Same')]) == [(('b', 'a'), 'Same'), ('a', 'Same')]
