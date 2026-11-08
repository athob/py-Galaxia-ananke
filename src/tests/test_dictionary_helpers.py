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
Tests of the small helpers that deal with dictionaries: common_entries and lexicalorder_dict.
"""
from collections import OrderedDict

from ..galaxia_ananke._builtin_utils import common_entries, lexicalorder_dict


def test_common_entries_zips_the_keys_that_all_the_dictionaries_have():
    assert sorted(common_entries({'a': 1, 'b': 2}, {'a': 3, 'c': 4})) == [('a', 1, 3)]
    assert sorted(common_entries({'a': 1, 'b': 2}, {'b': 3, 'a': 4}, {'a': 5, 'b': 6, 'c': 7})) == [('a', 1, 4, 5), ('b', 2, 3, 6)]


def test_common_entries_of_one_dictionary_are_its_items():
    assert sorted(common_entries({'a': 1, 'b': 2})) == [('a', 1), ('b', 2)]


def test_common_entries_of_nothing_or_of_disjoint_dictionaries_is_nothing():
    assert list(common_entries()) == []
    assert list(common_entries({'a': 1}, {'b': 2})) == []
    assert list(common_entries({}, {'b': 2})) == []


def test_common_entries_is_lazy():
    entries = common_entries({'a': 1}, {'a': 2})
    assert iter(entries) is entries


def test_lexicalorder_dict_sorts_the_keys():
    ordered = lexicalorder_dict({'b': 1, 'c': 2, 'a': 3})
    assert isinstance(ordered, OrderedDict)
    assert list(ordered.items()) == [('a', 3), ('b', 1), ('c', 2)]


def test_lexicalorder_dict_sorts_by_characters_not_by_meaning():
    assert list(lexicalorder_dict({'b': 1, 'B': 2, 'a': 3, '1': 4})) == ['1', 'B', 'a', 'b']


def test_lexicalorder_dict_leaves_its_argument_alone_and_keeps_the_values_themselves():
    value = [1, 2]
    original = {'b': value, 'a': 0}
    ordered = lexicalorder_dict(original)
    assert list(original) == ['b', 'a']
    assert ordered['b'] is value


def test_lexicalorder_dict_of_an_empty_dictionary_is_empty():
    assert lexicalorder_dict({}) == OrderedDict()


if __name__ == '__main__':
    pass
