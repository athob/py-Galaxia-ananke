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
Tests of the compact metadata, which Output uses to keep the logs of the backend (one per generation of the survey)
in one attribute of its HDF5 files without repeating what they have in common: create_compact_metadata,
get_majority_voting_ordered_common_keys, reconstruct_log and metadata_dicts_to_consolidated_json.
"""
import json
import random

import pytest

from ..galaxia_ananke._builtin_utils import (create_compact_metadata, get_majority_voting_ordered_common_keys,
                                             metadata_dicts_to_consolidated_json, reconstruct_log)


# -- the compact form ---------------------------------------------------------------

def test_nothing_makes_an_empty_compact_metadata():
    assert create_compact_metadata([]) == {}
    assert get_majority_voting_ordered_common_keys([]) == []


def test_what_all_the_dictionaries_share_is_kept_once_and_the_rest_varies():
    compact = create_compact_metadata([{'a': 1, 'b': 2, 'c': 3}, {'a': 1, 'b': 5, 'c': 3}], "logs")
    assert compact == {
        'metadata_name': 'logs',
        'schema_version': '1.0',
        'common_values': {'a': 1, 'c': 3},
        'varying_keys': ['b'],
        'metadata_count': 2,
        'variations': [{'b': 2}, {'b': 5}],
    }


def test_the_name_is_empty_by_default():
    assert create_compact_metadata([{'a': 1}])['metadata_name'] == ""


def test_a_key_that_not_all_the_dictionaries_have_varies():
    compact = create_compact_metadata([{'a': 1, 'x': 9}, {'a': 1}])
    assert compact['common_values'] == {'a': 1}
    assert compact['varying_keys'] == ['x']
    assert compact['variations'] == [{'x': 9}, {}]


def test_one_dictionary_is_all_common():
    compact = create_compact_metadata([{'a': 1, 'b': 'text'}])
    assert (compact['common_values'], compact['varying_keys'], compact['variations']) == ({'a': 1, 'b': 'text'}, [], [{}])
    assert compact['metadata_count'] == 1


def test_equal_dictionaries_have_no_variation():
    compact = create_compact_metadata([{'a': 1}, {'a': 1}, {'a': 1}])
    assert compact['variations'] == [{}, {}, {}] and compact['varying_keys'] == []


def test_the_keys_follow_the_order_most_dictionaries_give_them():
    compact = create_compact_metadata([{'a': 1, 'b': 2, 'c': 3}, {'a': 1, 'b': 7, 'c': 3}, {'c': 3, 'a': 1, 'b': 8}])
    assert list(compact['common_values']) + compact['varying_keys'] == ['a', 'c', 'b']      # b varies, so it comes after the common keys
    assert get_majority_voting_ordered_common_keys([{'a': 1, 'b': 2, 'c': 3}, {'a': 1, 'b': 2, 'c': 3}, {'c': 3, 'a': 1, 'b': 2}]) == ['a', 'b', 'c']


# -- the order of the keys ----------------------------------------------------------

def test_each_key_goes_where_most_dictionaries_put_it():
    dictionaries = [{'a': 1, 'b': 2, 'c': 3}, {'a': 1, 'b': 2, 'c': 3}, {'c': 3, 'a': 1, 'b': 2}]
    assert get_majority_voting_ordered_common_keys(dictionaries) == ['a', 'b', 'c']


def test_keys_that_tie_are_ordered_by_their_first_appearance():
    assert get_majority_voting_ordered_common_keys([{'a': 1}, {'b': 2, 'a': 1}]) == ['a', 'b']
    assert get_majority_voting_ordered_common_keys([{'b': 2, 'a': 1}, {'a': 1}]) == ['b', 'a']


def test_the_first_appearance_decides_between_keys_in_the_same_place():
    # a and b are both first in their dictionary: the one that was seen first is first, in either order
    assert get_majority_voting_ordered_common_keys([{'a': 1}, {'b': 2}]) == ['a', 'b']
    assert get_majority_voting_ordered_common_keys([{'b': 2}, {'a': 1}]) == ['b', 'a']


def test_every_key_of_every_dictionary_is_listed_once():
    keys = get_majority_voting_ordered_common_keys([{'a': 1, 'b': 2}, {'c': 3, 'a': 1}, {'d': 4}])
    assert sorted(keys) == ['a', 'b', 'c', 'd'] and len(keys) == 4


# -- going back ---------------------------------------------------------------------

def test_each_dictionary_is_reconstructed_from_the_compact_form():
    dictionaries = [{'a': 1, 'b': 2, 'x': 9}, {'a': 1, 'b': 5}]
    compact = create_compact_metadata(dictionaries)
    assert [reconstruct_log(compact, i) for i in range(2)] == dictionaries


def test_reconstructing_does_not_change_the_compact_form():
    compact = create_compact_metadata([{'a': 1, 'b': 2}, {'a': 1, 'b': 5}])
    before = json.dumps(compact)
    reconstruct_log(compact, 1)
    assert json.dumps(compact) == before


def test_an_index_past_the_last_dictionary_is_an_index_error():
    compact = create_compact_metadata([{'a': 1}, {'a': 2}])
    with pytest.raises(IndexError, match="Metadata index 2 out of range"):
        reconstruct_log(compact, 2)


def test_whatever_the_dictionaries_the_round_trip_gives_them_back():
    for _ in range(50):
        keys = random.sample(['a', 'b', 'c', 'd', 'e', 'f'], k=random.randint(1, 6))
        dictionaries = [{key: random.choice([1, 2, 'x'])
                         for key in random.sample(keys, k=random.randint(1, len(keys)))}
                        for _ in range(random.randint(1, 5))]
        compact = create_compact_metadata(dictionaries)
        assert compact['metadata_count'] == len(dictionaries)
        assert [reconstruct_log(compact, i) for i in range(len(dictionaries))] == dictionaries


# -- the json -----------------------------------------------------------------------

def test_the_json_is_the_compact_form_written_with_an_indent_of_two():
    dictionaries = [{'a': 1, 'b': 2}, {'a': 1, 'b': 5}]
    text = metadata_dicts_to_consolidated_json(dictionaries, "logs")
    assert json.loads(text) == create_compact_metadata(dictionaries, "logs")
    assert text == json.dumps(create_compact_metadata(dictionaries, "logs"), indent=2)


def test_the_json_takes_its_arguments_by_keyword_too():
    dictionaries = [{'a': 1}]
    assert (metadata_dicts_to_consolidated_json(metadata_dicts=dictionaries, metadata_name="logs")
            == metadata_dicts_to_consolidated_json(dictionaries, "logs"))


def test_the_json_of_nothing_is_an_empty_object():
    assert metadata_dicts_to_consolidated_json([]) == "{}"


if __name__ == '__main__':
    pass
