#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
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
import ast
import re

import pytest

from ..galaxia_ananke._builtin_utils import compare_given_and_required


GIVEN = ['A', 'B', 'C', 'D']
REQUIRED = ['A', 'C']
OPTIONAL = ['B', 'D']
EXTRA = ['E']


def keys_in(message, word):
    """The set that the message writes after the word (misses {...} or misincludes {...})."""
    return ast.literal_eval(re.search(rf"{word} (\{{[^}}]*\}})", message).group(1))


def test_given_keys_are_the_required_ones_and_any_optional_ones():
    compare_given_and_required(GIVEN, REQUIRED, OPTIONAL)
    compare_given_and_required(GIVEN, REQUIRED, OPTIONAL + EXTRA)
    compare_given_and_required(REQUIRED, REQUIRED, OPTIONAL)            # the optional ones may be absent
    compare_given_and_required([], [], OPTIONAL)


def test_a_wrong_set_of_keys_is_a_value_error():
    with pytest.raises(ValueError):
        compare_given_and_required(GIVEN)                               # nothing is accepted by default
    with pytest.raises(ValueError):
        compare_given_and_required(GIVEN, REQUIRED)                     # the optional ones are not allowed
    with pytest.raises(ValueError):
        compare_given_and_required(GIVEN, REQUIRED + EXTRA, OPTIONAL)   # E is required and missing
    with pytest.raises(ValueError):
        compare_given_and_required(GIVEN + EXTRA, REQUIRED, OPTIONAL)   # E is neither


def test_what_is_missing_is_named():
    with pytest.raises(ValueError) as error:
        compare_given_and_required(['A'], REQUIRED, OPTIONAL)
    assert str(error.value) == "Given particle data covers wrong set of keys: misses {'C'}"


def test_what_is_not_allowed_is_named():
    with pytest.raises(ValueError) as error:
        compare_given_and_required(REQUIRED + EXTRA, REQUIRED, OPTIONAL)
    assert str(error.value) == "Given particle data covers wrong set of keys: misincludes {'E'}"


def test_the_optional_keys_are_never_named_as_not_allowed():
    with pytest.raises(ValueError) as error:
        compare_given_and_required(GIVEN + EXTRA, REQUIRED, OPTIONAL)
    assert str(error.value) == "Given particle data covers wrong set of keys: misincludes {'E'}"


def test_what_is_missing_and_what_is_not_allowed_are_both_named():
    with pytest.raises(ValueError) as error:
        compare_given_and_required(['A', 'E'], REQUIRED, OPTIONAL)
    assert str(error.value) == "Given particle data covers wrong set of keys: misses {'C'} & misincludes {'E'}"


def test_several_keys_are_named_together():
    with pytest.raises(ValueError) as error:
        compare_given_and_required(['X', 'Y'], ['A', 'B', 'C'])
    message = str(error.value)
    assert keys_in(message, "misses") == {'A', 'B', 'C'}
    assert keys_in(message, "misincludes") == {'X', 'Y'}


def test_the_error_message_can_be_given():
    with pytest.raises(ValueError) as error:
        compare_given_and_required(['A'], REQUIRED, error_message="Wrong keys")
    assert str(error.value) == "Wrong keys: misses {'C'}"


def test_the_keys_can_be_given_as_any_iterables():
    compare_given_and_required({'A': 1, 'C': 2}, ('A', 'C'), {'B'})     # the keys of a dictionary, a tuple, a set
    compare_given_and_required((k for k in GIVEN), iter(REQUIRED), frozenset(OPTIONAL))


def test_the_arguments_are_not_changed():
    given, required, optional = list(GIVEN), list(REQUIRED), list(OPTIONAL)
    compare_given_and_required(given, required, optional)
    assert (given, required, optional) == (GIVEN, REQUIRED, OPTIONAL)


if __name__ == '__main__':
    pass
