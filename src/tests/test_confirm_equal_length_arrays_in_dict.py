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
import random

import pytest

from ..galaxia_ananke._builtin_utils import confirm_equal_length_arrays_in_dict


def make_random_list(length):
    return [random.randrange(100*length) for _ in range(length)]


def test_the_arrays_of_a_dictionary_have_the_length_of_the_control():
    length = 10
    control = 'control'
    key = 'prop{}'
    good_dict = {}
    bad_dict = {}
    good_dict[control] = bad_dict[control] = make_random_list(length)
    good_dict[key.format(1)] = make_random_list(length)
    good_dict[key.format(2)] = make_random_list(length)
    bad_dict[key.format(1)] = make_random_list(length)
    bad_dict[key.format(2)] = make_random_list(length+1)
    confirm_equal_length_arrays_in_dict(good_dict, control=control)
    with pytest.raises(ValueError):
        confirm_equal_length_arrays_in_dict(bad_dict, control=control)


def test_the_first_key_is_the_control_by_default():
    confirm_equal_length_arrays_in_dict({'first': [1, 2], 'second': [3, 4]})
    with pytest.raises(ValueError) as error:
        confirm_equal_length_arrays_in_dict({'first': [1, 2], 'second': [3]})
    assert str(error.value).endswith("not have the same length as property first.")


def test_an_empty_dictionary_and_a_single_array_are_fine():
    confirm_equal_length_arrays_in_dict({})
    confirm_equal_length_arrays_in_dict({'only': [1, 2, 3]})


def test_one_wrong_array_is_named_in_the_singular():
    with pytest.raises(ValueError) as error:
        confirm_equal_length_arrays_in_dict({'age': [1, 2], 'feh': [1]}, control='age')
    assert str(error.value) == ("Array representing property {'feh'} in the provided input dictionary "
                                "does not have the same length as property age.")


def test_several_wrong_arrays_are_named_in_the_plural():
    with pytest.raises(ValueError) as error:
        confirm_equal_length_arrays_in_dict({'age': [1, 2], 'feh': [1], 'mass': [1, 2, 3]}, control='age')
    message = str(error.value)
    assert message.startswith("Arrays representing properties {")
    assert message.endswith("in the provided input dictionary do not have the same length as property age.")
    assert "'feh'" in message and "'mass'" in message and "'age'" not in message.split("}")[0]


def test_the_dictionary_can_be_named_in_the_error():
    with pytest.raises(ValueError) as error:
        confirm_equal_length_arrays_in_dict({'age': [1, 2], 'feh': [1]}, control='age', error_message_dict_name='particles')
    assert "in the provided particles input dictionary" in str(error.value)


def test_any_sized_values_are_compared_by_their_length():
    confirm_equal_length_arrays_in_dict({'list': [1, 2], 'tuple': (3, 4), 'text': 'ab', 'range': range(2)})


def test_a_control_that_is_not_a_key_is_a_key_error():
    with pytest.raises(KeyError):
        confirm_equal_length_arrays_in_dict({'a': [1], 'b': [1]}, control='nothing')


if __name__ == '__main__':
    pass
