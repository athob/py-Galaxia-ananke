#!/usr/bin/env python
import pytest
import random

from ..Galaxia_ananke._builtin_utils import confirm_equal_length_arrays_in_dict


def make_random_list(length):
    return [random.randrange(100*length) for _ in range(length)]


def test():
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


if __name__ == '__main__':
    pass
