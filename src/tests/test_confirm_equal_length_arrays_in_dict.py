#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 (AGPL-3.0).
# 
# The full copyright notice, including terms governing use, modification,
# and redistribution, is contained in the files LICENSE and COPYRIGHT,
# which can be found at the root of the source code distribution tree:
# - LICENSE <https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE>
# - COPYRIGHT <https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT>
#
import pytest
import random

from ..galaxia_ananke._builtin_utils import confirm_equal_length_arrays_in_dict


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
