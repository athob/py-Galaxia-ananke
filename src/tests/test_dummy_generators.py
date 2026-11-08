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

from ..galaxia_ananke import Input
from ..galaxia_ananke._builtin_utils import confirm_equal_length_arrays_in_dict


def test_dummy_particles_input():
    n_parts = random.randint(10**2, 10**3)
    p = Input.make_dummy_particles_input(n_parts)
    assert set(p.keys()) == Input.all_possible_keys_in_particles
    confirm_equal_length_arrays_in_dict(p)
    assert all(len(array) == n_parts for array in p.values())

def test_dummy_kernels_input():
    n_parts = random.randint(10**2, 10**3)
    kernels = Input.make_dummy_kernels_input(n_parts)
    assert kernels.shape == (n_parts, 2)   # phase space kernels: one for the positions, one for the velocities
    assert (kernels > 0).all()             # they are exponentials


if __name__ == '__main__':
    pass
