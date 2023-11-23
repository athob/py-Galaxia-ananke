#!/usr/bin/env python
import random

from ..Galaxia_ananke import Input
from ..Galaxia_ananke.utils import confirm_equal_length_arrays_in_dict


n_parts = random.randint(10**3,10**4)


def test_dummy_particles_input():
    p = Input.make_dummy_particles_input(n_parts)
    assert set(p.keys()) == Input.all_possible_keys_in_particles
    confirm_equal_length_arrays_in_dict(p)

def test_dummy_densities_input():
    rho = Input.make_dummy_densities_input(n_parts)
    assert len(rho) == 2
    assert len(rho[0]) == len(rho[1])


if __name__ == '__main__':
    pass
