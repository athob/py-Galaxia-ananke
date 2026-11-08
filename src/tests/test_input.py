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
Tests of Input, which holds the particles and the kernels that Galaxia takes as input: how it is made and what it refuses, what
it completes, and what it says about itself. Its hash, the files it writes and how it is prepared for the backend are in
test_input_files.py. The inputs here are the smallest valid ones (see utils.minimal_particles), so that what is asserted is exact.
"""
import collections
import pathlib
import warnings

import numpy as np
import pytest

from ..galaxia_ananke import Input
from ..galaxia_ananke._defaults import GALAXIA_TMP
from .utils import minimal_kernels, minimal_particles

REQUIRED = {'pos3', 'vel3', 'massinit', 'age', 'feh'}
OPTIONAL = {'parentid', 'partitionid', 'id', 'dform', 'helium', 'carbon', 'nitrogen', 'oxygen', 'neon',
            'magnesium', 'silicon', 'sulphur', 'calcium', 'alpha'}


def make(**kwargs):
    return Input(minimal_particles(), minimal_kernels(), **kwargs)


# -- what a particle is made of -----------------------------------------------------

def test_the_required_and_the_optional_keys_of_the_particles():
    assert Input._required_keys_in_particles == REQUIRED
    assert Input._optional_keys_in_particles == OPTIONAL
    assert Input.all_possible_keys_in_particles == REQUIRED | OPTIONAL


def test_the_present_day_mass_is_not_a_key_of_the_particles():
    assert Input._mass == 'mass' and 'mass' not in Input.all_possible_keys_in_particles


def test_the_description_of_the_particles_names_every_key():
    description = Input.particles_dictionary_description
    for key in REQUIRED | OPTIONAL:
        assert f"via key ``{key}``" in description


def test_the_documentation_of_the_constructor_has_what_it_is_made_of_filled_in():
    documentation = Input.__init__.__doc__
    assert "{" not in documentation.replace("{}", "")                  # no placeholder left
    assert all(f"via key ``{key}``" in documentation for key in REQUIRED | OPTIONAL)


# -- the ways to make an Input ------------------------------------------------------

def test_an_input_is_made_from_the_particles_and_the_kernels():
    inp = make()
    assert inp.length == 2
    assert inp.particles['massinit'].tolist() == [1000.0, 2000.0]
    assert inp.kernels.tolist() == [[1.0, 2.0], [3.0, 4.0]]


@pytest.mark.parametrize("arguments, keywords", [
    ((), {}),
    ((), {'pname': 'a.ebf'}),                                              # a pname without its kname
    ((), {'kname': 'a_d6n64_den.ebf'}),
    ((), {'particles': minimal_particles(), 'kernels': minimal_kernels()}),    # the particles and the kernels are positional
])
def test_a_wrong_signature_is_a_value_error_that_says_so(arguments, keywords):
    with pytest.raises(ValueError, match="Wrong signature"):
        Input(*arguments, **keywords)


@pytest.mark.xfail(raises=RuntimeError, reason="a number of arguments that is not two reaches a bare `raise`, which fails "
                   "with 'No active exception to reraise' instead of saying what is wrong")
@pytest.mark.parametrize("arguments", [(1,), (1, 2, 3)])
def test_a_wrong_number_of_arguments_is_a_type_or_value_error(arguments):
    with pytest.raises((TypeError, ValueError)):
        Input(*arguments)


def test_existing_files_have_to_be_in_the_same_directory():
    with pytest.raises(ValueError, match="need to share the same parent directory"):
        Input(pname='one/sim_abcdef0.ebf', kname='another/sim_abcdef0_d6n64_den.ebf')


# -- what is refused ----------------------------------------------------------------

@pytest.mark.parametrize("key", sorted(REQUIRED))
def test_a_required_key_that_is_missing_is_named(key):
    particles = minimal_particles()
    del particles[key]
    with pytest.raises(ValueError, match=f"Given particle data covers wrong set of keys: misses {{'{key}'}}"):
        Input(particles, minimal_kernels())


def test_a_key_that_is_not_allowed_is_named():
    particles = {**minimal_particles(), 'luminosity': np.array([1.0, 2.0])}
    with pytest.raises(ValueError, match="misincludes {'luminosity'}"):
        Input(particles, minimal_kernels())


def test_the_present_day_mass_next_to_the_initial_mass_is_refused():
    particles = {**minimal_particles(), 'mass': np.array([700.0, 1400.0])}
    with pytest.raises(ValueError, match="misincludes {'mass'}"):
        Input(particles, minimal_kernels())


@pytest.mark.xfail(raises=ModuleNotFoundError, reason="the explanation of the change of the mass keys starts with an import "
                   "of __metadata__ that is not a relative import, so a script that gives only `mass` gets "
                   "'No module named __metadata__' instead of the explanation")
def test_the_present_day_mass_alone_is_refused_with_an_explanation_of_the_change():
    particles = minimal_particles()
    particles['mass'] = particles.pop('massinit')
    with pytest.raises(KeyError, match="BACKWARD INCOMPATIBILITY"):
        Input(particles, minimal_kernels())


@pytest.mark.parametrize("key", sorted(REQUIRED - {'massinit'}))
def test_an_array_that_is_not_as_long_as_the_initial_mass_is_named(key):
    particles = minimal_particles()
    particles[key] = particles[key][:1]
    with pytest.raises(ValueError, match=f"{{'{key}'}} in the provided particles input dictionary does not have the same length as property massinit"):
        Input(particles, minimal_kernels())


def test_kernels_with_three_columns_are_refused():
    with pytest.raises(ValueError, match=r"Given 2D-array kernels is \(Nx3\) but should be either \(Nx1\)"):
        Input(minimal_particles(), np.ones((2, 3)))


def test_kernels_that_are_not_one_or_two_dimensional_are_refused():
    with pytest.raises(ValueError, match="Given kernels ndim is 3 but should be either 1"):
        Input(minimal_particles(), np.ones((2, 2, 2)))


@pytest.mark.xfail(raises=pytest.fail.Exception, reason="nothing checks that there are as many kernels as particles, as the "
                   "documentation says there must be: the backend is given files that do not match")
def test_kernels_that_are_not_as_many_as_the_particles_are_refused():
    with pytest.raises(ValueError):
        Input(minimal_particles(), np.ones((3, 2)))


@pytest.mark.xfail(raises=AttributeError, reason="the kernels are used as an array whatever the documentation says (array_like): "
                   "a list has no ndim")
def test_kernels_given_as_a_list_are_accepted():
    assert Input(minimal_particles(), minimal_kernels().tolist()).hdim == 6


@pytest.mark.xfail(raises=AttributeError, reason="the particles are used as arrays whatever the documentation says (array_like): "
                   "a list has no shape")
def test_particles_given_as_lists_are_accepted():
    particles = {key: value.tolist() for key, value in minimal_particles().items()}
    assert Input(particles, minimal_kernels()).length == 2


# -- the kernels --------------------------------------------------------------------

def test_kernels_in_one_column_become_kernels_in_position_space():
    inp = Input(minimal_particles(), np.array([[1.0], [2.0]]))
    assert inp.kernels.tolist() == [1.0, 2.0] and inp.hdim == 3


def test_one_kernel_for_each_particle_is_kernels_in_position_space():
    inp = Input(minimal_particles(), np.array([1.0, 2.0]))
    assert inp.kernels.shape == (2,) and inp.hdim == 3


def test_two_kernels_for_each_particle_are_kernels_in_phase_space():
    assert make().hdim == 6


def test_the_kernels_are_scaled_by_the_factor_in_what_is_sent_to_the_backend_only():
    inp = make(k_factor=2.0)
    assert inp.kernels.tolist() == [[1.0, 2.0], [3.0, 4.0]]                # as given
    sent = dict(inp.parent_particles_data_generator)
    assert set(sent) == {'pos3', 'vel3', 'massinit', 'age', 'feh', 'h_cubic'}      # the parent index is not sent: it is what is made
    assert sent['h_cubic'].tolist() == [[2.0, 4.0], [6.0, 8.0]]            # as scaled


# -- what is completed --------------------------------------------------------------

def test_the_keys_are_in_lexical_order_and_the_ones_that_are_added_come_last():
    assert list(make().keys()) == ['age', 'feh', 'massinit', 'pos3', 'vel3', 'parentid', 'partitionid']


def test_the_parent_of_each_particle_is_itself_by_default():
    inp = make()
    assert inp.particle_parentids.tolist() == [0, 1]
    assert inp.particle_parentids.dtype.kind == 'i'


def test_the_partition_of_each_particle_is_the_first_by_default():
    inp = make()
    assert inp.particles['partitionid'].tolist() == [0, 0] and inp.particles['partitionid'].dtype.kind == 'i'


def test_a_parent_and_a_partition_that_are_given_are_kept():
    particles = {**minimal_particles(), 'parentid': np.array([7, 8]), 'partitionid': np.array([1, 0])}
    inp = Input(particles, minimal_kernels())
    assert inp.particle_parentids.tolist() == [7, 8] and inp.particles['partitionid'].tolist() == [1, 0]
    assert list(inp.keys()) == ['age', 'feh', 'massinit', 'parentid', 'partitionid', 'pos3', 'vel3']      # all in lexical order


def test_the_optional_keys_that_are_there_are_listed():
    particles = {**minimal_particles(), 'dform': np.array([0.0, 1.0]), 'alpha': np.array([0.1, 0.2])}
    assert set(Input(particles, minimal_kernels()).optional_keys()) == {'dform', 'alpha', 'parentid', 'partitionid'}
    assert set(make().optional_keys()) == {'parentid', 'partitionid'}


def test_the_particles_are_an_ordered_dictionary_that_is_not_the_one_given():
    inp = make()
    assert isinstance(inp.particles, collections.OrderedDict)
    given = minimal_particles()
    Input(given, minimal_kernels())
    assert list(given) == list(minimal_particles())                        # what was given was not changed (nor completed)
    assert inp.keys() == inp.particles.keys()


def test_the_arrays_are_the_ones_that_were_given_and_not_copies():
    # How it is today: changing an array after the Input was made changes the Input
    particles, kernels = minimal_particles(), minimal_kernels()
    inp = Input(particles, kernels)
    assert inp.particles['age'] is particles['age'] and inp.kernels is kernels
    particles['age'][0] = -999.0
    assert inp.particles['age'][0] == -999.0


# -- what an Input says about itself ------------------------------------------------

def test_the_defaults_of_an_input():
    inp = make()
    assert (inp.name, inp.k_factor, inp.caching, inp.append_hash) == ('sim', 1.0, False, True)
    assert inp._input_dir == pathlib.Path(GALAXIA_TMP)


@pytest.mark.xfail(raises=AssertionError, reason="the default number of neighbours is the name of the tag that holds it ('nres'), "
                   "not its value, 64: a survey of an Input that was not given ngb is written with `nres nres`, which "
                   "the backend cannot read (make_survey_from_particles always gives ngb, which hides it)")
def test_the_default_number_of_neighbours_is_64():
    assert make().ngb == 64


def test_what_is_given_is_what_the_input_says():
    inp = make(name='mine', ngb=32, k_factor=0.5, append_hash=False)
    assert (inp.name, inp.ngb, inp.k_factor, inp.append_hash) == ('mine', 32, 0.5, False)


@pytest.mark.parametrize("directory", ["somewhere", pathlib.Path("somewhere")])
def test_the_input_directory_is_a_path(directory):
    assert make(input_dir=directory)._input_dir == pathlib.Path("somewhere")


def test_caching_is_experimental_and_says_so():
    with pytest.warns(DeprecationWarning, match="caching mode"):
        inp = make(caching=True)
    assert inp.caching is True
    with pytest.warns(DeprecationWarning, match="caching mode"):
        inp.caching = True


def test_without_caching_nothing_is_said():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        inp = make(caching=False)
        inp.caching = False


# -- the dummy inputs ---------------------------------------------------------------

def test_the_dummy_inputs_are_a_valid_input():
    inp = Input(Input.make_dummy_particles_input(50), Input.make_dummy_kernels_input(50))
    assert inp.length == 50 and inp.hdim == 6
    assert set(inp.keys()) == Input.all_possible_keys_in_particles


def test_the_dummy_inputs_do_not_depend_on_anything_but_the_random_state():
    np.random.seed(5)
    first = Input.make_dummy_particles_input(10), Input.make_dummy_kernels_input(10)
    np.random.seed(5)
    second = Input.make_dummy_particles_input(10), Input.make_dummy_kernels_input(10)
    assert all(np.array_equal(first[0][key], second[0][key]) for key in first[0])
    assert np.array_equal(first[1], second[1])


if __name__ == '__main__':
    pass
