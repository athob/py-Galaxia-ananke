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
A survey run with the real backend, for the integration tests: a few hundred stars made from the dummy particles of the package by
the compiled Galaxia, once for all the tests, in a temporary directory (the files of the input of the backend that have to be in
the cache of the installation are removed afterwards). The numbers are not worked out here: what the tests check is what has to hold
whatever the numbers are.
"""
import shutil
import types

import numpy as np
import pytest

from ..galaxia_ananke import DEFAULT_PSYS, Input, Survey
from ..galaxia_ananke._constants import GALAXIA_FILENAMES, GALAXIA_NBODY1

PARTICLES = 1000
SEED = 2026
LIMITS = dict(fsample=0.01, abs_mag_lim_hi=6, app_mag_lim_hi=22)      # a few hundred stars, in a few seconds


def forget_input(input):
    """Removes from the cache of the installation what a run of the backend needed of this input (the links to its files)."""
    try:
        name_hash = input.name_hash
    except AttributeError:                      # the run did not get as far as naming its input: nothing was linked
        return
    shutil.rmtree(GALAXIA_NBODY1 / name_hash, ignore_errors=True)
    (GALAXIA_FILENAMES / f"{name_hash}.txt").unlink(missing_ok=True)


def forget_inputs_named(simname):
    """The same for the inputs of a name, for a run that failed before it could say which they were (the name has to be a test's)."""
    assert simname.startswith('itest'), f"{simname} is not the name that a test gives to its input: it is not for the tests to remove it"
    for entry in GALAXIA_NBODY1.glob(f"{simname}_*"):
        shutil.rmtree(entry, ignore_errors=True)
    for entry in GALAXIA_FILENAMES.glob(f"{simname}_*.txt"):
        entry.unlink(missing_ok=True)


def dummy_input():
    """The same particles and kernels each time."""
    np.random.seed(SEED)
    return Input.make_dummy_particles_input(PARTICLES), Input.make_dummy_kernels_input(PARTICLES)


def real_run(tmp_path_factory, name, photo_sys=DEFAULT_PSYS, **kwargs):
    """Runs a survey like that of real_survey, with other things that the test is about, and cleans up when it is done."""
    directory = tmp_path_factory.mktemp(name)
    (directory / 'in').mkdir()
    (directory / 'out').mkdir()
    particles, kernels = dummy_input()
    input = Input(particles, kernels, input_dir=directory / 'in', name=name, ngb=64)
    survey = Survey(input, photo_sys=photo_sys, surveyname=name, verbose=False)
    try:
        output = survey.make_survey(verbose=False, output_dir=directory / 'out', **{**LIMITS, **kwargs})
        yield types.SimpleNamespace(survey=survey, output=output, particles=particles, kernels=kernels, directory=directory)
    finally:
        forget_input(input)
        forget_inputs_named(name)


@pytest.fixture(scope='session')
def real_generations(tmp_path_factory):
    yield from real_run(tmp_path_factory, 'itest_gens', n_gens=2)


@pytest.fixture(scope='session')
def real_systems(tmp_path_factory):
    yield from real_run(tmp_path_factory, 'itest_systems', photo_sys=[DEFAULT_PSYS, 'padova/LSST'])


@pytest.fixture(scope='session')
def real_survey(tmp_path_factory):
    directory = tmp_path_factory.mktemp('integration')
    particles, kernels = dummy_input()
    (directory / 'in').mkdir()
    (directory / 'out').mkdir()
    input = Input(particles, kernels, input_dir=directory / 'in', name='itest', ngb=64)
    survey = Survey(input, surveyname='itest', verbose=False)
    try:
        output = survey.make_survey(verbose=False, output_dir=directory / 'out', **LIMITS)
        yield types.SimpleNamespace(survey=survey, output=output, particles=particles, kernels=kernels, directory=directory)
    finally:
        forget_input(input)


if __name__ == '__main__':
    pass
