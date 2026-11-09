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
Tests of make_survey_from_particles, the function that the package is documented to be used with: that it runs a survey end to end,
that the same input and the same parameters, or the same files, make the same survey, and what it says of a call that it cannot make
sense of. The runs use the compiled backend (and clean up what they leave in the cache of the installation).
"""
import sys

import numpy as np
import pytest

from ..galaxia_ananke import Input, make_survey_from_particles
from .real_backend import LIMITS, dummy_input, forget_input, real_survey          # noqa: F401 (the fixture)
from .utils import needs_backend

INPUT_MODULE = sys.modules[Input.__module__]


@pytest.fixture
def run(tmp_path, monkeypatch):
    """Runs the function with the input files in a temporary directory (the function has no argument for it) and cleans up."""
    (tmp_path / 'in').mkdir()
    (tmp_path / 'out').mkdir()
    monkeypatch.setattr(INPUT_MODULE, 'GALAXIA_TMP', tmp_path / 'in')
    outputs = []

    def run(*args, **kwargs):
        outputs.append(make_survey_from_particles(*args, **{'output_dir': tmp_path / 'out', 'verbose': False, **LIMITS, **kwargs}))
        return outputs[-1]
    yield run
    for output in outputs:
        forget_input(output.survey.input)


@pytest.mark.integration
@needs_backend
def test_the_function_runs_a_survey_and_gives_its_output(run):
    particles, kernels = dummy_input()
    output = run(particles, kernels, simname='itest_public')
    assert len(output._vaex) > 100 and output._hdf5s[0].is_file()
    assert output.survey.input.name == 'itest_public' and output.survey.surveyname == 'survey'


@pytest.mark.integration
@needs_backend
def test_the_same_input_and_parameters_make_the_same_survey(run, real_survey):
    particles, kernels = dummy_input()
    again = run(particles, kernels, simname='itest_again', surveyname='itest_again')
    first, second = real_survey.output._vaex, again._vaex
    assert len(second) == len(first)
    for key in ('px', 'py', 'pz', 'minit', 'gaia__dr2_g'):
        assert np.sort(second[key].to_numpy()).tolist() == np.sort(first[key].to_numpy()).tolist(), key


@pytest.mark.integration
@needs_backend
def test_the_files_of_an_input_make_the_same_survey_as_the_input(run, real_survey):
    input = real_survey.survey.input
    again = run(pname=input.pname, kname=input.kname)
    first, second = real_survey.output._vaex, again._vaex
    assert len(second) == len(first)
    assert np.sort(second['px'].to_numpy()).tolist() == np.sort(first['px'].to_numpy()).tolist()


@pytest.mark.xfail(raises=UnboundLocalError, reason="a call with the particles and no kernels is neither of the two documented calls, and "
                   "the function goes on with no input, which fails on a variable that was not made, where the message should say "
                   "what the calls are")
def test_a_call_that_is_none_of_the_documented_ones_says_what_the_calls_are(tmp_path):
    particles, _ = dummy_input()
    with pytest.raises((TypeError, ValueError), match="particles"):
        make_survey_from_particles(particles, output_dir=tmp_path)


if __name__ == '__main__':
    pass
