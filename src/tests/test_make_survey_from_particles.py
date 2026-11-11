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

import h5py
import numpy as np
import pytest

from ..galaxia_ananke import Input, make_survey_from_particles
from ..galaxia_ananke._constants import PARENTPARTICLE_GROUP, STARCATALOG_GROUP
from .real_backend import LIMITS, dummy_input, forget_input, forget_inputs_named, real_survey          # noqa: F401 (the fixture)
from .utils import needs_backend

INPUT_MODULE = sys.modules[Input.__module__]


@pytest.fixture
def run(tmp_path, monkeypatch):
    """Runs the function with the input files in a temporary directory (the function has no argument for it) and cleans up."""
    (tmp_path / 'in').mkdir()
    (tmp_path / 'out').mkdir()
    monkeypatch.setattr(INPUT_MODULE, 'GALAXIA_TMP', tmp_path / 'in')
    outputs, simnames = [], []

    def run(*args, **kwargs):
        if 'pname' not in kwargs:                    # the name of an input made from files is that of its files
            simname = kwargs.get('simname', 'sim')
            assert simname.startswith('itest'), "a test that runs the backend names its input itest...: it is cleaned up by that name"
            simnames.append(simname)
        outputs.append(make_survey_from_particles(*args, **{'output_dir': tmp_path / 'out', 'verbose': False, **LIMITS, **kwargs}))
        return outputs[-1]
    yield run
    for output in outputs:
        forget_input(output.survey.input)
    for simname in simnames:                        # a run that failed gave no output to ask which input it linked
        forget_inputs_named(simname)


@pytest.mark.integration
@needs_backend
def test_the_function_runs_a_survey_and_gives_its_output(run):
    particles, kernels = dummy_input()
    output = run(particles, kernels, simname='itest_public')
    assert len(output._vaex) > 100 and output._hdf5s[0].is_file()
    assert output.survey.input.name == 'itest_public' and output.survey.surveyname == 'survey'


@pytest.mark.integration
@needs_backend
def test_an_input_with_a_dot_in_its_name_makes_a_survey_too(run):
    particles, kernels = dummy_input()
    output = run(particles, kernels, simname='itest.dotted')
    assert output.survey.input.name == 'itest.dotted' and len(output._vaex) > 100
    assert output._hdf5s[0].name.startswith(f"{output.survey.surveyname_hash}.itest.dotted_")


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


@pytest.mark.integration
@needs_backend
def test_a_survey_without_a_single_star_has_an_empty_catalogue_with_the_columns_and_the_types_of_any_other(run, real_survey):
    particles, kernels = dummy_input()
    with pytest.warns(UserWarning, match="no star in any partition"):
        output = run(particles, kernels, simname='itest_empty', app_mag_lim_hi=-30)       # no star is that bright
    other = real_survey.output
    assert len(output._vaex) == 0 and list(output._hdf5s) == [-1]
    assert set(output._vaex.get_column_names()) == set(other._vaex.get_column_names())

    # The types are those of the files: vaex gives a column that has no row the type that numpy gives an empty array, float64, whatever
    # the type of the column is, so `to_numpy()` of an empty catalogue says nothing of the types of its columns.
    def types(path, group):
        with h5py.File(path, 'r') as f5:
            return {key: str(f5[group][key].dtype) for key in f5[group]}
    for group in (STARCATALOG_GROUP, PARENTPARTICLE_GROUP):
        assert types(output._hdf5s[-1], group) == types(other._hdf5s[0], group), group
    with h5py.File(output._hdf5s[-1], 'r') as empty:
        assert all(len(empty[PARENTPARTICLE_GROUP][key]) == 0 for key in empty[PARENTPARTICLE_GROUP])


@pytest.mark.parametrize("args, keywords", [((0,), {}), ((0, 1, 2), {}), ((), {'pname': 'sim_abcdef0.ebf'}), ((), {})],
                         ids=["the particles only", "three arguments", "pname only", "nothing"])
def test_a_call_that_is_none_of_the_documented_ones_says_what_the_calls_are(args, keywords):
    with pytest.raises(TypeError, match="takes the particles and the kernels as its two arguments, or the keywords pname and kname together"):
        make_survey_from_particles(*args, **keywords)


if __name__ == '__main__':
    pass
