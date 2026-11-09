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
Tests of the types of the synthetic output of the backend, against the types of the files of the real backend (they were recorded
from a run of it: test_integration.py compares them again with each run, where the backend is built).
"""
import ebf
import numpy as np
import pytest

from .synthetic import dtype_of, prepared_survey, star_columns, write_backend_output

REAL = {**{key: 'float64' for key in ('age', 'alpha', 'calcium', 'carbon', 'dform', 'dmod', 'feh', 'helium', 'magnesium', 'minit', 'neon',
                                       'nitrogen', 'oxygen', 'px', 'py', 'pz', 'rad', 'silicon', 'sulphur', 'vx', 'vy', 'vz')},
        **{key: 'float32' for key in ('gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp', 'grav', 'lum', 'mact', 'mtip', 'teff')},
        **{key: 'int64' for key in ('parentid', 'partid', 'partitionid', 'satid')}}
MAGNITUDES = ('gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp')


def test_the_types_of_the_synthetic_keys_are_those_that_the_real_backend_gave():
    assert len(REAL) == 34
    assert {key: str(dtype_of(key, MAGNITUDES)) for key in REAL} == REAL


def test_a_magnitude_is_single_precision_only_when_it_is_said_to_be_one():
    assert dtype_of('lsst_g', ('lsst_g',)) == np.float32 and dtype_of('lsst_g') == np.float64


def test_the_columns_that_are_made_have_the_types_of_their_keys():
    columns = star_columns(list(REAL), 4, magnitudes=MAGNITUDES)
    assert {key: str(column.dtype) for key, column in columns.items()} == REAL


def test_the_files_that_are_written_have_the_types_of_the_real_ones(tmp_path, monkeypatch):
    survey = prepared_survey(tmp_path, monkeypatch)
    path, _ = write_backend_output(survey, n_stars=3)[0]
    assert set(survey.output.export_keys) <= set(REAL)
    for key in survey.output.export_keys:
        assert str(ebf.getHeader(str(path), f"/{key}").get_dtype()) == REAL[key], key


if __name__ == '__main__':
    pass
