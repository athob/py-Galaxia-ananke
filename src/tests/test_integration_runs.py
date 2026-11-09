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
Integration tests of the variants of a survey with the real backend: several generations, a second photometric system, partitioning
rules, caching, the photometric systems of the installation. Each variant is run once for all the tests (see real_backend.py), with
the same small input as the survey of test_integration.py.
"""
import json

import ebf
import h5py
import numpy as np
import pytest

from ..galaxia_ananke._constants import STARCATALOG_GROUP
from .real_backend import real_generations, real_systems                # noqa: F401 (the fixtures)
from .synthetic import dtype_of
from .utils import needs_backend

pytestmark = [pytest.mark.integration, needs_backend]


def in_files(output, key):
    """The column of a key in each file of the backend."""
    return [ebf.read(str(path), f"/{key}") for path in output._ebfs]


def in_partition(path, key):
    with h5py.File(path, 'r') as f5:
        return f5[STARCATALOG_GROUP][key][...]


# -- several generations ------------------------------------------------------------

def test_each_generation_has_its_file_of_the_backend(real_generations):
    output = real_generations.output
    assert real_generations.survey.n_gens == [0, 1]
    assert [path.name.split('.')[-2] for path in output._ebfs] == ['0', '1']
    assert all(len(px) > 100 for px in in_files(output, 'px'))


def test_the_generations_are_different_samples_of_the_same_input(real_generations):
    px, parentid = in_files(real_generations.output, 'px'), in_files(real_generations.output, 'parentid')
    assert not np.array_equal(np.sort(px[0]), np.sort(px[1]))
    assert set(parentid[0].tolist()) & set(parentid[1].tolist())            # the stars of both come from the same particles
    assert all(ids.min() >= 0 and ids.max() < len(real_generations.particles['massinit']) for ids in parentid)


def test_the_stars_of_all_the_generations_are_in_the_file_of_the_partition_in_the_order_of_the_generations(real_generations):
    output = real_generations.output
    px = in_files(output, 'px')
    assert list(output._hdf5s) == [0]
    assert in_partition(output._hdf5s[0], 'px').tolist() == np.concatenate(px).tolist()
    assert len(output._vaex) == len(px[0]) + len(px[1])


def test_the_logs_of_all_the_generations_are_kept_and_say_what_varies(real_generations):
    with h5py.File(real_generations.output._hdf5s[0], 'r') as f5:
        logs = json.loads(f5.attrs['Galaxia_ebf_output_logs'])
    assert logs['metadata_name'] == 'Galaxia_ebf_output_logs' and logs['metadata_count'] == 2
    assert {'outputFile', 'ngen'} <= set(logs['varying_keys'])
    assert [variation['ngen'] for variation in logs['variations']] == ['0', '1']


# -- a second photometric system ----------------------------------------------------

GAIA = ['gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp']
LSST = ['lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y']


def test_the_keys_of_the_second_system_are_those_of_the_catalogue_with_the_keys_of_the_first(real_systems):
    output = real_systems.output
    assert len(output.export_keys) == 34 + 6 and len(output.catalogue_keys) == 44 + 6
    assert set(GAIA + LSST) <= set(output.export_keys)
    assert set(output._vaex.get_column_names()) == set(output.catalogue_keys)


def test_the_magnitudes_of_both_systems_are_single_precision_and_all_finite(real_systems):
    for key in GAIA + LSST:
        column = real_systems.output._vaex[key].to_numpy()
        assert column.dtype == np.float32 and np.isfinite(column).all(), key


def test_the_types_of_the_keys_of_two_systems_are_those_of_the_synthetic_output(real_systems):
    output, path = real_systems.output, str(real_systems.output._ebfs[0])
    magnitudes = {key for photosystem in output.photosystems for key in photosystem.to_export_keys}
    assert magnitudes == set(GAIA + LSST)
    for key in output.export_keys:
        assert ebf.getHeader(path, f"/{key}").get_dtype() == dtype_of(key, magnitudes), key


def test_the_attributes_of_the_file_say_which_systems_made_the_survey(real_systems):
    with h5py.File(real_systems.output._hdf5s[0], 'r') as f5:
        attrs = dict(f5.attrs)
    assert (attrs['SurveyParam_photo_categ'], attrs['SurveyParam_photo_sys']) == ('padova', 'GAIA__DR2')
    assert (attrs['SurveyParam_photo_categ_1'], attrs['SurveyParam_photo_sys_1']) == ('padova', 'LSST')


def test_the_second_system_is_appended_to_each_file_of_the_backend_and_the_catalogue_has_all_its_stars(real_systems):
    output = real_systems.output
    assert len(output._ebfs) == 1
    assert all(len(ebf.read(str(output._ebfs[0]), f"/{key}")) == len(output._vaex) for key in LSST)


if __name__ == '__main__':
    pass
