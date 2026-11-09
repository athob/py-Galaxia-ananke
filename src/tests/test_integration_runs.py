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
from .real_backend import real_generations                              # noqa: F401 (the fixtures)
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


if __name__ == '__main__':
    pass
