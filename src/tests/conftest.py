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
What every test of the suite starts from.

The random state of NumPy and of Python is seeded before each test, with a number that the name of the
test determines: a run draws the same values as the one before it, two tests do not draw the same values, and
the seed of a test that fails is in its report. The dummy inputs of the package (`Input.make_dummy_particles_input`
and `make_dummy_kernels_input`) draw from the global state of NumPy, which nothing else seeds.

GALAXIA_TEST_SEED=N shifts all the seeds, to try other draws.

A test works in its own temporary directory, and not in the directory that pytest was run from (where the package makes its
files by default): a test that leaves a file there fails, with the name of the file in its report, and the file is removed.
"""
import os
import pathlib
import random
import zlib

import pytest


def seed_of(nodeid, shift=0):
    """The seed of the test of that name: a 32 bit number that only the name and the shift determine."""
    return (zlib.crc32(nodeid.encode("utf-8")) + shift) % 2**32


def left_in(directory, before):
    """The entries that are in the directory and were not there before (their names), leaving out the hidden ones, which pytest and
    coverage keep there, and the caches of Python."""
    return sorted(path for path in pathlib.Path(directory).iterdir()
                  if path.name not in before and not path.name.startswith('.') and path.name != '__pycache__')


@pytest.fixture(autouse=True)
def nothing_is_left_where_pytest_runs_from(request):
    directory = pathlib.Path(request.config.invocation_params.dir)
    before = {path.name for path in directory.iterdir()}
    yield
    left = left_in(directory, before)
    for path in left:
        if path.is_file() or path.is_symlink():
            path.unlink()                # so that the next test that leaves the same file is reported too, and the directory stays clean
    if left:
        pytest.fail(f"this test left {', '.join(path.name for path in left)} in {directory}, where pytest runs from: "
                    f"a test works in its own temporary directory (tmp_path)", pytrace=False)


@pytest.fixture(autouse=True)
def seeded_random_state(request):
    import numpy as np   # here, so that this file does not need numpy to be read
    shift = int(os.environ.get("GALAXIA_TEST_SEED", "0"))
    seed = seed_of(request.node.nodeid, shift)
    random.seed(seed)
    np.random.seed(seed)
    request.node.add_report_section("setup", "random seed", f"{seed} (the name of the test, shifted by GALAXIA_TEST_SEED={shift})")
