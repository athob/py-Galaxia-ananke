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
Tests of what importing the package does that nothing else checks: it patches vaex (how it reads HDF5 files and opens many
files, and what it logs), and it warns, on purpose, that it does. The warning is tested in a Python of its own, because
the package is imported once per process and the warning only comes then.
"""
import logging
import pathlib
import subprocess
import sys

import vaex

from ..galaxia_ananke import _patch_vaex

ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_importing_the_package_warns_about_the_patch_it_makes_to_vaex():
    done = subprocess.run([sys.executable, "-W", "always", "-c", "import src.galaxia_ananke"],
                          cwd=ROOT, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert "RuntimeWarning: Importing galaxia_ananke monkey-patches vaex's HDF5 reading" in done.stderr


def test_importing_the_package_replaces_how_vaex_reads_and_opens_files():
    assert vaex.hdf5.dataset.Hdf5MemoryMapped._load is _patch_vaex.new_hdf5_load
    assert vaex.open_many is _patch_vaex.new_open_many


def test_importing_the_package_silences_the_logging_of_vaex():
    assert logging.getLogger("vaex").level == logging.CRITICAL


if __name__ == '__main__':
    pass
