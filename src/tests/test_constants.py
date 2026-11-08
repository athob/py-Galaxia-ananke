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
Tests of the constants of the package: the names of the groups and of the keys that the files of a survey are written with,
and the places that the cache of the backend is searched at, which follow the prefix of the installation.
"""
import importlib
import os
import sys

import pytest

from ..galaxia_ananke import _constants
from .utils import fresh_package


def constants_with_prefix(monkeypatch, prefix):
    """The constants as a fresh import reads them: the prefix is that of the environment variable, or sys.prefix without it."""
    if prefix is None:
        monkeypatch.delenv("ANANKE_SYSTEM_PREFIX", raising=False)
    else:
        monkeypatch.setenv("ANANKE_SYSTEM_PREFIX", str(prefix))
    return importlib.import_module(f"{PACKAGE}._constants")


@pytest.fixture
def package():
    global PACKAGE
    with fresh_package() as PACKAGE:
        yield PACKAGE


# -- the names that the files of a survey are written with --------------------------

@pytest.mark.parametrize("name, value", [
    ("NAME", "galaxia_ananke"),
    ("GALAXIA_EXEC", "galaxia"),
    ("GALAXIA_SUBMODULE_NAME", "Galaxia-ananke"),    # the installer finds the backend by this name, in .gitmodules
    ("NBODY1", "nbody1"),
    ("LEGACY_PHOTOCAT", "padova"),
    ("CUSTOM_PHOTOCAT", "py_custom"),
    ("HASH_EXT", "hash"),
    ("HASH_ENCODING", "ascii"),
    ("STARCATALOG_GROUP", "StarCatalog"),            # the groups and the keys of the HDF5 files: whoever reads them directly...
    ("PARENTPARTICLE_GROUP", "ParentParticles"),
    ("PARTITIONING_GROUP", "Partitioning"),
    ("METADATA_HEADER_KEY", "Header"),
    ("PARAMETER_METADATA_PREFIX", "SurveyParam_"),
    ("INPUT_METADATA_PREFIX", "Input_"),
    ("SURVEY_METADATA_PREFIX", "Survey_"),            # ...depends on them
])
def test_the_names_of_the_files_of_a_survey_are_what_they_were_released_with(name, value):
    assert getattr(_constants, name) == value


# -- the cache of the backend -------------------------------------------------------

def test_the_cache_is_in_the_prefix_given_by_the_environment(tmp_path, monkeypatch, package):
    constants = constants_with_prefix(monkeypatch, tmp_path)
    assert constants.PREFIX == tmp_path
    assert constants.CACHE == (tmp_path / ".cache" / "galaxia_ananke").resolve()


def test_without_the_environment_variable_the_prefix_is_that_of_python(monkeypatch, package):
    constants = constants_with_prefix(monkeypatch, None)
    assert constants.PREFIX == __import__("pathlib").Path(sys.prefix)


def test_a_relative_prefix_is_relative_to_the_working_directory(tmp_path, monkeypatch, package):
    monkeypatch.chdir(tmp_path)
    constants = constants_with_prefix(monkeypatch, "relative/prefix")
    assert constants.PREFIX == __import__("pathlib").Path("relative/prefix")
    assert constants.CACHE == (tmp_path / "relative" / "prefix" / ".cache" / "galaxia_ananke").resolve()


def test_everything_the_backend_needs_is_in_the_cache(tmp_path, monkeypatch, package):
    constants = constants_with_prefix(monkeypatch, tmp_path)
    cache = (tmp_path / ".cache" / "galaxia_ananke").resolve()
    assert constants.GALAXIA == cache / "bin" / "galaxia"
    assert constants.GALAXIA_DATA == cache / "GalaxiaData"
    assert constants.ISOCHRONES_PATH == cache / "GalaxiaData" / "Isochrones"
    assert constants.GALAXIA_NBODY1 == cache / "GalaxiaData" / "nbody1"
    assert constants.GALAXIA_FILENAMES == cache / "GalaxiaData" / "nbody1" / "filenames"
    assert constants.GALAXIA_LOG == cache / "log"


def test_a_prefix_that_is_a_symlink_is_resolved_for_the_cache(tmp_path, monkeypatch, package):
    (tmp_path / "real").mkdir()
    (tmp_path / "link").symlink_to(tmp_path / "real")
    constants = constants_with_prefix(monkeypatch, tmp_path / "link")
    assert constants.PREFIX == tmp_path / "link"
    assert constants.CACHE == (tmp_path / "real" / ".cache" / "galaxia_ananke").resolve()


@pytest.mark.xfail(raises=pytest.fail.Exception, reason="the module builds a PermissionError for a prefix that cannot be "
                   "written to and never raises it, so the installation goes on and fails later")
def test_a_prefix_that_cannot_be_written_to_is_refused(tmp_path, monkeypatch, package):
    prefix = tmp_path / "read-only"
    prefix.mkdir()
    prefix.chmod(0o555)
    try:
        with pytest.raises(PermissionError):
            constants_with_prefix(monkeypatch, prefix)
    finally:
        prefix.chmod(0o755)


if __name__ == '__main__':
    pass
