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


def test_a_prefix_that_cannot_be_written_to_does_not_stop_the_import(tmp_path, monkeypatch, package):
    # what the installation needs to write is checked by the build, which is where it is needed: the package can be imported
    # from a prefix that is read-only, to read what a survey made or to list the photometric systems
    prefix = tmp_path / "read-only"
    prefix.mkdir()
    prefix.chmod(0o555)
    try:
        assert constants_with_prefix(monkeypatch, prefix).PREFIX == prefix
    finally:
        prefix.chmod(0o755)


# -- what the import says when the backend is missing ---------------------------------

URL = 'https://github.com/athob/py-Galaxia-ananke'


def test_the_message_of_a_missing_backend_gives_where_it_should_be_which_follows_the_prefix(tmp_path, monkeypatch, package):
    constants = constants_with_prefix(monkeypatch, tmp_path)
    first = constants.backend_missing_message(URL).split('\n')[0]
    executable = (tmp_path / '.cache' / 'galaxia_ananke' / 'bin' / 'galaxia').resolve()
    assert first == (f"Galaxia's executable {executable} doesn't exist, which means that galaxia_ananke's installation "
                     "didn't build the backend Galaxia-ananke submodule appropriately.")


def test_the_message_of_a_missing_backend_says_what_uv_and_pip_do_with_it_and_what_to_do(tmp_path, monkeypatch, package):
    constants = constants_with_prefix(monkeypatch, tmp_path)
    lines = constants.backend_missing_message(URL).split('\n')
    assert len(lines) == 3
    assert ("uv builds in a temporary environment (give ANANKE_SYSTEM_PREFIX the path of the environment to install into)"
            in lines[1])
    assert "uv and pip reuse a wheel that they built before, without building it again" in lines[1]
    assert "`uv cache clean galaxia-ananke` or `pip install --no-cache-dir` make them build it" in lines[1]
    assert lines[2] == f"Please consult {URL}/#troubleshooting-installation for troubleshooting."


def test_the_variable_that_the_message_names_is_the_one_that_gives_the_prefix(tmp_path, monkeypatch, package):
    constants = constants_with_prefix(monkeypatch, tmp_path)
    assert f"give {constants.PREFIX_ENV_VAR} the path" in constants.backend_missing_message(URL)
    assert os.environ[constants.PREFIX_ENV_VAR] == str(tmp_path)


if __name__ == '__main__':
    pass
