#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
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
Package parameters
"""
import os
import sys
import pathlib

from ._name import *

__all__ = ['NAME', 'GALAXIA_SUBMODULE_NAME', 'GALAXIA_EXEC', 'SRC_DIR', 'LOG_DIR', 'NBODY1', 'LEGACY_PHOTOCAT', 'CUSTOM_PHOTOCAT', 'HASH_EXT', 'HASH_ENCODING', 'GALAXIA_DATA', 'GALAXIA_NBODY1', 'GALAXIA_FILENAMES', 'GALAXIA_LOG', 'GALAXIA', 'ISOCHRONES_PATH', 'CACHE', 'STARCATALOG_GROUP', 'PARENTPARTICLE_GROUP', 'PARTITIONING_GROUP', 'METADATA_HEADER_KEY', 'PARAMETER_METADATA_PREFIX', 'INPUT_METADATA_PREFIX', 'SURVEY_METADATA_PREFIX']

GALAXIA_SUBMODULE_NAME = 'Galaxia-ananke'
GALAXIA_EXEC = 'galaxia'
SRC_DIR = 'src'
BIN_DIR = 'bin'
LOG_DIR = 'log'
TMP_DIR = 'tmp'
GALDATA_DIR = 'GalaxiaData'

NBODY1 = 'nbody1'
FILENAMES = 'filenames'
LEGACY_PHOTOCAT = 'padova'
CUSTOM_PHOTOCAT = 'py_custom'
HASH_EXT = 'hash'
HASH_ENCODING = 'ascii'

PREFIX_ENV_VAR = "ANANKE_SYSTEM_PREFIX"
if PREFIX_ENV_VAR in os.environ:
    PREFIX = pathlib.Path(os.environ[PREFIX_ENV_VAR])
else:
    PREFIX = pathlib.Path(sys.prefix)

GLOBAL_CACHE = PREFIX / '.cache'

CACHE = pathlib.Path(GLOBAL_CACHE).expanduser().resolve() / NAME
GALAXIA_DATA = CACHE / GALDATA_DIR
GALAXIA_NBODY1 = GALAXIA_DATA / NBODY1
GALAXIA_FILENAMES = GALAXIA_NBODY1 / FILENAMES
GALAXIA_LOG = CACHE / LOG_DIR
GALAXIA = CACHE / BIN_DIR / GALAXIA_EXEC

ISOCHRONES_PATH = GALAXIA_DATA / 'Isochrones'

STARCATALOG_GROUP = 'StarCatalog'
PARENTPARTICLE_GROUP = 'ParentParticles'
PARTITIONING_GROUP = 'Partitioning'
METADATA_HEADER_KEY = 'Header'
PARAMETER_METADATA_PREFIX = 'SurveyParam_'
INPUT_METADATA_PREFIX = 'Input_'
SURVEY_METADATA_PREFIX = 'Survey_'


def backend_missing_message(url):
    """What the AssertionError of the import says when the executable of the backend is not where the installation puts it."""
    distribution = NAME.replace('_', '-')
    return (f"Galaxia's executable {GALAXIA} doesn't exist, which means that {NAME}'s installation didn't build the backend "
            f"{GALAXIA_SUBMODULE_NAME} submodule appropriately.\n"
            "The backend is built when the installer builds the wheel, which is not always for the environment that it "
            f"installs into: uv builds in a temporary environment (give {PREFIX_ENV_VAR} the path of the environment to "
            "install into), and uv and pip reuse a wheel that they built before, without building it again "
            f"(`uv cache clean {distribution}` or `pip install --no-cache-dir` make them build it).\n"
            f"Please consult {url}/#troubleshooting-installation for troubleshooting.")


if __name__ == '__main__':
    raise NotImplementedError()
