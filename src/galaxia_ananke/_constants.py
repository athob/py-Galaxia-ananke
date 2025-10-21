#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 (AGPL-3.0).
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

__all__ = ['NAME', 'GALAXIA_SUBMODULE_NAME', 'GALAXIA_URL', 'GALAXIA_EXEC', 'SRC_DIR', 'LOG_DIR', 'NBODY1', 'LEGACY_PHOTOCAT', 'CUSTOM_PHOTOCAT', 'HASH_EXT', 'HASH_ENCODING', 'GALAXIA_DATA', 'GALAXIA_NBODY1', 'GALAXIA_FILENAMES', 'GALAXIA_LOG', 'GALAXIA', 'ISOCHRONES_PATH', 'CACHE', 'ATTRIBUTE_HEADER_KEY', 'PARAMETER_ATTRIBUTE_PREFIX']

GALAXIA_SUBMODULE_NAME = 'Galaxia-ananke'
GALAXIA_URL = 'https://sourceforge.net/projects/galaxia/files/galaxia-0.7.2.tar.gz/download'
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
if not os.access(PREFIX, os.W_OK):
    PermissionError(f"Installation cannot complete: to proceed, please give write permission to directory {PREFIX} or define a custom system prefix via the environment variable {PREFIX_ENV_VAR}")  # TODO make the variable conda-dependent?

GLOBAL_CACHE = PREFIX / '.cache'

CACHE = pathlib.Path(GLOBAL_CACHE).expanduser().resolve() / NAME
GALAXIA_DATA = CACHE / GALDATA_DIR
GALAXIA_NBODY1 = GALAXIA_DATA / NBODY1
GALAXIA_FILENAMES = GALAXIA_NBODY1 / FILENAMES
GALAXIA_LOG = CACHE / LOG_DIR
GALAXIA = CACHE / BIN_DIR / GALAXIA_EXEC

ISOCHRONES_PATH = GALAXIA_DATA / 'Isochrones'

ATTRIBUTE_HEADER_KEY = 'Header'
PARAMETER_ATTRIBUTE_PREFIX = 'Survey_'


if __name__ == '__main__':
    raise NotImplementedError()
