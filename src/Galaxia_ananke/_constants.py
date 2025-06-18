#!/usr/bin/env python
"""
Package parameters
"""
import sys
import pathlib

from ._name import *

__all__ = ['NAME', 'GALAXIA_SUBMODULE_NAME', 'GALAXIA_URL', 'GALAXIA_EXEC', 'SRC_DIR', 'LOG_DIR', 'NBODY1', 'LEGACY_PHOTOCAT', 'CUSTOM_PHOTOCAT', 'HASH_EXT', 'HASH_ENCODING', 'GALAXIA_DATA', 'GALAXIA_NBODY1', 'GALAXIA_FILENAMES', 'GALAXIA_LOG', 'GALAXIA', 'ISOCHRONES_PATH', 'CACHE']

GALAXIA_SUBMODULE_NAME = 'galaxia-ananke'
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
if PREFIX_ENV_VAR in pathlib.os.environ:
    PREFIX = pathlib.Path(pathlib.os.environ[PREFIX_ENV_VAR])
else:
    PREFIX = pathlib.Path(sys.prefix)
if not pathlib.os.access(PREFIX, pathlib.os.W_OK):
    PermissionError(f"Installation cannot complete: to proceed, please give write permission to directory {PREFIX} or define a custom system prefix via the environment variable {PREFIX_ENV_VAR}")  # TODO make the variable conda-dependent?

GLOBAL_CACHE = PREFIX / '.cache'

CACHE = pathlib.Path(GLOBAL_CACHE).expanduser().resolve() / NAME
GALAXIA_DATA = CACHE / GALDATA_DIR
GALAXIA_NBODY1 = GALAXIA_DATA / NBODY1
GALAXIA_FILENAMES = GALAXIA_NBODY1 / FILENAMES
GALAXIA_LOG = CACHE / LOG_DIR
GALAXIA = CACHE / BIN_DIR / GALAXIA_EXEC

ISOCHRONES_PATH = GALAXIA_DATA / 'Isochrones'


if __name__ == '__main__':
    raise NotImplementedError()
