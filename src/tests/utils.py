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
import os
import sys
import io
import tempfile
import random
import string
import pathlib
import stat
import types
import itertools
import contextlib
import functools

import numpy as np
import pytest


FOURTHIRDPI = 4*np.pi/3

def make_random_string(length=10):
    return "".join([
        random.choice(string.ascii_letters+string.digits)
        for i in range(length)
        ])

@contextlib.contextmanager
def tmp_wd():
    """
    Changes working directory to a temporary one
    and returns to previous on exit.
    Credit to https://stackoverflow.com/a/42441759
    """
    prev_cwd = pathlib.Path.cwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        try:
            yield
        finally:
            os.chdir(prev_cwd)

def in_tmp_wd(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        with tmp_wd():
            func(*args, **kwargs)
    return wrapped_func

@contextlib.contextmanager
def fake_commands(commands):
    """
    Makes the commands of the dictionary, and only them, available on the PATH until the end of the block:
    each one prints its text (and exits with its status, if the value is a tuple of the text and the status),
    whatever it is given. Needs a Unix system, like the package.
    """
    with tempfile.TemporaryDirectory() as directory:
        for name, output in commands.items():
            text, status = output if isinstance(output, tuple) else (output, 0)
            quoted = text.replace("'", "'\\''")                          # printf is a builtin: no PATH needed
            script = pathlib.Path(directory, name)
            script.write_text(f"#!/bin/sh\nprintf '%s\\n' '{quoted}'\nexit {status}\n")
            script.chmod(script.stat().st_mode | stat.S_IEXEC)
        previous = os.environ['PATH']
        os.environ['PATH'] = directory
        try:
            yield
        finally:
            os.environ['PATH'] = previous

PACKAGE_DIR = pathlib.Path(__file__).resolve().parents[1] / "galaxia_ananke"
_fresh_counter = itertools.count()

@contextlib.contextmanager
def fresh_package():
    """
    A copy of the package under a name of its own until the end of the block, made of its files and without running
    its __init__ (which imports everything): the block imports the modules that it needs from it,
    `importlib.import_module(f"{package}._constants")` for one, and gets them as they are when they are first
    imported, which is when what they read (the environment, the version) is what the test set up.
    """
    name = f"galaxia_ananke_under_test_{next(_fresh_counter)}"
    package = types.ModuleType(name)
    package.__path__ = [str(PACKAGE_DIR)]
    sys.modules[name] = package
    try:
        yield name
    finally:
        for key in [key for key in sys.modules if key == name or key.startswith(name + ".")]:
            del sys.modules[key]

def patch_cache(monkeypatch, tmp_path):
    """
    Puts the places of the cache of the backend (where a build installs it, with its data and its logs) under a
    temporary directory, in _build_utils, which took them from the constants when it was imported.
    """
    from .. import _build_utils as build
    root = tmp_path / '.cache' / 'galaxia_ananke'
    places = {'CACHE': root, 'GALAXIA': root / 'bin' / 'galaxia', 'GALAXIA_DATA': root / 'GalaxiaData',
              'GALAXIA_LOG': root / 'log', 'ISOCHRONES_PATH': root / 'GalaxiaData' / 'Isochrones'}
    for name, path in places.items():
        monkeypatch.setattr(build, name, path)
    return types.SimpleNamespace(**places)

def minimal_particles():
    """The smallest valid particles dictionary: two particles, the five required keys, values that a float holds exactly."""
    return {'age': np.array([9.5, 10.0]),
            'feh': np.array([-0.5, -1.0]),
            'massinit': np.array([1000.0, 2000.0]),
            'pos3': np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
            'vel3': np.array([[0.5, 1.5, 2.5], [3.5, 4.5, 5.5]])}


def minimal_kernels():
    """Phase space kernels (Nx2) for minimal_particles."""
    return np.array([[1.0, 2.0], [3.0, 4.0]])

def patch_nbody1(monkeypatch, tmp_path):
    """
    Puts the directory where the backend finds its inputs under a temporary directory, in Input, which took it from the
    constants when it was imported: so that preparing an input never writes in the cache of the installation.
    """
    from ..galaxia_ananke import Input
    module = sys.modules[Input.__module__]
    directory = tmp_path / 'GalaxiaData' / 'nbody1'
    (directory / 'filenames').mkdir(parents=True)
    monkeypatch.setattr(module, 'GALAXIA_NBODY1', directory)
    monkeypatch.setattr(module, 'GALAXIA_FILENAMES', directory / 'filenames')
    return directory

def write_isochrone_file(path, header, rows):
    """An isochrone file as the padova ones are: a commented header with the names of the columns, then a row for each star."""
    path = pathlib.Path(path)
    path.write_text("# " + " ".join(header) + "\n" + "".join(" ".join(str(value) for value in row) + "\n" for row in rows))
    return path

SVO_DATA = pathlib.Path(__file__).resolve().parent / "data"


def custom_isochrone_data():
    """Data for an isochrone made in Python: a table for each metallicity that is needed, with the magnitudes before the rest."""
    import pandas as pd
    from ..galaxia_ananke.photometry.Isochrone import Isochrone
    table = pd.DataFrame({'G': [5.0, 3.0, 5.1], 'BP': [5.5, 3.2, 5.6], 'Age': [9.0, 9.0, 10.0], 'M_ini': [1.0, 2.0, 1.0],
                          'M_act': [0.9, 1.8, 0.95], 'Lum': [0.0, 1.0, 0.1], 'T_eff': [3.7, 3.8, 3.7], 'Grav': [4.4, 4.0, 4.3]})
    return {feh: table.copy() for feh in sorted(Isochrone._required_metallicities)}


def replay_svo(monkeypatch, **files):
    """
    Makes the SVO filter profile service answer from tables that were saved, so that nothing here uses the network: each keyword
    is a facility and gives the file (in tests/data) of the table that astroquery returned for it, or a table. Returns the list of
    the facilities that were asked about, in order. A table is made again at each question, as the code that asks changes what it gets.
    """
    from astropy.table import Table
    from ..galaxia_ananke.photometry.SvoFpsDriver import SvoFpsDriver
    module = sys.modules[SvoFpsDriver.__module__]
    asked = []

    def get_filter_list(facility=None, **kwargs):
        asked.append(facility)
        source = files[facility]
        return source.copy() if isinstance(source, Table) else Table.read(SVO_DATA / source)
    monkeypatch.setattr(module.SvoFps, 'get_filter_list', get_filter_list)
    return asked

def write_catalogue_hdf5(path, columns, group='StarCatalog', attrs=None):
    """An HDF5 file with a group of one dataset for each column, as the files of a catalogue have them."""
    import h5py
    with h5py.File(path, 'w') as f5:
        for name, values in columns.items():
            f5.require_group(group).create_dataset(name=name, data=values)
        for key, value in (attrs or {}).items():
            f5.attrs[key] = value
    return pathlib.Path(path)

def backend_is_built():
    from ..galaxia_ananke._constants import GALAXIA
    return GALAXIA.exists() and os.access(GALAXIA, os.X_OK)


needs_backend = pytest.mark.skipif(not backend_is_built(), reason="the backend is not built")
needs_network = pytest.mark.skipif(not os.environ.get('GALAXIA_TEST_NETWORK'), reason="the network is only used where GALAXIA_TEST_NETWORK is set")

class list_stdout(list):
    """
    Credit to https://stackoverflow.com/a/16571630
    """
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = io.StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout


if __name__ == '__main__':
    pass
