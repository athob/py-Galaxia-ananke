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
