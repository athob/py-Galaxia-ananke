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
Tests of the steps that build and install the backend (configure, make, make install, make distclean), and of the whole
build, with a toy backend in place of Galaxia: a directory with a configure script that writes a Makefile, a Makefile that
makes an executable and installs it in the prefix, and some data to copy. They run the real make.
"""
import pathlib
import shutil
import subprocess
import textwrap
from unittest import mock

import pytest

pytest.importorskip('setuptools')
pytest.importorskip('distutils.errors')

from distutils.errors import CompileError

from .. import _build_utils as build
from .utils import in_tmp_wd, patch_cache

pytestmark = pytest.mark.skipif(shutil.which('make') is None, reason='needs make')

CONFIGURE = """#!/bin/sh
echo "$@" > configure.args
echo configuring
echo "a warning" >&2
{configure_fails}
for argument in "$@"; do
  case "$argument" in --prefix=*) prefix="${{argument#--prefix=}}" ;; esac
done
{{ echo "prefix = $prefix"; cat Makefile.template; }} > Makefile
"""

MAKEFILE = """all:
\t@echo building
\t@{make_fails}
\t@printf '#!/bin/sh\\necho toy galaxia\\n' > galaxia
\t@chmod +x galaxia

install: all
\t@{install_fails}
\t@mkdir -p $(prefix)/bin
\t@cp galaxia $(prefix)/bin/galaxia

distclean:
\t@{distclean_fails}
\t@rm -f Makefile galaxia configure.args
"""


def toy_galaxia(base, fails=None):
    """The toy backend in base/Galaxia-ananke, in which the step called `fails` (if any) exits with an error."""
    toy = pathlib.Path(base) / 'Galaxia-ananke'
    toy.mkdir(parents=True)
    (toy / 'configure').write_text(CONFIGURE.format(configure_fails='exit 1' if fails == 'configure' else ''))
    (toy / 'configure').chmod(0o755)
    (toy / 'Makefile.template').write_text(MAKEFILE.format(**{f'{step}_fails': 'false' if fails == step else 'true'
                                                              for step in ('make', 'install', 'distclean')}))
    data = toy / 'GalaxiaData' / 'Isochrones' / 'padova' / 'GAIA__DR2'
    data.mkdir(parents=True)
    (data / 'iso.dat').write_text('isochrones')
    return toy


@pytest.fixture
def cache(tmp_path, monkeypatch):
    places = patch_cache(monkeypatch, tmp_path / 'prefix')
    places.GALAXIA_LOG.mkdir(parents=True)               # build_and_install_galaxia makes it; the steps alone expect it
    return places


@pytest.fixture
def toy(tmp_path):
    return toy_galaxia(tmp_path / 'source')


@pytest.fixture(autouse=True)
def quiet():
    with mock.patch.object(build, 'say', lambda text: None):
        yield


def run_toy(path):
    return subprocess.run([str(path)], capture_output=True, text=True).stdout


# -- the steps ----------------------------------------------------------------------

def test_configure_is_given_the_prefix_and_the_directory_of_the_data(cache, toy):
    build.configure_galaxia(toy)
    assert (toy / 'configure.args').read_text().split() == [f'--prefix={cache.CACHE}', f'--datadir={cache.GALAXIA_DATA}']
    assert (toy / 'Makefile').is_file()


def test_the_log_of_a_step_keeps_what_it_says_on_both_of_its_streams(cache, toy):
    build.configure_galaxia(toy)
    assert sorted((cache.GALAXIA_LOG / 'Galaxia-configure.log').read_text().split('\n')) == ['', 'a warning', 'configuring']


def test_make_makes_the_executable_and_its_log_says_what_it_did(cache, toy):
    build.configure_galaxia(toy)
    build.make_galaxia(toy)
    assert run_toy(toy / 'galaxia') == 'toy galaxia\n'
    assert 'building' in (cache.GALAXIA_LOG / 'Galaxia-make.log').read_text()


def test_make_install_installs_the_executable_in_the_prefix_and_copies_the_data(cache, toy):
    build.configure_galaxia(toy)
    build.make_galaxia(toy)
    build.make_install_galaxia(toy)
    assert run_toy(cache.GALAXIA) == 'toy galaxia\n'
    assert (cache.GALAXIA_DATA / 'Isochrones' / 'padova' / 'GAIA__DR2' / 'iso.dat').read_text() == 'isochrones'
    assert (cache.GALAXIA_LOG / 'Galaxia-make-install.log').is_file()


def test_make_install_does_not_overwrite_data_that_is_already_there(cache, toy):
    # the build removes the earlier installation first: data that is still there means that something went wrong
    cache.GALAXIA_DATA.mkdir(parents=True)
    build.configure_galaxia(toy)
    build.make_galaxia(toy)
    with pytest.raises(FileExistsError):
        build.make_install_galaxia(toy)


def test_make_distclean_cleans_what_configure_and_make_made(cache, toy):
    build.configure_galaxia(toy)
    build.make_galaxia(toy)
    build.make_distclean_galaxia(toy)
    assert not (toy / 'Makefile').exists() and not (toy / 'galaxia').exists()
    assert (cache.GALAXIA_LOG / 'Galaxia-make-distclean.log').is_file()


def test_each_step_has_a_log_of_its_own(cache, toy):
    build.configure_galaxia(toy)
    build.make_galaxia(toy)
    build.make_install_galaxia(toy)
    build.make_distclean_galaxia(toy)
    assert sorted(path.name for path in cache.GALAXIA_LOG.iterdir()) == [
        'Galaxia-configure.log', 'Galaxia-make-distclean.log', 'Galaxia-make-install.log', 'Galaxia-make.log']


# -- the whole build ----------------------------------------------------------------

def test_the_build_installs_the_backend_with_its_data_and_leaves_the_source_clean(cache, toy):
    build.build_and_install_galaxia(toy)
    assert run_toy(cache.GALAXIA) == 'toy galaxia\n'
    assert (cache.GALAXIA_DATA / 'Isochrones' / 'padova' / 'GAIA__DR2' / 'iso.dat').is_file()
    assert sorted(path.name for path in cache.GALAXIA_LOG.iterdir()) == [
        'Galaxia-configure.log', 'Galaxia-make-distclean.log', 'Galaxia-make-install.log', 'Galaxia-make.log']
    assert not (toy / 'Makefile').exists() and not (toy / 'configure.args').exists()      # distclean ran


def test_the_build_says_what_it_is_doing(cache, toy):
    said = []
    with mock.patch.object(build, 'say', said.append):
        build.build_and_install_galaxia(toy)
    assert said == ['\nBuilding Galaxia', '\n\tConfiguring', '\n\tRunning make', '\n\tRunning make install',
                    '\n\tRunning make distclean', '\n\tCleaning temporary', '\n']


def test_the_build_makes_the_directories_it_needs(tmp_path, monkeypatch, toy):
    places = patch_cache(monkeypatch, tmp_path / 'not-made-yet')
    assert not places.CACHE.exists()
    build.build_and_install_galaxia(toy)
    assert places.GALAXIA.is_file() and places.GALAXIA_LOG.is_dir()


def test_the_steps_are_taken_in_order(cache, toy):
    taken = []
    for name in ('remove_existing_galaxia', 'touch_autotools_files', 'configure_galaxia', 'make_galaxia',
                 'make_install_galaxia', 'make_distclean_galaxia', 'clean_up_temporary'):
        original = getattr(build, name)
        mock.patch.object(build, name, lambda *args, _name=name, _original=original: (taken.append(_name), _original(*args))[1]).start()
    try:
        build.build_and_install_galaxia(toy)
    finally:
        mock.patch.stopall()
    assert taken == ['remove_existing_galaxia', 'touch_autotools_files', 'configure_galaxia', 'make_galaxia',
                     'make_install_galaxia', 'make_distclean_galaxia', 'clean_up_temporary']


def test_the_earlier_installation_is_replaced(cache, toy):
    (cache.CACHE / 'bin').mkdir(parents=True)
    (cache.CACHE / 'bin' / 'galaxia').write_text('old')
    (cache.CACHE / 'stale.txt').write_text('left by the earlier installation')
    build.build_and_install_galaxia(toy)
    assert not (cache.CACHE / 'stale.txt').exists()
    assert run_toy(cache.GALAXIA) == 'toy galaxia\n'


def test_the_custom_photometric_systems_survive_a_new_installation(cache, toy):
    mine = cache.ISOCHRONES_PATH / build.CUSTOM_PHOTOCAT / 'Mine'
    mine.mkdir(parents=True)
    (mine / 'isochrone.dat').write_text('my own system')
    build.build_and_install_galaxia(toy)
    assert (cache.ISOCHRONES_PATH / build.CUSTOM_PHOTOCAT / 'Mine' / 'isochrone.dat').read_text() == 'my own system'
    assert (cache.ISOCHRONES_PATH / 'padova' / 'GAIA__DR2' / 'iso.dat').is_file()         # next to the new data


def test_without_custom_photometric_systems_the_new_installation_has_none(cache, toy):
    (cache.CACHE / 'bin').mkdir(parents=True)
    build.build_and_install_galaxia(toy)
    assert not (cache.ISOCHRONES_PATH / build.CUSTOM_PHOTOCAT).exists()


def test_the_backend_to_build_can_be_given_as_a_relative_name_or_a_string(cache, tmp_path, monkeypatch):
    toy_galaxia(tmp_path / 'there')
    monkeypatch.chdir(tmp_path / 'there')
    build.build_and_install_galaxia('Galaxia-ananke')               # as build_ext gives it
    assert run_toy(cache.GALAXIA) == 'toy galaxia\n'


@pytest.mark.parametrize('step', ['configure', 'make', 'install', 'distclean'])
@pytest.mark.xfail(raises=pytest.fail.Exception, reason="the exit status of the steps of the build is not looked at: a build "
                   "that fails goes on, and the package installs without a backend, which is only found out at the first import")
def test_a_step_that_fails_is_an_error(step, tmp_path, monkeypatch):
    places = patch_cache(monkeypatch, tmp_path / 'prefix')
    with pytest.raises((OSError, CompileError, subprocess.CalledProcessError)):
        build.build_and_install_galaxia(toy_galaxia(tmp_path / 'source', fails=step))


if __name__ == '__main__':
    pass
