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
Tests of the tools that build the package (src/_build_utils.py): what the system must have, the lookup of the commit of the
backend that GitHub is asked for, what is done with the cache of a previous installation, and the commands that setuptools runs.
Tests that build a backend are in test_build_galaxia.py, and those that fetch one with git in test_check_galaxia_submodule.py.
"""
import contextlib
import errno
import io
import json
import os
import pathlib
import re
import urllib.error
import urllib.request
from unittest import mock

import pytest

pytest.importorskip('setuptools')
pytest.importorskip('distutils.errors')

from .. import _build_utils as build
from .utils import in_tmp_wd, fake_commands, patch_cache

GOOD_COMMANDS = {'git': 'git version 2.43.0', 'make': 'GNU Make 4.3', 'g++': 'g++ (GCC) 11.2.0'}


@contextlib.contextmanager
def quiet():
    """The tools print to the terminal, even when the output is redirected: not in the tests."""
    with mock.patch.object(build, 'say', lambda text: None):
        yield


# -- the system has what is needed --------------------------------------------------

def test_a_system_with_everything_is_accepted():
    with fake_commands(GOOD_COMMANDS):
        build.verify_system_dependencies()


@pytest.mark.parametrize('missing, words', [('git', 'git'), ('make', 'gnumake'), ('g++', 'C++ compiler')])
def test_a_command_that_is_missing_is_refused_by_name(missing, words):
    with fake_commands({name: output for name, output in GOOD_COMMANDS.items() if name != missing}):
        with pytest.raises(OSError, match=f'does not have .*{re.escape(words)}'):
            build.verify_system_dependencies()


def test_the_system_is_checked_for_git_first():
    with fake_commands({}):
        with pytest.raises(OSError, match='does not have git'):
            build.verify_system_dependencies()


@pytest.mark.parametrize('make, accepted', [('4.2', True), ('4.1', False), ('4.2.1', True), ('4.4.1', True), ('3.82', False), ('4.10', True)])
def test_the_oldest_make_is_4_2(make, accepted):
    with fake_commands({**GOOD_COMMANDS, 'make': f'GNU Make {make}'}):
        if accepted:
            build.verify_system_dependencies()
        else:
            with pytest.raises(OSError, match=f'gnumake v{make} installed, but galaxia_ananke requires v4.2'):
                build.verify_system_dependencies()


@pytest.mark.parametrize('gpp, accepted', [('8.5.0', True), ('8.4.9', False), ('13.2.0', True), ('7.5.0', False), ('10.0.0', True)])
def test_the_oldest_gpp_is_8_5(gpp, accepted):
    with fake_commands({**GOOD_COMMANDS, 'g++': f'g++ (GCC) {gpp}'}):
        if accepted:
            build.verify_system_dependencies()
        else:
            with pytest.raises(OSError, match=f'g\\+\\+ v{gpp} installed, but galaxia_ananke requires v8.5'):
                build.verify_system_dependencies()


def test_the_oldest_versions_are_the_ones_the_readme_gives():
    assert (str(build.MIN_MAKE_VERSION), str(build.MIN_GPP_VERSION)) == ('4.2', '8.5')


# -- the helpers --------------------------------------------------------------------

def test_say_writes_the_text_and_a_space(capsys):
    # no terminal here: with one, say would also write to it, in the middle of the output of the tests
    with mock.patch.object(pathlib.Path, 'exists', return_value=False):
        build.say('Building')
    assert capsys.readouterr().out == 'Building '


def test_say_also_writes_to_the_terminal_when_the_output_is_redirected():
    with mock.patch.object(pathlib.Path, 'exists', return_value=True), mock.patch('builtins.open', mock.mock_open()) as opened:
        build.say('Building')
    opened.assert_called_once_with('/dev/tty', 'w')
    opened().write.assert_called_once_with('Building ')


def test_say_does_without_a_terminal():
    with mock.patch.object(pathlib.Path, 'exists', return_value=False), mock.patch('builtins.open', mock.mock_open()) as opened:
        build.say('Building')
    opened.assert_not_called()


@pytest.mark.parametrize('error', [OSError('no tty'), PermissionError('not allowed')])
def test_say_does_without_a_terminal_that_cannot_be_written_to(error):
    with mock.patch.object(pathlib.Path, 'exists', return_value=True), mock.patch('builtins.open', side_effect=error):
        build.say('Building')


def test_say_does_not_write_to_the_terminal_when_the_output_is_one():
    class Terminal(io.StringIO):
        def isatty(self):
            return True
    with mock.patch('sys.stdout', Terminal()), mock.patch('builtins.open', mock.mock_open()) as opened:
        build.say('Building')
    opened.assert_not_called()


def test_all_files_lists_the_files_of_the_directories_and_their_subdirectories_from_the_base(tmp_path):
    for name in ('a/x.txt', 'a/b/y.txt', 'c/z.txt', 'd/ignored.txt'):
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_text('')
    assert sorted(build.all_files('a', 'c', basedir=tmp_path)) == ['a/b/y.txt', 'a/x.txt', 'c/z.txt']


def test_all_files_of_a_directory_that_does_not_exist_is_nothing(tmp_path):
    assert build.all_files('nowhere', basedir=tmp_path) == []


@in_tmp_wd
def test_all_files_takes_the_base_to_be_the_working_directory():
    pathlib.Path('a').mkdir()
    pathlib.Path('a/x.txt').write_text('')
    assert build.all_files('a') == ['a/x.txt']


def test_the_package_data_is_the_license_of_the_package_and_builds_nothing(monkeypatch):
    calls = []
    monkeypatch.chdir(build.ROOT_DIR)
    with mock.patch.object(build, 'verify_system_dependencies', lambda: calls.append('verify')), \
            mock.patch.object(build, 'build_and_install_galaxia', lambda *a: calls.append('build')):
        data = build.make_package_data()
    assert calls == []                                                  # reading setup.py builds nothing
    assert list(data) == [build.NAME]
    assert {'__license__/LICENSE', '__license__/__init__.py'} <= set(data[build.NAME])
    assert all(name.startswith('__license__/') for name in data[build.NAME])


def test_the_autotools_files_are_made_and_made_recent(tmp_path):
    # so that make does not try to remake them with autotools, which the system may not have
    build.touch_autotools_files(tmp_path)
    names = ['Makefile.am', 'configure.ac', 'aclocal.m4', 'Makefile.in', 'configure', 'config.h.in']
    assert all((tmp_path / name).is_file() for name in names)
    old = 1_000_000_000
    for name in names:
        os.utime(tmp_path / name, (old, old))
    build.touch_autotools_files(tmp_path)
    assert all((tmp_path / name).stat().st_mtime > old + 1_000 for name in names)


# -- the commit of the backend that GitHub knows ------------------------------------

class Response:
    """What urlopen gives: a context manager that reads the body and has headers."""
    def __init__(self, body, headers=None):
        self.body, self.headers = json.dumps(body).encode('utf-8'), headers or {}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self):
        return self.body

    def info(self):
        return mock.Mock(get_content_charset=lambda default: default)


def lookup(tree, headers=None, **call):
    """get_submodule_commit with GitHub answering with this tree: (what it returns, the url it asked)."""
    asked = []

    def urlopen(request, context=None):
        asked.append(request.full_url)
        return Response({'tree': tree}, headers)
    arguments = dict(repo_owner='athob', repo_name='py-Galaxia-ananke', commit_sha='abc123', submodule_path='Galaxia-ananke')
    arguments.update(call)
    with mock.patch.object(urllib.request, 'urlopen', urlopen):
        return build.get_submodule_commit(**arguments), asked[0]


TREE = [{'path': 'README.md', 'type': 'blob', 'sha': '1' * 40},
        {'path': 'Galaxia-ananke-notes.md', 'type': 'blob', 'sha': '3' * 40},       # the name of the submodule is the start of its name
        {'path': 'Galaxia-ananke', 'type': 'commit', 'sha': 'c' * 40},
        {'path': 'src', 'type': 'tree', 'sha': '2' * 40}]


def test_the_commit_of_the_submodule_is_the_one_that_github_gives_for_that_commit_of_the_repository():
    sha, url = lookup(TREE)
    assert sha == 'c' * 40
    assert url == 'https://api.github.com/repos/athob/py-Galaxia-ananke/git/trees/abc123?recursive=1'


def test_a_path_that_is_not_a_submodule_is_a_value_error():
    with pytest.raises(ValueError, match="'src' is not a submodule \\(type: tree\\)"):
        lookup(TREE, submodule_path='src')


@pytest.mark.parametrize('tree', [TREE, []])
def test_a_path_that_is_not_in_the_tree_is_a_file_not_found_error(tree):
    with pytest.raises(FileNotFoundError, match="Submodule path 'nowhere' not found"):
        lookup(tree, submodule_path='nowhere')


def test_a_response_that_says_that_the_limit_is_reached_is_an_error_that_says_when_it_resets():
    # even though this response itself is good: the one after it would not be
    with pytest.raises(RuntimeError, match='rate limit exceeded. Reset at 1700000000'):
        lookup(TREE, headers={'X-RateLimit-Remaining': '0', 'X-RateLimit-Reset': '1700000000'})


def test_a_response_with_requests_left_is_used():
    assert lookup(TREE, headers={'X-RateLimit-Remaining': '59'})[0] == 'c' * 40


def http_error(code, reason):
    return urllib.error.HTTPError('https://api.github.com', code, reason, {}, None)


def test_a_refusal_for_the_rate_limit_is_an_error_that_says_to_try_later():
    with mock.patch.object(urllib.request, 'urlopen', side_effect=http_error(403, 'API Rate Limit exceeded')):
        with pytest.raises(RuntimeError, match='Try again later'):
            build.get_submodule_commit('athob', 'py-Galaxia-ananke', 'abc123', 'Galaxia-ananke')


@pytest.mark.parametrize('code, reason', [(403, 'Forbidden'), (404, 'Not Found'), (500, 'Server Error')])
def test_any_other_http_error_is_raised_as_it_is(code, reason):
    with mock.patch.object(urllib.request, 'urlopen', side_effect=http_error(code, reason)):
        with pytest.raises(urllib.error.HTTPError) as error:
            build.get_submodule_commit('athob', 'py-Galaxia-ananke', 'abc123', 'Galaxia-ananke')
    assert error.value.code == code


# -- the cache of a previous installation -------------------------------------------

@pytest.fixture
def cache(tmp_path, monkeypatch):
    """The places of the cache, under a temporary directory: where the build puts what it makes."""
    return patch_cache(monkeypatch, tmp_path)


def make_old_installation(cache, custom=True):
    (cache.CACHE / 'bin').mkdir(parents=True)
    (cache.CACHE / 'bin' / 'galaxia').write_text('old')
    (cache.ISOCHRONES_PATH / 'padova').mkdir(parents=True)
    (cache.ISOCHRONES_PATH / 'padova' / 'iso.dat').write_text('old data')
    if custom:
        (cache.ISOCHRONES_PATH / build.CUSTOM_PHOTOCAT / 'Mine').mkdir(parents=True)
        (cache.ISOCHRONES_PATH / build.CUSTOM_PHOTOCAT / 'Mine' / 'isochrone.dat').write_text('my own system')


def test_without_an_earlier_installation_there_is_nothing_to_remove(cache, tmp_path):
    build.remove_existing_galaxia(tmp_path / 'kept')
    assert not cache.CACHE.exists() and not (tmp_path / 'kept').exists()


def test_the_earlier_installation_is_removed(cache, tmp_path):
    make_old_installation(cache, custom=False)
    build.remove_existing_galaxia(tmp_path / 'kept')
    assert not cache.CACHE.exists() and not (tmp_path / 'kept').exists()


def test_the_custom_photometric_systems_are_put_aside_when_the_earlier_installation_is_removed(cache, tmp_path):
    make_old_installation(cache)
    build.remove_existing_galaxia(tmp_path / 'kept')
    assert not cache.CACHE.exists()
    assert (tmp_path / 'kept' / 'Mine' / 'isochrone.dat').read_text() == 'my own system'


def test_the_custom_photometric_systems_are_put_back_in_the_new_installation(cache, tmp_path):
    make_old_installation(cache)
    build.remove_existing_galaxia(tmp_path / 'kept')
    cache.ISOCHRONES_PATH.mkdir(parents=True)                           # what the new installation made
    build.clean_up_temporary(tmp_path / 'kept')
    assert (cache.ISOCHRONES_PATH / build.CUSTOM_PHOTOCAT / 'Mine' / 'isochrone.dat').read_text() == 'my own system'
    assert not (tmp_path / 'kept').exists()


def test_cleaning_up_when_nothing_was_put_aside_does_nothing(cache, tmp_path):
    cache.ISOCHRONES_PATH.mkdir(parents=True)
    build.clean_up_temporary(tmp_path / 'kept')
    assert list(cache.ISOCHRONES_PATH.iterdir()) == []


@pytest.mark.xfail(raises=OSError, reason="the custom photometric systems are moved with a rename, which fails between two "
                   "file systems (the temporary directory and the cache often are): the installation stops, and has not "
                   "removed the old cache yet, so nothing is lost but it cannot be done")
def test_the_custom_photometric_systems_can_be_put_aside_between_two_file_systems(cache, tmp_path):
    make_old_installation(cache)

    def across_file_systems(*arguments, **keywords):
        raise OSError(errno.EXDEV, 'Invalid cross-device link')
    # pathlib calls os.rename in the recent versions of python, and the function that it kept when it was imported in the older
    # ones (3.9 and 3.10): the two are changed, so that nothing can be renamed across file systems, as shutil.move can cope with
    with mock.patch('os.rename', across_file_systems), mock.patch.object(pathlib.Path, 'rename', across_file_systems):
        build.remove_existing_galaxia(tmp_path / 'kept')
    assert (tmp_path / 'kept' / 'Mine' / 'isochrone.dat').read_text() == 'my own system'


# -- the commands that setuptools runs ----------------------------------------------

def build_ext_of_the_package():
    """The build_ext command of the package, without the layer versioneer puts on it."""
    command = build.make_cmdclass()['build_ext']
    return next(c for c in command.__mro__ if c.__name__ == '_build_ext')


def run_build_ext(**replaced):
    """Run the command with the steps replaced: the order in which they were called, and what they were given."""
    from setuptools import Distribution
    from setuptools.command.build_ext import build_ext
    calls = []
    steps = {'verify_system_dependencies': lambda: calls.append(('verify',)),
             'check_galaxia_submodule': lambda root: calls.append(('check', root)),
             'build_and_install_galaxia': lambda path: calls.append(('build', path))}
    steps.update(replaced)
    with contextlib.ExitStack() as stack:
        stack.enter_context(mock.patch.object(build_ext, 'run', lambda self: calls.append(('extensions',))))
        for name, step in steps.items():
            stack.enter_context(mock.patch.object(build, name, step))
        build_ext_of_the_package()(Distribution()).run()
    return calls


def test_build_ext_goes_through_the_extensions_the_system_the_submodule_and_the_backend_in_that_order():
    assert run_build_ext() == [('extensions',), ('verify',), ('check', build.ROOT_DIR), ('build', build.GALAXIA_SUBMODULE_NAME)]


def test_build_ext_stops_when_the_system_lacks_what_is_needed():
    def refuse():
        raise OSError('no compiler')
    calls = []
    with pytest.raises(OSError, match='no compiler'):
        calls = run_build_ext(verify_system_dependencies=refuse)
    assert calls == []


def test_build_ext_does_not_build_a_backend_that_cannot_be_used():
    built = []

    def refuse(root):
        raise PermissionError('install-sh cannot be run')
    with pytest.raises(PermissionError):
        run_build_ext(check_galaxia_submodule=refuse, build_and_install_galaxia=lambda path: built.append(path))
    assert built == []


@pytest.mark.xfail(raises=AssertionError, reason="the backend is built from its name, which is taken from the working "
                   "directory, while it was checked at the root of the repository: the two are the same only when "
                   "the build runs at the root")
def test_build_ext_builds_the_backend_it_checked_wherever_it_runs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calls = run_build_ext()
    checked = [root for step, *rest in calls if step == 'check' for root in rest]
    built = [path for step, *rest in calls if step == 'build' for path in rest]
    assert pathlib.Path(built[0]).resolve() == (pathlib.Path(checked[0]) / build.GALAXIA_SUBMODULE_NAME).resolve()


def test_the_commands_of_the_build_are_those_of_versioneer_and_a_test_command_that_does_nothing():
    from setuptools import Distribution
    from setuptools.command.build_ext import build_ext
    commands = build.make_cmdclass()
    assert 'sdist' in commands and 'build_py' in commands              # what versioneer adds
    assert issubclass(commands['build_ext'], build_ext)
    test = commands['test'](Distribution())
    assert test.user_options == [] and test.description == 'run tests'
    assert test.initialize_options() is None and test.finalize_options() is None and test.run() is None


if __name__ == '__main__':
    pass
