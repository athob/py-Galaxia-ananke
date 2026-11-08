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
Tests of how the backend is fetched (check_galaxia_submodule and clone_and_checkout_submodules), with the real git, on local
repositories laid out as they are on GitHub: <owner>/py-Galaxia-ananke, whose .gitmodules gives the backend by a relative url,
next to <owner>/Galaxia-ananke.git. The .gitmodules is the one of this repository. What GitHub says about the commit that a
version pins is replaced, since nothing here may use the network.
"""
import configparser
import os
import pathlib
import shutil
import subprocess
from unittest import mock

import pytest

pytest.importorskip('setuptools')
pytest.importorskip('distutils.errors')

from distutils.errors import CompileError

from .. import _build_utils as build
from .utils import fake_commands

pytestmark = pytest.mark.skipif(shutil.which('git') is None, reason='needs git')

REVISION = 'f' * 40         # the revision of this repository that the version being installed was made from


def git(*arguments, cwd):
    return subprocess.run(['git', '-c', 'user.email=test@example.org', '-c', 'user.name=Test', '-c', 'commit.gpgsign=false',
                           '-c', 'init.defaultBranch=main', *arguments],
                          cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


@pytest.fixture(autouse=True)
def clean_git(tmp_path, monkeypatch):
    """Git as it is for the tests, whoever runs them: no configuration of the user, nothing above the temporary directory."""
    (tmp_path / 'home').mkdir()
    for name, value in {'HOME': tmp_path / 'home', 'XDG_CONFIG_HOME': tmp_path / 'home', 'GIT_CONFIG_NOSYSTEM': '1',
                        'GIT_CONFIG_GLOBAL': tmp_path / 'home' / 'gitconfig', 'GIT_CEILING_DIRECTORIES': tmp_path,
                        'GIT_TERMINAL_PROMPT': '0'}.items():
        monkeypatch.setenv(name, str(value))
    for name in ('GIT_DIR', 'GIT_WORK_TREE', 'GIT_INDEX_FILE'):
        monkeypatch.delenv(name, raising=False)


@pytest.fixture
def allow_local_repositories(monkeypatch):
    """Git does not take a submodule from a repository on the same machine unless told to."""
    monkeypatch.setenv('GIT_CONFIG_COUNT', '1')
    monkeypatch.setenv('GIT_CONFIG_KEY_0', 'protocol.file.allow')
    monkeypatch.setenv('GIT_CONFIG_VALUE_0', 'always')


@pytest.fixture(autouse=True)
def quiet():
    with mock.patch.object(build, 'say', lambda text: None):
        yield


def declared_url():
    """The relative url of the backend that .gitmodules gives, and who owns it."""
    modules = configparser.ConfigParser()
    modules.read(build.ROOT_DIR / '.gitmodules')
    url = modules[f'submodule "{build.GALAXIA_SUBMODULE_NAME}"']['url']
    return url, url.split('/')[-2]


class World:
    """The backend with two commits, a repository that pins the first, and the layout of GitHub around them."""
    def __init__(self, base):
        url, self.owner = declared_url()
        self.base = base
        self.backend = base / self.owner / url.split('/')[-1]               # <owner>/Galaxia-ananke.git, a bare repository
        self.package_url = base / self.owner / 'py-Galaxia-ananke'          # where the repository would be on GitHub
        work = base / 'backend-work'
        work.mkdir(parents=True)
        git('init', cwd=work)
        (work / 'build-aux').mkdir()
        (work / 'build-aux' / 'install-sh').write_text('#!/bin/sh\n')
        (work / 'build-aux' / 'install-sh').chmod(0o755)
        (work / 'README').write_text('first')
        git('add', '-A', cwd=work)
        git('commit', '-m', 'first', cwd=work)
        self.first = git('rev-parse', 'HEAD', cwd=work)
        (work / 'README').write_text('second')
        git('commit', '-am', 'second', cwd=work)
        self.second = git('rev-parse', 'HEAD', cwd=work)
        self.backend.parent.mkdir(parents=True, exist_ok=True)
        git('init', '--bare', str(self.backend), cwd=base)
        git('push', str(self.backend), 'HEAD:refs/heads/main', cwd=work)

    def sources(self, name='repository', git_repository=True):
        """A copy of the repository that pins the first commit: a clone, or just files (as an sdist unpacks)."""
        root = self.base / name
        root.mkdir()
        shutil.copy(build.ROOT_DIR / '.gitmodules', root / '.gitmodules')
        if git_repository:
            git('init', cwd=root)
            git('remote', 'add', 'origin', str(self.package_url), cwd=root)
            git('add', '.gitmodules', cwd=root)
            git('update-index', '--add', '--cacheinfo', f'160000,{self.first},{build.GALAXIA_SUBMODULE_NAME}', cwd=root)
            git('commit', '-m', 'the repository', cwd=root)
        return root

    def github_says(self, monkeypatch, commit=None):
        """The repository is where GitHub is, and GitHub says that the revision pins `commit` (the first one by default)."""
        asked = []
        monkeypatch.setattr(build, '__url__', 'file://' + str(self.package_url))
        monkeypatch.setattr(build._version, 'get_versions', lambda: {'full-revisionid': REVISION})
        monkeypatch.setattr(build, 'get_submodule_commit', lambda *arguments: asked.append(arguments) or (commit or self.first))
        return asked


@pytest.fixture
def world(tmp_path):
    return World(tmp_path / 'world')


# -- the repository is a git repository ---------------------------------------------

def test_the_url_of_gitmodules_leads_to_the_backend_next_to_the_repository(world, allow_local_repositories):
    root = world.sources()
    build.check_galaxia_submodule(root)
    backend = root / build.GALAXIA_SUBMODULE_NAME
    assert (backend / '.git').exists()
    assert git('rev-parse', 'HEAD', cwd=backend) == world.first                       # the one that the repository pins, not the latest
    assert os.access(backend / 'build-aux' / 'install-sh', os.X_OK)


def test_the_root_of_the_sources_can_be_given_as_a_string(world, allow_local_repositories):
    root = world.sources()
    build.check_galaxia_submodule(str(root))
    assert (root / build.GALAXIA_SUBMODULE_NAME / 'build-aux' / 'install-sh').exists()


def test_a_backend_that_is_already_there_is_left_as_it_is(world, allow_local_repositories):
    root = world.sources()
    build.check_galaxia_submodule(root)
    build.check_galaxia_submodule(root)
    assert git('rev-parse', 'HEAD', cwd=root / build.GALAXIA_SUBMODULE_NAME) == world.first


def test_the_install_script_of_the_backend_has_to_be_executable(world, allow_local_repositories):
    root = world.sources()
    build.check_galaxia_submodule(root)
    script = root / build.GALAXIA_SUBMODULE_NAME / 'build-aux' / 'install-sh'
    script.chmod(0o644)
    with pytest.raises(PermissionError, match=f'user-execute permission to file {script}'):
        build.check_galaxia_submodule(root)


@pytest.mark.xfail(raises=AssertionError, reason="when git cannot initialise the submodule, for any reason but that the "
                   "directory is not a repository, nothing is said about git: the error is about the install script, "
                   "which is not there for want of the backend")
def test_a_submodule_that_git_cannot_initialise_is_an_error_about_git(world):
    # here git refuses the repository on the same machine, as it would one that cannot be reached
    root = world.sources()
    with pytest.raises(OSError, match='git submodule'):
        build.check_galaxia_submodule(root)


def test_git_that_is_not_installed_is_an_error_that_says_so(world):
    root = world.sources()
    with fake_commands({}):
        with pytest.raises(OSError, match='does not have git installed'):
            build.check_galaxia_submodule(root)


# -- the sources are not a git repository, as when an sdist is unpacked --------------

def test_sources_that_are_not_a_git_repository_get_the_backend_by_cloning_it(world):
    root = world.sources('sdist', git_repository=False)
    cloned = []

    def clone(directory, names):
        cloned.append((directory, names))
        script = directory / build.GALAXIA_SUBMODULE_NAME / 'build-aux' / 'install-sh'
        script.parent.mkdir(parents=True)
        script.write_text('#!/bin/sh\n')
        script.chmod(0o755)
    with mock.patch.object(build, 'clone_and_checkout_submodules', clone):
        build.check_galaxia_submodule(root)
    assert cloned == [(root, [build.GALAXIA_SUBMODULE_NAME])]


def test_the_progress_of_the_fetch_is_told_to_the_user(world):
    said = []
    root = world.sources('sdist', git_repository=False)
    with mock.patch.object(build, 'say', said.append), \
            mock.patch.object(build, 'clone_and_checkout_submodules', lambda *a: None), \
            mock.patch.object(os, 'access', return_value=True):
        build.check_galaxia_submodule(root)
    assert said == ['\nChecking submodule Galaxia, running git...',
                    '\n\tFailed to git submodule init Galaxia, attempting to clone...', '\n\tClone was successful']


def test_the_backend_is_cloned_at_the_commit_that_the_revision_pins(world, monkeypatch):
    asked = world.github_says(monkeypatch)
    root = world.sources('sdist', git_repository=False)
    build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])
    backend = root / build.GALAXIA_SUBMODULE_NAME
    assert git('rev-parse', 'HEAD', cwd=backend) == world.first != world.second       # not the latest commit of the backend
    assert asked == [(world.owner, 'py-Galaxia-ananke', REVISION, build.GALAXIA_SUBMODULE_NAME)]


def test_what_git_says_while_it_clones_is_told_to_the_user(world, monkeypatch):
    world.github_says(monkeypatch)
    root = world.sources('sdist', git_repository=False)
    said = []
    with mock.patch.object(build, 'say', said.append):
        build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])
    assert any("Cloning into 'Galaxia-ananke'" in line for line in said)


def test_the_url_that_is_cloned_is_that_of_the_backend_next_to_the_repository_on_github(world, monkeypatch):
    # git is not run here: the command it would be given is what is looked at, with the real address of the package
    clone = []
    monkeypatch.setattr(build._version, 'get_versions', lambda: {'full-revisionid': REVISION})
    monkeypatch.setattr(build, 'get_submodule_commit', lambda *arguments: world.first)
    monkeypatch.setattr(subprocess, 'Popen', lambda command, **kwargs: clone.append((command, kwargs['cwd'])) or mock.Mock(
        stdout=mock.Mock(readline=lambda: ''), wait=lambda: 0))
    monkeypatch.setattr(subprocess, 'call', lambda command, **kwargs: 0)
    root = world.sources('sdist', git_repository=False)
    build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])
    assert clone == [(['git', 'clone', 'https://github.com/athob/Galaxia-ananke.git'], root)]


def test_sources_that_are_not_a_git_repository_end_up_with_the_backend_at_the_pinned_commit(world, monkeypatch):
    world.github_says(monkeypatch)
    root = world.sources('sdist', git_repository=False)
    build.check_galaxia_submodule(root)                                                # the whole path of an installation from an sdist
    backend = root / build.GALAXIA_SUBMODULE_NAME
    assert git('rev-parse', 'HEAD', cwd=backend) == world.first
    assert os.access(backend / 'build-aux' / 'install-sh', os.X_OK)


@pytest.mark.xfail(raises=pytest.fail.Exception, reason="the status of git checkout is not looked at: a commit that cannot be "
                   "checked out leaves the backend at the latest commit of its default branch, and the package is built "
                   "from it without a word")
def test_a_commit_that_cannot_be_checked_out_is_an_error(world, monkeypatch):
    world.github_says(monkeypatch, commit='0' * 40)
    root = world.sources('sdist', git_repository=False)
    with pytest.raises((OSError, CompileError, subprocess.CalledProcessError)):
        build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])


@pytest.mark.xfail(raises=FileNotFoundError, reason="a clone that fails is not noticed (the exception that is caught is the "
                   "wrong one, and the status of git is not looked at): the checkout that follows fails on a directory "
                   "that is not there, with an error that says nothing of the clone")
def test_a_clone_that_fails_is_a_compile_error(world, monkeypatch, tmp_path):
    world.github_says(monkeypatch)
    monkeypatch.setattr(build, '__url__', 'file://' + str(tmp_path / 'nowhere' / 'py-Galaxia-ananke'))
    root = world.sources('sdist', git_repository=False)
    with pytest.raises(CompileError):
        build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])


if __name__ == '__main__':
    pass
