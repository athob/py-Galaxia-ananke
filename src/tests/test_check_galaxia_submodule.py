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
import gc
import http.client
import json
import os
import pathlib
import shutil
import socket
import ssl
import subprocess
import urllib.error
import warnings
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


def test_a_submodule_that_git_cannot_initialise_is_an_error_about_git(world):
    # here git refuses the repository on the same machine, as it would one that cannot be reached
    root = world.sources()
    with pytest.raises(OSError, match=r'git submodule update --init --recursive failed in .* \(exit status [1-9]'):
        build.check_galaxia_submodule(root)


def test_a_git_repository_that_git_cannot_read_is_an_error_about_git_and_is_not_cloned_over(world):
    root = world.sources('broken', git_repository=False)
    (root / '.git').mkdir()                                  # it is there, and git does not take it for a repository
    with mock.patch.object(build, 'clone_and_checkout_submodules') as clone:
        with pytest.raises(OSError, match='git submodule update --init --recursive failed.*exit status 128'):
            build.check_galaxia_submodule(root)
    clone.assert_not_called()


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


def test_a_commit_that_cannot_be_checked_out_is_an_error(world, monkeypatch):
    world.github_says(monkeypatch, commit='0' * 40)
    root = world.sources('sdist', git_repository=False)
    with pytest.raises(CompileError, match=f"Error checking out the commit {'0' * 40} of {build.GALAXIA_SUBMODULE_NAME} "
                                           r"\(git exit status [1-9]"):
        build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])


def test_a_clone_that_fails_is_a_compile_error_that_gives_what_git_said(world, monkeypatch, tmp_path):
    world.github_says(monkeypatch)
    monkeypatch.setattr(build, '__url__', 'file://' + str(tmp_path / 'nowhere' / 'py-Galaxia-ananke'))
    root = world.sources('sdist', git_repository=False)
    pattern = rf"Error cloning {build.GALAXIA_SUBMODULE_NAME} from file://.*\(git exit status 128\)"
    with pytest.raises(CompileError, match=pattern) as error:
        build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])
    assert 'fatal' in str(error.value) and str(error.value).endswith('aborting...\n')


def test_the_error_of_a_clone_gives_the_last_lines_of_what_git_said_and_no_more(world, monkeypatch):
    world.github_says(monkeypatch)
    root = world.sources('sdist', git_repository=False)
    said = [f'line {i}\n' for i in range(1, 51)]
    process = mock.Mock(stdout=mock.Mock(readline=mock.Mock(side_effect=said + [''])), wait=lambda: 128)
    monkeypatch.setattr(subprocess, 'Popen', lambda command, **kwargs: process)
    with pytest.raises(CompileError) as error:
        build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])
    lines = str(error.value).split('\n')
    assert lines[1:-2] == [f'line {i}' for i in range(31, 51)] and lines[-2:] == ['aborting...', '']


def test_the_pipe_of_the_clone_is_closed_when_it_is_over(world, monkeypatch):
    world.github_says(monkeypatch)
    root = world.sources('sdist', git_repository=False)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])
        gc.collect()
    assert [str(w.message) for w in caught if issubclass(w.category, ResourceWarning)] == []


def test_the_rate_limit_of_github_is_a_compile_error_that_says_what_was_being_cloned(world, monkeypatch):
    world.github_says(monkeypatch)

    def rate_limit(*arguments):
        raise RuntimeError("GitHub API rate limit exceeded. Try again later.")
    monkeypatch.setattr(build, 'get_submodule_commit', rate_limit)
    root = world.sources('sdist', git_repository=False)
    pattern = rf"rate limit exceeded.*\nError finding the commit of {build.GALAXIA_SUBMODULE_NAME} to clone"
    with pytest.raises(CompileError, match=pattern):
        build.clone_and_checkout_submodules(root, [build.GALAXIA_SUBMODULE_NAME])


# -- GitHub cannot be asked for the commit to clone -----------------------------------

def lookup_that_raises(world, monkeypatch, error):
    """Sources that are not a git repository, to which GitHub answers with `error` when asked for the commit to clone."""
    world.github_says(monkeypatch)

    def raises(*arguments):
        raise error
    monkeypatch.setattr(build, 'get_submodule_commit', raises)
    monkeypatch.setattr(subprocess, 'Popen', lambda *a, **k: pytest.fail("the backend must not be cloned without its commit"))
    return world.sources('sdist', git_repository=False)


NAME = build.GALAXIA_SUBMODULE_NAME


@pytest.mark.parametrize('error', [
    urllib.error.URLError(socket.gaierror(-2, 'Name or service not known')),
    urllib.error.URLError(ssl.SSLCertVerificationError('CERTIFICATE_VERIFY_FAILED')),
    ConnectionResetError('Connection reset by peer'),
    ConnectionRefusedError('Connection refused'),
    TimeoutError('timed out'),
    socket.timeout('timed out'),
    ssl.SSLError('handshake failure'),
    http.client.IncompleteRead(b'abc', 10),
    http.client.RemoteDisconnected('Remote end closed connection without response'),
], ids=['URLError-gaierror', 'URLError-certificate', 'ConnectionResetError', 'ConnectionRefusedError', 'TimeoutError',
        'socket.timeout', 'SSLError', 'IncompleteRead', 'RemoteDisconnected'])
def test_a_github_that_cannot_be_reached_is_a_compile_error_that_says_what_to_check(world, monkeypatch, error):
    root = lookup_that_raises(world, monkeypatch, error)
    with pytest.raises(CompileError) as raised:
        build.clone_and_checkout_submodules(root, [NAME])
    lines = str(raised.value).split('\n')
    assert lines[0].startswith(f"GitHub could not be reached ({type(error).__name__}: ")
    assert "api.github.com and github.com" in lines[0]
    assert "https_proxy" in lines[0] and "SSL_CERT_FILE or SSL_CERT_DIR" in lines[0]
    assert lines[1:] == [f"Error finding the commit of {NAME} to clone, aborting...", ""]


@pytest.mark.parametrize('code, reason', [(403, 'Forbidden'), (404, 'Not Found'), (500, 'Internal Server Error')])
def test_a_refusal_of_github_is_a_compile_error_that_says_what_it_answered_and_for_what(world, monkeypatch, code, reason):
    root = lookup_that_raises(world, monkeypatch, urllib.error.HTTPError('https://api.github.com', code, reason, {}, None))
    with pytest.raises(CompileError) as raised:
        build.clone_and_checkout_submodules(root, [NAME])
    assert str(raised.value) == (f"GitHub answered {code} ({reason}) when asked for the commit of {NAME} that the revision "
                                 f"{REVISION} of this package pins\n"
                                 f"Error finding the commit of {NAME} to clone, aborting...\n")


def test_an_answer_that_is_not_json_is_a_compile_error_that_says_that_something_else_may_have_answered(world, monkeypatch):
    with pytest.raises(json.JSONDecodeError) as bad:
        json.loads('<html>Sign in to the network</html>')
    root = lookup_that_raises(world, monkeypatch, bad.value)
    with pytest.raises(CompileError) as raised:
        build.clone_and_checkout_submodules(root, [NAME])
    lines = str(raised.value).split('\n')
    assert lines[0].startswith(f"The answer to the request for the commit of {NAME} is not what GitHub sends (Expecting value")
    assert lines[0].endswith("a proxy or a portal may have answered in its place")
    assert lines[1:] == [f"Error finding the commit of {NAME} to clone, aborting...", ""]


@pytest.mark.parametrize('error', [FileNotFoundError("Submodule path 'Galaxia-ananke' not found in repository"),
                                   ValueError("Path 'Galaxia-ananke' is not a submodule (type: tree)")],
                         ids=['FileNotFoundError', 'ValueError'])
def test_what_is_not_about_reaching_github_is_raised_as_it_is(world, monkeypatch, error):
    root = lookup_that_raises(world, monkeypatch, error)
    with pytest.raises(type(error)) as raised:
        build.clone_and_checkout_submodules(root, [NAME])
    assert raised.value is error


if __name__ == '__main__':
    pass
