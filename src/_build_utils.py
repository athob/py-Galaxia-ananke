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
Contains the galaxia_ananke module building utility tools. Credit to
https://github.com/GalacticDynamics-Oxford/Agama/blob/master/setup.py.
"""
import platform
import os
import pathlib
import shutil
import sys
import subprocess
import socket
import ssl
import http.client
import urllib.error
import urllib.request
import urllib.parse
import json
import configparser
from typing import Any, ClassVar, Dict, IO, List, Optional, Tuple, Union, cast
from warnings import warn
from distutils.errors import CompileError
from setuptools.command.build_py import build_py
from setuptools import Command
from packaging.version import Version

from ._builtin_utils import get_version_of_command
from ._constants import *
from ._constants import PREFIX_ENV_VAR
from .__metadata__ import *
from . import _version, versioneer

__all__ = ['make_package_data', 'make_cmdclass']

if platform.system() == "Windows":
    raise OSError(f"Windows compatibility is not currently supported by {NAME}. We apologize for the inconvenience.")

ROOT_DIR = pathlib.Path(__file__).parent.parent
MIN_GPP_VERSION = Version("8.5")
MIN_MAKE_VERSION = Version("4.2")
# what a request to GitHub raises when GitHub cannot be reached (an HTTPError, which is a URLError too, is GitHub that
# answers: it comes first). socket.timeout is not the TimeoutError of Python 3.9, and is the same one from Python 3.10
GITHUB_UNREACHABLE = (urllib.error.URLError, ConnectionError, TimeoutError, socket.timeout, ssl.SSLError,
                      http.client.HTTPException)


# force printing to the terminal even if stdout was redirected
def say(text: str) -> None:
    """
    Write a text, followed by a space, to the standard output and, if that is redirected, to the terminal too.

    The terminal is the one of ``/dev/tty``, if there is one that can be written to: nothing is said if there is none.

    Parameters
    ----------
    text : str
        The text.
    """
    text += ' '
    sys.stdout.write(text)
    sys.stdout.flush()
    if not sys.stdout.isatty():
        # output was redirected, but we still try to send the message to the terminal
        try:
            if pathlib.Path('/dev/tty').exists():
                with open('/dev/tty', 'w') as out:
                    out.write(text)
                    out.flush()
        except (OSError, PermissionError):
            # /dev/tty may not exist or may not be writable!
            pass


# get the list of all files in the given directories (including those in nested directories)
def all_files(*paths: str, basedir: Union[str, pathlib.Path] = '.') -> List[str]:
    """
    List the files of directories, those in nested directories included.

    Parameters
    ----------
    *paths : str
        The directories, relative to ``basedir``.
    basedir : str or path, optional
        Where the paths are relative to. Default to the working directory.

    Returns
    -------
    files : list of str
        The files, as paths relative to ``basedir``. A directory that does not exist has none.
    """
    basedir = pathlib.Path(basedir)
    return [str(pathlib.Path(dirpath, f).relative_to(basedir))
            for path in paths
            for dirpath, dirnames, files in os.walk(basedir / path)
            for f in files]


def verify_system_dependencies() -> None:
    """
    Check that the system has what building the backend needs: git, GNU make and a C++ compiler.

    The versions have to be at least ``MIN_MAKE_VERSION`` for make and ``MIN_GPP_VERSION`` for g++.

    Raises
    ------
    OSError
        If a tool is missing, or too old: the message says which, and what to do.
    """
    try:
        git_version = Version(get_version_of_command("git"))
    except FileNotFoundError:
        raise OSError("Your system does not have git installed. Please install git before proceeding")
    try:
        make_version = Version(get_version_of_command("make"))
    except FileNotFoundError:
        raise OSError("Your system does not have the utility gnumake installed. Please install one before proceeding")
    if make_version < MIN_MAKE_VERSION:
        raise OSError(f"Your system has gnumake v{make_version} installed, but {NAME} requires v{MIN_MAKE_VERSION}")
    try:
        gpp_version = Version(get_version_of_command("g++"))
    except FileNotFoundError:
        raise OSError("Your system does not have a C++ compiler installed. Please install one before proceeding")
    if gpp_version < MIN_GPP_VERSION:
        raise OSError(f"Your system has g++ v{gpp_version} installed, but {NAME} requires v{MIN_GPP_VERSION}")


def verify_cache_can_be_written() -> None:
    """
    Refuse, before anything is built, a cache that the installation cannot make. The earlier installation is put aside next to
    the cache and the new one is built in its place, so the directory that the cache is in has to be writable: the nearest one
    that exists, if the cache is not made yet.
    """
    directory = CACHE.parent
    while not directory.exists() and directory != directory.parent:
        directory = directory.parent
    if not os.access(directory, os.W_OK | os.X_OK):
        # TODO make the variable conda-dependent?
        raise PermissionError("Installation cannot complete: to proceed, please give write permission to directory "
                              f"{directory} or define a custom system prefix via the environment variable {PREFIX_ENV_VAR}")


def get_submodule_commit(repo_owner: str, repo_name: str, commit_sha: str, submodule_path: str) -> str:
    """
    Give the commit of a submodule that a commit of a GitHub repository pins.

    The repository tree of that commit is asked of the GitHub API.

    Parameters
    ----------
    repo_owner : str
        The owner of the repository.
    repo_name : str
        The name of the repository.
    commit_sha : str
        The commit of the repository.
    submodule_path : str
        The path of the submodule in the repository.

    Returns
    -------
    commit : str
        The SHA of the commit of the submodule.

    Raises
    ------
    RuntimeError
        If the rate limit of the GitHub API is exceeded.
    urllib.error.HTTPError
        If GitHub refuses the request for another reason.
    ValueError
        If the path is in the repository, but is not a submodule.
    FileNotFoundError
        If there is nothing at the path of the submodule.
    """
    # Create SSL context
    ctx = ssl.create_default_context()
    # Step 1: Get repository tree recursively for the specific commit
    tree_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/git/trees/{commit_sha}?recursive=1"
    try:
        # Make HTTP request
        req = urllib.request.Request(tree_url)
        with urllib.request.urlopen(req, context=ctx) as response:
            # Check rate limit headers
            remaining = response.headers.get('X-RateLimit-Remaining')
            if remaining == '0':
                reset_time = response.headers.get('X-RateLimit-Reset')
                raise RuntimeError(f"GitHub API rate limit exceeded. Reset at {reset_time}")
            # Parse JSON response
            data = response.read()
            encoding = response.info().get_content_charset('utf-8')
            tree_data = json.loads(data.decode(encoding))
    except urllib.error.HTTPError as e:
        if e.code == 403 and 'rate limit' in e.reason.lower():
            raise RuntimeError("GitHub API rate limit exceeded. Try again later.")
        raise  # Re-raise other HTTP errors
    # Step 2: Find the submodule entry in the tree
    for item in tree_data.get('tree', []):
        if item['path'] == submodule_path:
            if item['type'] == 'commit':
                return item['sha']
            raise ValueError(f"Path '{submodule_path}' is not a submodule (type: {item['type']})")
    raise FileNotFoundError(f"Submodule path '{submodule_path}' not found in repository")


def commit_lookup_error(reason: str, submodule_name: str) -> CompileError:
    """
    Make the error of the build for when the commit of a submodule to clone cannot be found, whatever the reason.

    Parameters
    ----------
    reason : str
        What went wrong, which starts the message.
    submodule_name : str
        The name of the submodule.

    Returns
    -------
    error : distutils.errors.CompileError
        The error, which is not raised.
    """
    return CompileError(f"{reason}\nError finding the commit of {submodule_name} to clone, aborting...\n")


def clone_and_checkout_submodules(root_dir: pathlib.Path, submodule_names: List[str]) -> None:
    """
    Clone submodules from GitHub, and check out the commit that this version of the package pins.

    This is for sources that are not a git repository, as an unpacked sdist is. The address of each submodule is read in
    ``.gitmodules``, and its commit is the one that the revision of this package pins, which is asked of GitHub. What git
    prints is said as it comes.

    Parameters
    ----------
    root_dir : pathlib.Path
        The root of the sources, where ``.gitmodules`` is and where the submodules are cloned.
    submodule_names : list of str
        The names of the submodules.

    Raises
    ------
    distutils.errors.CompileError
        If the commit of a submodule cannot be found (the message says why, as GitHub answered or could not be reached), or
        if the clone or the check out fail (the message gives what git said).
    """
    parsed_url = urllib.parse.urlparse(__url__)
    url_path = pathlib.Path(parsed_url.path)
    gitmodules = configparser.ConfigParser()
    gitmodules.read(root_dir / '.gitmodules')
    for submodule_name in submodule_names:
        galaxia_url = parsed_url._replace(path=str((url_path / gitmodules[f'submodule "{submodule_name}"']['url']).resolve())).geturl()
        revision = _version.get_versions()['full-revisionid']
        try:
            galaxia_commit = get_submodule_commit(*(url_path.parts[-2:] + (revision, submodule_name)))
        except RuntimeError as e:  # the rate limit of GitHub
            raise commit_lookup_error(str(e), submodule_name)
        except urllib.error.HTTPError as e:  # GitHub answers, and refuses
            raise commit_lookup_error(f"GitHub answered {e.code} ({e.reason}) when asked for the commit of {submodule_name} "
                                      f"that the revision {revision} of this package pins", submodule_name)
        except GITHUB_UNREACHABLE as e:
            raise commit_lookup_error(f"GitHub could not be reached ({type(e).__name__}: {e}): the installation needs "
                                      "api.github.com and github.com, check the network, the proxy variables (https_proxy) "
                                      "and, for a certificate error, SSL_CERT_FILE or SSL_CERT_DIR", submodule_name)
        except json.JSONDecodeError as e:  # something that is not GitHub answers
            raise commit_lookup_error(f"The answer to the request for the commit of {submodule_name} is not what GitHub sends "
                                      f"({e}): a proxy or a portal may have answered in its place", submodule_name)
        output = []
        process = subprocess.Popen(['git', 'clone', galaxia_url], cwd=root_dir,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        stdout = cast(IO[str], process.stdout)
        for line in iter(stdout.readline, ""):
            say(line)
            output.append(line.rstrip('\n'))
        stdout.close()
        status = process.wait()
        if status:
            raise CompileError('\n'.join([f"Error cloning {submodule_name} from {galaxia_url} (git exit status {status})"]
                                         + output[-20:] + ["aborting...\n"]))
        status = subprocess.call(['git', 'checkout', galaxia_commit], cwd=root_dir / submodule_name)
        if status:
            raise CompileError(f"Error checking out the commit {galaxia_commit} of {submodule_name} "
                               f"(git exit status {status}), see what git says above\naborting...\n")


def check_galaxia_submodule(root_dir: Union[str, pathlib.Path]) -> None:
    """
    Check that the sources of the backend are there, and fetch them if they are not.

    Git is asked to update the submodule. Where the sources are not a git repository, which git says with the status 128,
    the submodule is cloned instead (see ``clone_and_checkout_submodules``).

    Parameters
    ----------
    root_dir : str or path
        The root of the sources of the package.

    Raises
    ------
    OSError
        If git is not installed, or fails in another way.
    PermissionError
        If the script ``install-sh`` of the backend cannot be executed.
    distutils.errors.CompileError
        If the clone fails (see ``clone_and_checkout_submodules``).
    """
    root_dir = pathlib.Path(root_dir)
    say("\nChecking submodule Galaxia, running git...")
    try:
        _temp = subprocess.call(['git', 'submodule', 'update', '--init', '--recursive'], cwd=root_dir)
    except FileNotFoundError:
        raise OSError("Your system does not have git installed. Please install git before proceeding")
    # 128 is what git says of sources that are not a git repository, as when an sdist is unpacked, but also of other fatal
    # errors: only the first are cloned
    if _temp == 128 and not (root_dir / '.git').exists():
        say("\n\tFailed to git submodule init Galaxia, attempting to clone...")
        clone_and_checkout_submodules(root_dir, [GALAXIA_SUBMODULE_NAME])
        say("\n\tClone was successful")
    elif _temp:
        raise OSError(f"git submodule update --init --recursive failed in {root_dir} (exit status {_temp}), "
                      "see what git says above")
    install_sh_path = root_dir / GALAXIA_SUBMODULE_NAME / 'build-aux' / 'install-sh'
    if not os.access(install_sh_path, os.X_OK):
        raise PermissionError(f"Installation cannot complete: to proceed, please give user-execute permission to file {install_sh_path}")


def installation_aside() -> pathlib.Path:
    """
    Give where the earlier installation waits while the new one is built.

    It is next to the cache, so on its file system, where a directory can be renamed whatever it holds.

    Returns
    -------
    path : pathlib.Path
        The path.
    """
    return CACHE.with_name(CACHE.name + '.aside')


def installation_to_delete() -> pathlib.Path:
    """
    Give where the earlier installation goes, once the new one is complete, to be deleted: a rename cannot be left half
    done.

    Returns
    -------
    path : pathlib.Path
        The path.
    """
    return CACHE.with_name(CACHE.name + '.delete')


def restore_earlier_installation() -> None:
    """
    Put the earlier installation back in the place of what a build that stopped made, which is removed but for its logs: they
    take the place of the earlier ones, where the error of the build says that they are. Without an earlier installation, as in
    a first build, what the build made is left, logs included.
    """
    aside = installation_aside()
    if not aside.is_dir():
        return
    if GALAXIA_LOG.is_dir():
        logs = aside / GALAXIA_LOG.relative_to(CACHE)
        if logs.is_dir():
            shutil.rmtree(logs)
        logs.parent.mkdir(parents=True, exist_ok=True)
        GALAXIA_LOG.rename(logs)
    if CACHE.exists():
        shutil.rmtree(CACHE)
    aside.rename(CACHE)


def recover_from_a_stopped_build() -> None:
    """
    Put the cache back as it was before a build that was killed in the middle of its work, if there was one. The earlier
    installation waits aside during a build and is only given up when the new one is complete: one that is aside belongs to a
    build that did not finish, and is put back; one that is to delete belongs to a build that did, and is deleted.
    """
    to_delete = installation_to_delete()
    if to_delete.is_dir():
        shutil.rmtree(to_delete)
    restore_earlier_installation()


def put_earlier_installation_aside() -> None:
    """Rename the cache, if there is one, to where the earlier installation waits while the new one is built."""
    if CACHE.is_dir():
        CACHE.rename(installation_aside())


def touch_autotools_files(galaxia_dir: pathlib.Path) -> None:
    """
    Make the files that autotools generates, and make them recent.

    This is so that make does not try to remake them with autotools, which the system may not have.

    Parameters
    ----------
    galaxia_dir : pathlib.Path
        The sources of the backend.
    """
    (galaxia_dir / 'Makefile.am').touch()
    (galaxia_dir / 'configure.ac').touch()
    (galaxia_dir / 'aclocal.m4').touch()
    (galaxia_dir / 'Makefile.in').touch()
    (galaxia_dir / 'configure').touch()
    (galaxia_dir / 'config.h.in').touch()


def run_step(name: str, command: List[str], log: str, galaxia_dir: pathlib.Path, fatal: bool = True, tail: int = 20) -> None:
    """
    Run a step of the build with its output in its log.

    If it fails, an error is raised that gives the last lines of the log, or only a warning if the step is not ``fatal``.

    Parameters
    ----------
    name : str
        What the step is called in the messages.
    command : list of str
        The command, as a list of its arguments.
    log : str
        The name of the file of the log, in the folder of the logs of the backend.
    galaxia_dir : pathlib.Path
        The sources of the backend, where the command is run.
    fatal : bool, optional
        Whether a failure stops the build. Default to True.
    tail : int, optional
        The number of the last lines of the log that the error gives. Default to 20.

    Raises
    ------
    distutils.errors.CompileError
        If the command exited with an error, and the step is fatal.
    """
    with (GALAXIA_LOG / log).open('w') as f:
        status = subprocess.call(command, cwd=galaxia_dir, stdout=f, stderr=f)
    if status:
        message = f"Error running {name} for Galaxia (exit status {status}), see the log {GALAXIA_LOG / log}"
        if fatal:
            last_lines = (GALAXIA_LOG / log).read_text(errors='replace').splitlines()[-tail:]
            raise CompileError('\n'.join([message] + last_lines + ["aborting...\n"]))
        warn(message + ": the installation goes on", RuntimeWarning)


def configure_galaxia(galaxia_dir: pathlib.Path) -> None:
    """
    Run the configure script of the backend, so that it is installed in the cache, with its data in the data directory of
    the installation.

    Parameters
    ----------
    galaxia_dir : pathlib.Path
        The sources of the backend.
    """
    run_step('configure', ["./configure", f"--prefix={CACHE}", f"--datadir={GALAXIA_DATA}"],
             'Galaxia-configure.log', galaxia_dir)


def make_galaxia(galaxia_dir: pathlib.Path) -> None:
    """
    Run make in the sources of the backend.

    Parameters
    ----------
    galaxia_dir : pathlib.Path
        The sources of the backend.
    """
    run_step('make', ["make"], 'Galaxia-make.log', galaxia_dir)


def make_install_galaxia(galaxia_dir: pathlib.Path) -> None:
    """
    Run make install in the sources of the backend, and copy its data directory to the one of the installation.

    Parameters
    ----------
    galaxia_dir : pathlib.Path
        The sources of the backend.
    """
    run_step('make install', ["make", "install"], 'Galaxia-make-install.log', galaxia_dir)
    shutil.copytree(galaxia_dir / GALAXIA_DATA.name, GALAXIA_DATA)


def make_distclean_galaxia(galaxia_dir: pathlib.Path) -> None:
    """
    Run make distclean in the sources of the backend, to clean them. A failure only warns: it only cleans.

    Parameters
    ----------
    galaxia_dir : pathlib.Path
        The sources of the backend.
    """
    run_step('make distclean', ["make", "distclean"], 'Galaxia-make-distclean.log', galaxia_dir, fatal=False)  # it only cleans


def give_up_earlier_installation() -> None:
    """
    The new installation is complete: take the custom photometric systems of the earlier one, and give it up. The systems are
    copied before the earlier installation is renamed to be deleted, so that a build that is killed before the rename leaves
    them in it.
    """
    aside = installation_aside()
    if not aside.is_dir():
        return
    custom_photocat = aside / ISOCHRONES_PATH.relative_to(CACHE) / CUSTOM_PHOTOCAT
    if custom_photocat.is_dir():
        shutil.copytree(custom_photocat, ISOCHRONES_PATH / CUSTOM_PHOTOCAT, dirs_exist_ok=True)
    to_delete = installation_to_delete()
    aside.rename(to_delete)
    try:
        shutil.rmtree(to_delete)
    except OSError as error:      # the new installation is complete: the next build deletes it again
        warn(f"Could not delete the earlier installation of Galaxia in {to_delete}: {error}", RuntimeWarning)


def build_and_install_galaxia(galaxia_dir: pathlib.Path) -> None:
    """
    Build the backend, and install it in the cache.

    The earlier installation, if there is one, is put aside, and put back if the build fails or is interrupted.

    Parameters
    ----------
    galaxia_dir : pathlib.Path
        The sources of the backend.
    """
    galaxia_dir = pathlib.Path(galaxia_dir).resolve()
    recover_from_a_stopped_build()
    put_earlier_installation_aside()
    say("\nBuilding Galaxia")
    try:
        GALAXIA.parent.mkdir(parents=True, exist_ok=True)
        GALAXIA_LOG.mkdir(parents=True, exist_ok=True)
        touch_autotools_files(galaxia_dir)
        say("\n\tConfiguring")
        configure_galaxia(galaxia_dir)
        say("\n\tRunning make")
        make_galaxia(galaxia_dir)
        say("\n\tRunning make install")
        make_install_galaxia(galaxia_dir)
        say("\n\tRunning make distclean")
        make_distclean_galaxia(galaxia_dir)
        say("\n\tCleaning up")
        give_up_earlier_installation()
    except BaseException:  # a failure or an interruption: the earlier installation is put back
        restore_earlier_installation()
        raise
    say("\n")


def make_package_data() -> Dict[str, List[str]]:
    """
    Give what the package ships besides its modules: the marker of its type hints and the files of its license.

    Reading this builds nothing.

    Returns
    -------
    package_data : dict
        The files to ship, relative to the package, by the name of the package.
    """
    for_all_files = ('__license__', )
    return {NAME: ['py.typed'] + all_files(*for_all_files,
                                           basedir=pathlib.Path(SRC_DIR, NAME))}


def make_cmdclass() -> Dict[str, Any]:
    """
    Give the commands of setup.py: those of versioneer, with a ``build_py`` that builds the backend and a ``test`` that does
    nothing.

    The ``build_py`` command runs the standard one, and then checks the system and the cache, fetches the sources of the
    backend if they are not there, and builds and installs the backend.

    Returns
    -------
    cmdclass : dict
        The commands, by name.
    """
    class _build_py(build_py):
        def run(self) -> None:
            build_py.run(self)
            verify_system_dependencies()
            verify_cache_can_be_written()
            check_galaxia_submodule(ROOT_DIR)
            build_and_install_galaxia(ROOT_DIR / GALAXIA_SUBMODULE_NAME)


    class _test(Command):
        description = 'run tests'
        user_options: ClassVar[List[Tuple[str, Optional[str], str]]] = []

        def initialize_options(self) -> None: pass

        def finalize_options(self) -> None: pass

        def run(self) -> None: pass

    return versioneer.get_cmdclass({'build_py': _build_py, 'test': _test})
