#!/usr/bin/env python
"""
Contains the galaxia_ananke module building utility tools. Credit to
https://github.com/GalacticDynamics-Oxford/Agama/blob/master/setup.py.
"""
import os
import pathlib
import shutil
import sys
import subprocess
import ssl
import urllib.request
import urllib.parse
import json
import configparser
import tempfile
from distutils.errors import CompileError
from setuptools.command.build_ext import build_ext
from setuptools import Command
from packaging.version import Version

from ._builtin_utils import get_version_of_command
from ._constants import *
from .__metadata__ import *
from . import _version, versioneer

__all__ = ['make_package_data', 'make_cmdclass']

ROOT_DIR = pathlib.Path(__file__).parent.parent
MIN_GPP_VERSION = Version("8.5")
MIN_MAKE_VERSION = Version("4.2")


# force printing to the terminal even if stdout was redirected
def say(text):
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
def all_files(*paths, basedir='.'):
    basedir = pathlib.Path(basedir)
    return [str(pathlib.Path(dirpath, f).relative_to(basedir))
            for path in paths
            for dirpath, dirnames, files in os.walk(basedir / path)
            for f in files]


def verify_system_dependencies():
    try:
        git_version = Version(get_version_of_command("git"))
    except FileNotFoundError:
        raise OSError("Your system does not have git installed. Please install git before proceeding")
    try:
        make_version = Version(get_version_of_command("make"))
    except FileNotFoundError:
        raise OSError("Your system does not have the utility gnumake installed. Please install one before proceeding")
    if make_version < MIN_MAKE_VERSION:
        raise OSError(f"Your system has gnumake v{make_version} installed, but galaxia_ananke requires v{MIN_MAKE_VERSION}")
    try:
        gpp_version = Version(get_version_of_command("g++"))
    except FileNotFoundError:
        raise OSError("Your system does not have a C++ compiler installed. Please install one before proceeding")
    if gpp_version < MIN_GPP_VERSION:
        raise OSError(f"Your system has g++ v{gpp_version} installed, but galaxia_ananke requires v{MIN_GPP_VERSION}")


def download_galaxia(galaxia_dir):
    say("\nDownloading Galaxia")
    tarfile = galaxia_dir.with_suffix('.tar.gz')
    try:
        urllib.request.urlretrieve(GALAXIA_URL, filename=tarfile)
        if tarfile.is_file():
            say("\nUnpacking Galaxia")
            subprocess.call(['tar', 'xzvf', tarfile, '-C', tarfile.parent],
                            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            tarfile.unlink()
            if not galaxia_dir.is_dir():
                raise RuntimeError("Error unpacking Galaxia")
        else:
            raise RuntimeError("Cannot find downloaded file")
    except RuntimeError as e:
        raise CompileError(str(e) + "\nError downloading Galaxia, aborting...\n")


def get_submodule_commit(repo_owner, repo_name, commit_sha, submodule_path):
    """
    Get submodule commit hash for a specific GitHub repository commit.

    Args:
        repo_owner: Repository owner
        repo_name: Repository name
        commit_sha: Commit SHA in the main repository
        submodule_path: Relative path to the submodule

    Returns:
        Submodule commit SHA string
    """
    # Create SSL context (bypass verification if needed)
    ctx = ssl.create_default_context()
    # Uncomment below line if you encounter certificate errors
    # ctx.check_hostname = False
    # ctx.verify_mode = ssl.CERT_NONE
    # Step 1: Get repository tree recursively for the specific commit
    tree_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/git/trees/{commit_sha}?recursive=1"
    ################################################################
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
    ################################################################
    # response = requests.get(tree_url)
    # if response.status_code == 403 and 'rate limit' in response.text.lower():
    #     raise RuntimeError("GitHub API rate limit exceeded. Try again later.")
    # response.raise_for_status()
    # tree_data = response.json()
    ################################################################
    # Step 2: Find the submodule entry in the tree
    for item in tree_data.get('tree', []):
        if item['path'] == submodule_path:
            if item['type'] == 'commit':
                return item['sha']
            raise ValueError(f"Path '{submodule_path}' is not a submodule (type: {item['type']})")
    raise FileNotFoundError(f"Submodule path '{submodule_path}' not found in repository")


def clone_and_checkout_submodules(root_dir, submodule_names):
    parsed_url = urllib.parse.urlparse(__url__)
    url_path = pathlib.Path(parsed_url.path)
    gitmodules = configparser.ConfigParser()
    gitmodules.read(root_dir / '.gitmodules')
    for submodule_name in submodule_names:
        galaxia_url = parsed_url._replace(path=str((url_path / gitmodules[f'submodule "{submodule_name}"']['url']).resolve())).geturl()
        galaxia_commit = get_submodule_commit(
            *(url_path.parts[-2:] + (_version.get_versions()['full-revisionid'], submodule_name))
            )
        try:
            _temp = subprocess.Popen(['git', 'clone', galaxia_url], cwd=root_dir,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            for line in iter(_temp.stdout.readline, ""):
                say(line)
            _temp.wait()
        except RuntimeError as e:
            raise CompileError(str(e) + f"\nError cloning {submodule_name}, aborting...\n")
        _temp = subprocess.call(['git', 'checkout', galaxia_commit], cwd=root_dir / submodule_name)
        # git archive -vvv --remote=<REPO_URL> <REF> | tar -xvf -


def check_galaxia_submodule(root_dir):
    root_dir = pathlib.Path(root_dir)
    # if not os.listdir(GALAXIA_SUBMODULE_NAME):
    say("\nChecking submodule Galaxia, running git...")
    try:
        _temp = subprocess.call(['git', 'submodule', 'update', '--init', '--recursive'], cwd=root_dir)
    except FileNotFoundError:
        raise OSError("Your system does not have git installed. Please install git before proceeding")
    if _temp == 128:
        say("\n\tFailed to git submodule init Galaxia, attempting to clone...")
        clone_and_checkout_submodules(root_dir, [GALAXIA_SUBMODULE_NAME])
        say("\n\tClone was successful")
        # raise OSError(f"The repository from which you are attempting to install this package is not a git repository.\nPlease follow the online instructions for proper installation ({__url__}/#installation).")
    install_sh_path = root_dir / GALAXIA_SUBMODULE_NAME / 'build-aux' / 'install-sh'
    if not os.access(install_sh_path, os.X_OK):
        raise PermissionError(f"Installation cannot complete: to proceed, please give user-execute permission to file {install_sh_path}")


def remove_existing_galaxia(temp_photocat: pathlib.Path):
    if CACHE.is_dir():
        custom_photocat = ISOCHRONES_PATH / CUSTOM_PHOTOCAT
        if custom_photocat.is_dir():
            custom_photocat.rename(temp_photocat)
        shutil.rmtree(CACHE)


def touch_autotools_files(galaxia_dir: pathlib.Path):
    (galaxia_dir / 'configure').touch()
    (galaxia_dir / 'config.h.in').touch()
    (galaxia_dir / 'Makefile.in').touch()
    (galaxia_dir / 'aclocal.m4').touch()
    (galaxia_dir / 'build-aux' / 'install-sh').touch()


def configure_galaxia(galaxia_dir: pathlib.Path):
    with (GALAXIA_LOG / 'Galaxia-configure.log').open('w') as f:
        subprocess.call([f"./configure",
                        f"--prefix={CACHE}",
                        f"--datadir={GALAXIA_DATA}"],
                        cwd=galaxia_dir, stdout=f, stderr=f)


def make_galaxia(galaxia_dir: pathlib.Path):
    with (GALAXIA_LOG / 'Galaxia-make.log').open('w') as f:
        subprocess.call(["make"],
                        cwd=galaxia_dir, stdout=f, stderr=f)


def make_install_galaxia(galaxia_dir: pathlib.Path):
    with (GALAXIA_LOG / 'Galaxia-make-install.log').open('w') as f:
        subprocess.call(["make", "install"],
                        cwd=galaxia_dir, stdout=f, stderr=f)
    shutil.copytree(galaxia_dir / GALAXIA_DATA.name, GALAXIA_DATA)


def make_distclean_galaxia(galaxia_dir: pathlib.Path):
    with (GALAXIA_LOG / 'Galaxia-make-distclean.log').open('w') as f:
        subprocess.call(["make", "distclean"],
                        cwd=galaxia_dir, stdout=f, stderr=f)


def clean_up_temporary(temp_photocat: pathlib.Path):
    if temp_photocat.is_dir():
        temp_photocat.rename(ISOCHRONES_PATH / CUSTOM_PHOTOCAT)


def build_and_install_galaxia(galaxia_dir: pathlib.Path):
    galaxia_dir = pathlib.Path(galaxia_dir).resolve()
    temp_dir = tempfile.TemporaryDirectory()
    temp_photocat = pathlib.Path(temp_dir.name) / CUSTOM_PHOTOCAT
    remove_existing_galaxia(temp_photocat)
    say("\nBuilding Galaxia")
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
    say("\n\tCleaning temporary")
    clean_up_temporary(temp_photocat)
    say("\n")


def make_package_data():
    for_all_files = ('__license__', )
    return {NAME: all_files(*for_all_files,
                            basedir=pathlib.Path(SRC_DIR, NAME))}


def make_cmdclass():
    """
    """
    class _build_ext(build_ext):
        def run(self):
            build_ext.run(self)
            verify_system_dependencies()
            check_galaxia_submodule(ROOT_DIR)
            build_and_install_galaxia(GALAXIA_SUBMODULE_NAME)


    class _test(Command):
        description = 'run tests'
        user_options = []

        def initialize_options(self): pass

        def finalize_options(self): pass

        def run(self): pass

    return versioneer.get_cmdclass({'build_ext': _build_ext, 'test': _test})
