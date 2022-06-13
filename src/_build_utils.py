#!/usr/bin/env python
"""
Docstring
"""
import pathlib
import shutil
import sys
import subprocess
import urllib.request
from distutils.errors import CompileError

from .Galaxia.constants import *

__all__ = ['say', 'all_files', 'download_galaxia', 'check_galaxia_submodule', 'build_and_install_galaxia']

ROOT_DIR = pathlib.Path(__file__).parent.parent


# force printing to the terminal even if stdout was redirected
def say(text):
    text += ' '
    sys.stdout.write(text)
    sys.stdout.flush()
    if not sys.stdout.isatty():
        # output was redirected, but we still try to send the message to the terminal
        try:
            with open('/dev/tty', 'w') as out:
                out.write(text)
                out.flush()
        except PermissionError:
            # /dev/tty may not exist or may not be writable!
            pass


# get the list of all files in the given directories (including those in nested directories)
def all_files(*paths, basedir='.'):
    basedir = pathlib.Path(basedir)
    return [str(pathlib.Path(dirpath, f).relative_to(basedir))
            for path in paths
            for dirpath, dirnames, files in pathlib.os.walk(basedir / path)
            for f in files]


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


def check_galaxia_submodule(root_dir):
    if not pathlib.os.listdir(GALAXIA_SUBMODULE_NAME):
        say("\nEmpty submodule Galaxia, running git...")
        subprocess.call(['git', 'submodule', 'init'], cwd=root_dir)
        subprocess.call(['git', 'submodule', 'update'], cwd=root_dir)


def remove_existing_galaxia():
    shutil.rmtree(CACHE)


def configure_galaxia(galaxia_dir):
    with (GALAXIA_LOG / 'Galaxia-configure.log').open('w') as f:
        subprocess.call([f"./configure",
                        f"--prefix={CACHE}",
                        f"--datadir={GALAXIA_DATA}"],
                        cwd=galaxia_dir, stdout=f, stderr=f)


def make_galaxia(galaxia_dir):
    with (GALAXIA_LOG / 'Galaxia-make.log').open('w') as f:
        subprocess.call(["make"],
                        cwd=galaxia_dir, stdout=f, stderr=f)


def make_install_galaxia(galaxia_dir):
    with (GALAXIA_LOG / 'Galaxia-make-install.log').open('w') as f:
        subprocess.call(["make", "install"],
                        cwd=galaxia_dir, stdout=f, stderr=f)
    shutil.copytree(galaxia_dir / GALAXIA_DATA.name, GALAXIA_DATA)


def make_distclean_galaxia(galaxia_dir):
    with (GALAXIA_LOG / 'Galaxia-make-distclean.log').open('w') as f:
        subprocess.call(["make", "distclean"],
                        cwd=galaxia_dir, stdout=f, stderr=f)


def build_and_install_galaxia(galaxia_dir):
    galaxia_dir = pathlib.Path(galaxia_dir).resolve()
    remove_existing_galaxia()
    say("\nBuilding Galaxia")
    GALAXIA.parent.mkdir(parents=True, exist_ok=True)
    GALAXIA_LOG.mkdir(parents=True, exist_ok=True)
    say("\n\tConfiguring")
    configure_galaxia(galaxia_dir)
    say("\n\tRunning make")
    make_galaxia(galaxia_dir)
    say("\n\tRunning make install")
    make_install_galaxia(galaxia_dir)
    say("\n\tRunning make distclean")
    make_distclean_galaxia(galaxia_dir)
    say("\n")
