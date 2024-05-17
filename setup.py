#!/usr/bin/env python
import pathlib
import setuptools
from setuptools.command.build_ext import build_ext
from setuptools import setup, Command

from src._build_utils import *
from src.constants import NAME, GALAXIA_SUBMODULE_NAME, SRC_DIR
from src.__metadata__ import *

ROOT_DIR = pathlib.Path(__file__).parent

long_description = ""


class MyBuildExt(build_ext):
    def run(self):
        check_galaxia_submodule(ROOT_DIR)
        build_and_install_galaxia(GALAXIA_SUBMODULE_NAME)


class MyTest(Command):
    description = 'run tests'
    user_options = []

    def initialize_options(self): pass

    def finalize_options(self): pass

    def run(self): pass


setup(name=NAME,
      version=__version__,
      author=__author__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url=__url__,
      description=f"{__project__}: {__description__}",
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=__classifiers__,
      python_requires='>=3.8,<3.11',
      packages=[NAME, f"{NAME}.photometry"],
      package_dir={'': SRC_DIR},
      install_requires=['numpy>=1.22,<2', 'pandas>=2,<3', 'vaex>=4.17,<5', 'astropy>=5,<7', 'h5py>=3.6,<4', 'ebfpy>=0.0.20,<1', 'astroquery>=0.4.2,<1'],
      ext_modules=[setuptools.extension.Extension('', [])],
      cmdclass={'build_ext': MyBuildExt, 'test': MyTest},
      )
