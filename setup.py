#!/usr/bin/env python
import pathlib
import setuptools
from setuptools import setup

from src._build_utils import *
from src._constants import NAME, SRC_DIR
from src.__metadata__ import *

ROOT_DIR = pathlib.Path(__file__).parent

for_all_files = ('__license__', )

long_description = ""

package_data = {NAME: all_files(*for_all_files,
                                basedir=pathlib.Path(SRC_DIR, NAME))}

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
      package_data=package_data,
      include_package_data=True,
      install_requires=['numpy>=1.22,<2', 'pandas>=2,<3', 'vaex>=4.17,<5', 'astropy>=5,<7', 'h5py>=3.6,<4', 'ebfpy>=0.0.20,<1', 'astroquery>=0.4.2,<1'],
      ext_modules=[setuptools.extension.Extension('', [])],
      cmdclass=make_cmdclass(ROOT_DIR),
      )
