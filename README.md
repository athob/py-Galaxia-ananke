# py-Galaxia-ananke

[![Documentation Status](https://readthedocs.org/projects/py-galaxia-ananke/badge/?version=latest)](https://py-galaxia-ananke.readthedocs.io/en/latest/?badge=latest)[![DOI](https://zenodo.org/badge/501369954.svg)](https://zenodo.org/badge/latestdoi/501369954)

Python wrapper for a modified version of Galaxia ([Sharma et al. 2011](http://ascl.net/1101.007)).

## Getting started

`py-Galaxia-ananke` is compatible with Python versions above 3.7.12 and below 3.11. The project is organized into three branches: [main](https://github.com/athob/py-Galaxia-ananke/tree/main), [stable](https://github.com/athob/py-Galaxia-ananke/tree/stable), and [develop](https://github.com/athob/py-Galaxia-ananke/tree/develop). The main branch contains the latest released version, while the stable and develop branches host versions currently in development, with stable being the most recent stable version. `py-Galaxia-ananke` is linked to a separate repository hosting the C++ backend software, [`galaxia-ananke`](https://github.com/athob/galaxia-ananke), a modified version of [`Galaxia`](http://ascl.net/1101.007). It is worth noting that [`galaxia-ananke`](https://github.com/athob/galaxia-ananke) incorporates several pre-installed photometric systems, represented by sets of isochrones generated from the [CMD web interface](http://stev.oapd.inaf.it/cgi-bin/cmd) (commonly referred to as Padova isochrones). Among the available options are HST, GAIA, Euclid, Rubin, JWST & Roman.

### Installation

To install `py-Galaxia-ananke`, you can use the following pip command, which pulls the latest version directly from the repository's main branch:

    pip install git+https://github.com/athob/py-Galaxia-ananke@main

or

    python -m pip install git+https://github.com/athob/py-Galaxia-ananke@main

Alternatively, if you prefer, you may clone the repository to your local machine and then install `py-Galaxia-ananke` using the following pip command, which installs it from your local copy of the repository:

    git clone https://github.com/athob/py-Galaxia-ananke
    cd py-Galaxia-ananke
    pip install .

Please note that the command with flag `pip install . --no-cache-dir` may be necessary due to some dependencies issues.

***Warning: DO NOT download the repository as a ZIP archive with intention to install it this way, the installation requires the git set up of the repository to propertly install its submodule dependencies.***

After installation, the module can be imported in Python under the name `galaxia_ananke` and be ran as such.

### Troubleshooting installation

You may find yourself in a situation after installation where importing the package module errors out in an `AssertionError`. The installation compiles and installs the backend C++ submodule galaxia-ananke which is required, this `AssertionError` means that process failed in some way at installation. When installing the galaxia-ananke submodule, `galaxia_ananke`'s setup write log files in a cache location. The `AssertionError` at import that calls for the missing Galaxia executable gives the `bin` path where that executable should be located. The parent directory for that `bin` path should contain also a `log` directory, where those log files can be found and can help troubleshooting the missing executable. Below are some potential situations:

#### galaxia-ananke submodule didn't pull appropriately

The installation of `galaxia_ananke` is supposed to automatically pull the galaxia-ananke git submodule. However, if the directory of that submodule is empty, it means that the pull failed. Try to manually run `git submodule update --init` from the root of this repository before installing.

#### build-aux/install-sh: Permission denied

In the log files, check if `Galaxia-make-install.log` contains a mention regarding a file named `build-aux/install-sh` with permission denied. This file is an executable ran by the installer, and it may need executable permission to be ran. It is located in the `galaxia-ananke` submodule.

#### no writing permission in sys.prefix directory

When you run in a python terminal the following `import sys; sys.prefix`, the resulting path is the path of the directory where the galaxia-ananke cached data is meant to be stored. If this directory doesn't have write permission, the installation will not complete. It is ideal to let the installation use that directory, so troubleshooting that missing write permission should be the priority. That said in last resort, it is possible to set a custom prefix directory by exporting its full path in the environment variable `ANANKE_SYSTEM_PREFIX`.
