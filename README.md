# py-Galaxia-ananke

[![PyPI - Version](https://img.shields.io/pypi/v/galaxia_ananke)](https://pypi.org/project/galaxia-ananke/)
[![Documentation Status](https://readthedocs.org/projects/py-galaxia-ananke/badge/?version=latest)](https://py-galaxia-ananke.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/501369954.svg)](https://zenodo.org/badge/latestdoi/501369954)
[![License](https://img.shields.io/badge/license-AGPL--3.0--or--later-blue.svg)](https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE)
[![Slack](https://img.shields.io/badge/slack-chat-green.svg)](https://join.slack.com/t/ananke-users/shared_invite/zt-37wmyn7ki-kmroxUul2W_VdFlEt0qCig)

Python wrapper for a modified version of Galaxia ([Sharma et al. 2011](http://ascl.net/1101.007)), to be used by [py-ananke](https://github.com/athob/py-ananke).

## Getting started

`py-Galaxia-ananke` is compatible with Python versions 3.9 to 3.12 (see Requirements below). The project is organized into three branches: [main](https://github.com/athob/py-Galaxia-ananke/tree/main), [stable](https://github.com/athob/py-Galaxia-ananke/tree/stable), and [develop](https://github.com/athob/py-Galaxia-ananke/tree/develop). The main branch contains the latest released version, while the stable and develop branches host versions currently in development, with stable being the most recent stable version. `py-Galaxia-ananke` is linked to a separate repository hosting the C++ backend software, [`Galaxia-ananke`](https://github.com/athob/Galaxia-ananke), a modified version of [`Galaxia`](http://ascl.net/1101.007). It is worth noting that [`Galaxia-ananke`](https://github.com/athob/Galaxia-ananke) incorporates several pre-installed photometric systems, represented by sets of isochrones generated from the [CMD web interface](http://stev.oapd.inaf.it/cgi-bin/cmd) (commonly referred to as Padova isochrones). Among the available options are HST, GAIA, Euclid, Rubin, JWST & Roman.

### Installation

`py-Galaxia-ananke` is available on the PyPI, so it may be installed using the command:

```bash
pip install galaxia_ananke
```

Alternatively, if you wish to run the latest version on the repository, you can use the following pip command, which pulls it directly from the repository's main branch:

```bash
pip install git+https://github.com/athob/py-Galaxia-ananke@main
```

You may also change the branch to use in the above command by replacing the `main` that follows the `@` symbol. If you prefer, you may clone the repository to your local machine and then install `py-Galaxia-ananke` using the following pip command, which installs it from your local copy of the repository:

```bash
git clone https://github.com/athob/py-Galaxia-ananke
cd py-Galaxia-ananke
pip install .
```

Please note that the command with flag `pip install . --no-cache-dir` may be necessary due to some dependencies issues.

<!-- ***Warning: DO NOT download the repository as a ZIP archive with intention to install it this way, the installation requires the git set up of the repository to propertly install its submodule dependencies.*** -->

After installation, the module can be imported in Python under the name `galaxia_ananke` and be ran as such.

### Requirements

Python 3.9 to 3.12 (the vaex dependency does not support newer versions yet), on a Unix system: Windows is not supported. Installing `py-Galaxia-ananke` builds its C++ backend, [`Galaxia-ananke`](https://github.com/athob/Galaxia-ananke), but the package does not carry that source: the installation fetches it from GitHub (it is a git submodule of this repository, pinned to one commit). It therefore needs network access and git, as well as GNU make (4.2 or later) and g++ (8.5 or later) to build it. The `Galaxia-ananke` repository is large, about 5.8 GB at the time of writing (most of it is data, such as the isochrones), and the installation copies its data into a cache directory under `sys.prefix`: make sure that there is enough free disk space before installing.

### Troubleshooting installation

You may find yourself in a situation after installation where importing the package module errors out in an `AssertionError`. The installation compiles and installs the backend C++ submodule Galaxia-ananke which is required, this `AssertionError` means that process failed in some way at installation. When installing the Galaxia-ananke submodule, `galaxia_ananke`'s setup write log files in a cache location. The `AssertionError` at import that calls for the missing Galaxia executable gives the `bin` path where that executable should be located. The parent directory for that `bin` path should contain also a `log` directory, where those log files can be found and can help troubleshooting the missing executable. Below are some potential situations:

#### Galaxia-ananke submodule didn't pull appropriately

The installation of `galaxia_ananke` is supposed to automatically pull the Galaxia-ananke git submodule. However, if the directory of that submodule is empty, it means that the pull failed. Try to manually run `git submodule update --init` from the root of this repository before installing.

#### build-aux/install-sh: Permission denied

In the log files, check if `Galaxia-make-install.log` contains a mention regarding a file named `build-aux/install-sh` with permission denied. This file is an executable ran by the installer, and it may need executable permission to be ran. It is located in the `Galaxia-ananke` submodule.

#### no writing permission in sys.prefix directory

When you run in a python terminal the following `import sys; sys.prefix`, the resulting path is the path of the directory where the Galaxia-ananke cached data is meant to be stored. If this directory doesn't have write permission, the installation will not complete. It is ideal to let the installation use that directory, so troubleshooting that missing write permission should be the priority. That said in last resort, it is possible to set a custom prefix directory by exporting its full path in the environment variable `ANANKE_SYSTEM_PREFIX`.

## Usage

Visit the [jupyter](https://github.com/athob/py-Galaxia-ananke/tree/main/jupyter) directory for example uses, and the [documentation](https://py-galaxia-ananke.readthedocs.io) for the API reference.

## Citing

If you use py-Galaxia-ananke in your work, please cite, in this order:

1. The py-ananke article ([Thob et al. 2024](https://doi.org/10.21105/joss.06234)): this is the main citation. GitHub's "Cite this repository" button gives it in the usual formats, from [`CITATION.cff`](https://github.com/athob/py-Galaxia-ananke/blob/main/CITATION.cff).
2. The Zenodo record of py-Galaxia-ananke itself ([10.5281/zenodo.13942945](https://doi.org/10.5281/zenodo.13942945)). The article alone could suggest that the whole py-ananke pipeline was used; this one makes clear that py-Galaxia-ananke was the software. That DOI stands for all versions; each release also has a DOI of its own on Zenodo, which is the one to cite to name the version you used.
3. Galaxia ([Sharma et al. 2011](https://doi.org/10.1088/0004-637X/730/1/3), also [ASCL 1101.007](http://ascl.net/1101.007)), of which Galaxia-ananke, the C++ code this package wraps, is a modified version.

## Contributing

Bug reports and feature requests go to the [issue tracker](https://github.com/athob/py-Galaxia-ananke/issues), which has a form for each. Contributions are welcome: [CONTRIBUTING.md](https://github.com/athob/py-Galaxia-ananke/blob/main/CONTRIBUTING.md) says how to set up a copy, what the automatic checks do and how releases work, and everyone taking part is expected to follow the [code of conduct](https://github.com/athob/py-Galaxia-ananke/blob/main/CODE_OF_CONDUCT.md). Please report security problems privately, as described in [SECURITY.md](https://github.com/athob/py-Galaxia-ananke/blob/main/SECURITY.md).

## License

py-Galaxia-ananke is licensed under the GNU Affero General Public License v3.0 or later: see [LICENSE](https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE) and [COPYRIGHT](https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT). Galaxia-ananke, the C++ code the package builds, is licensed under the GNU Affero General Public License v3.0 too: see [its repository](https://github.com/athob/Galaxia-ananke).
