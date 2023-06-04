# py-Galaxia-ananke

Python wrapper for a modified version of Galaxia ([Sharma et al. 2011](http://ascl.net/1101.007)).

## Getting started

The Python package py-Galaxia-ananke requires galaxia-ananke which is included in this repository as a submodule. After cloning this repository and before any installation attempt, make sure to run `git submodule update --init` to pull the remote repository content. The package py-Galaxia-ananke can eventually be installed using `pip install .` from the root of this repository. The command with flag `pip install . --no-cache-dir` may be necessary.

After installation, the module can be imported in Python as `Galaxia_ananke` and be ran as such.

## Troubleshooting installation

You may find yourself in a situation after installation where importing the package module errors out in an `AssertionError`. The installation compiles and installs the backend C++ submodule galaxia-ananke which is required, this `AssertionError` means that process failed in some way at installation. When installing the galaxia-ananke submodule, `Galaxia_ananke`'s setup write log files in a cache location. The `AssertionError` at import that calls for the missing Galaxia executable gives the `bin` path where that executable should be located. The parent directory for that `bin` path should contain also a `log` directory, where those log files can be found and can help troubleshooting the missing executable. Below are some potential situations:

### galaxia-ananke submodule didn't pull appropriately

The installation of `Galaxia_ananke` is supposed to automatically pull the galaxia-ananke git submodule. However, if the directory of that submodule is empty, it means that the pull failed. Try to manually run `git submodule update --init` from the root of this repository before installing.

### build-aux/install-sh: Permission denied

In the log files, check if `Galaxia-make-install.log` contains a mention regarding a file named `build-aux/install-sh` with permission denied. This file is an executable ran by the installer, and it may need executable permission to be ran. It is located in the `galaxia-ananke` submodule.

### no writing permission in sys.prefix directory

When you run in a python terminal the following `import sys; sys.prefix`, the resulting path is the path of the directory where the galaxia-ananke cached data is meant to be stored. If this directory doesn't have write permission, the installation will not complete. It is ideal to let the installation use that directory, so troubleshooting that missing write permission should be the priority. That said in last resort, it is possible to set a custom prefix directory by exporting its full path in the environment variable `ANANKE_SYSTEM_PREFIX`.
