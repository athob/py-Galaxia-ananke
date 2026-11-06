# Changelog

All notable changes to py-Galaxia-ananke are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project follows [Semantic Versioning](https://semver.org/) as far as `0.x` allows:
anything may change between minor versions (the sections say so when the interface did).
A section is written for every release when it is made; those below were written
afterwards, from the release notes and the history. The date of a version is the day it
was published: its upload to PyPI, which began with 0.6.0, or, for 0.4.0 and 0.5.0, its
GitHub release (in UTC). The package did not start with 0.4.0: it had been in development
since June 2022, in pre-releases that are summarized at the bottom, under
[Before 0.4.0](#before-040). The C++ backend, Galaxia-ananke, is a git submodule that has
no changelog of its own: its changes are told here, under the version that pins them, as
far as the commits of this package and of the backend say.

## [Unreleased]

### Added

- Project files: a security policy, a code of conduct, issue forms and a pull request
  template, an editorconfig, a gitattributes that normalizes the line endings, a
  pre-commit configuration, a markdownlint configuration, and Dependabot for the GitHub
  Actions.
- `CITATION.cff`, with the py-ananke article as the preferred citation, the Zenodo record
  of this package, and Galaxia. The Zenodo record now points back to the article.
- A check in CI: ruff. CI also runs on pull requests, every week and on demand.
- A README that gives the installation requirements, the size of the C++ backend among
  them, and says how to cite the package.

### Changed

- **Breaking:** the list of all the partition ids of a survey, which the output files held
  in their metadata, is a dataset of its own, `partition_id`, in a new HDF5 group,
  `Partitioning`. The attribute `All_partition_ids` now says where to find it.
- `append_hash` defaults to `True`, where it defaulted to the value of `caching`: the names
  of the input and of the survey carry the first seven characters of their hash.
- The source file headers say AGPL v3.0 or later (AGPL-3.0-or-later), as the `COPYRIGHT`
  file and the package metadata do.
- The `CI/CD` workflow is now two: `CI` and `CD`. `stable` is moved forward to a commit
  that passed CI by a fast-forward, the actions are pinned to commit SHAs, and the
  workflows only get the permissions they need.
- The install job of CI fetches the C++ backend in its checkout, shallow (only the commit
  this repository pins), instead of leaving the installer to clone it whole, and runs
  Python 3.9 on `ubuntu-22.04`.

### Fixed

- The README said that Python 3.7.12 to 3.11 was supported: the package supports 3.9 to
  3.12.
- The Zenodo record declared the license GPL-3.0-or-later: the package is
  AGPL-3.0-or-later.

## [0.8.0] - 2026-08-10

### Added

- A grouped structure in the output files: the star catalog is in an HDF5 group of its own,
  `StarCatalog`, next to `ParentParticles`, the group of the parent particles.
- The data of the parent particles, their kernels included, are exported with the mock
  star catalogs, with a record of the mass budgeting properties of the parent stars. The
  export takes the caching into account.
- A monkey-patch of vaex that makes it open the HDF5 groups, which the grouped structure
  needs. Importing the package says so with a warning, and silences the log messages of
  vaex.
- Support of astropy 8.

### Changed

- **Breaking:** `Input` takes the initial mass of the particles under a new key,
  `massinit`, in place of the single mass key: `mass` now stands for the present-day mass.
  A script that passes only `mass` fails with an informative error. The key of the initial
  mass in the output catalog is renamed too, in the backend.
- **Breaking:** the output files hold their catalog in a group within the HDF5 structure,
  rather than at the root of the file: whoever reads them directly will not find the former
  flat layout.
- vaex is pinned to one version (4.19), which the monkey-patch requires.

## [0.7.0] - 2026-04-15

### Added

- Support of output surveys that are empty: the pipeline goes through, and the backend
  warns, instead of failing, when it appends magnitudes to an empty star population.
- Metadata in the output files: the parameters of the survey, and the logs of the Galaxia
  runs, consolidated in a compact JSON form. The post-processing flags the EBF files it
  modifies in their logs, and can update the metadata when it flushes.
- Python 3.11 and 3.12 are supported, with numpy 2 and astropy 7, and the installation
  works on Python 3.9 and 3.10: a `pyproject.toml` declares the build requirements.
- A GitHub Actions workflow that publishes to TestPyPI and PyPI. For release candidates it
  waits for the CI run of the tagged commit.

### Changed

- **Breaking:** `Input` receives the sizes of the kernels directly, with the argument
  `kernels`, instead of the densities (`rho_pos` and `rho_vel`) it used to compute them
  from: the densities are no longer accepted, and `make_dummy_kernels_input` replaces
  `make_dummy_densities_input`.
- The `former_kernel` option, which kept the former determination of the smoothing
  lengths, is removed, and `k_factor` is plainly a scale factor applied to the kernel
  lengths.
- The license is the GNU Affero General Public License v3 or later, where it had been the
  GNU General Public License v3 or later, and every source file carries a copyright
  notice. The backend, a modified version of Galaxia, is under the original AGPL-3.0 of
  Galaxia as well: the GPL-3.0 text that had been added to it was superfluous and is
  removed, it carries copyright notices, and files of the original Galaxia were restored
  with updated headers.
- The default spectral unit of the photometric zero points is the jansky.
- The backend no longer uses `unary_function` and `binary_function`, which C++17
  removed: it builds with recent compilers.

### Removed

- Python 3.8 is no longer supported.

## [0.6.1] - 2025-06-28

### Added

- The installation stops with an error on Windows, which is not supported.
- Instructions to install from PyPI, in the README.

### Changed

- The installation no longer calls autotools when it compiles the backend: the generated
  files are touched, so that `make` does not try to regenerate them.
- The git submodule is `Galaxia-ananke`, with a capital, and the backend is updated, with
  its new name in particular.

### Fixed

- The photometric filters are found whatever the case of their names, among the filters
  of the SVO Filter Profile Service.

The release candidate 0.6.1rc0 was yanked from PyPI.

## [0.6.0] - 2025-06-20

### Changed

- **Breaking:** the module has a new name, `galaxia_ananke`, in lower case to comply with
  [PEP 8](https://peps.python.org/pep-0008/): imports have to be updated. The name of the
  package follows, in lower case, ahead of the first upload to PyPI.
- The build was reworked to work with versioneer, and so that the git submodule can be
  cloned on the fly from GitHub, at the commit this version pins, when the sources are not
  a git repository, as in an archive or a built tarball. It takes what was learned with
  enbid-ananke.

## [0.5.0] - 2025-05-28

### Added

- Caching: a survey can reuse the output of a previous run of the same input and
  parameters, which it recognizes by a hash of the input and a hash of the survey
  parameters (`caching`; `append_hash` adds the first characters of the hash to the names).
  It is marked as experimental.
- `partitioning_rule` and `input_sorter`, to choose how the particles and the output are
  split in partitions: the partitions are redefined after Galaxia has run, a partition id
  has to be positive, and a negative one marks stars that are left out of the output.
- Post-processing in a process pool instead of threads, with `max_pp_workers` and
  `pp_auto_flush`, which is on by default: it uses pathos, a new dependency, and
  `classproperty`, a decorator of the package that can be pickled.
- `n_gens`, which can also be an iterable of the numbers of the generations to run, and
  replaces `n_jobs`, now deprecated.
- Isochrone tracks: a `Formatting` class that handles the formats of the tracks, a query
  of the track of a given age and metallicity on `PhotoSystem`, and the data of an
  `Isochrone` consolidated as the inputs of a general photometric interpolator.

### Changed

- The spherical coordinates are computed by the post-processing, where they used to be
  expected properties of the output.
- `Input` copes with `rho_vel` being `None`, and the optional key `dform` is truly
  optional.
- The backend gets a new interpolation scheme, a fix of the HST isochrones that the
  package adapts to, and a more consistent determination of the number of stars.

### Fixed

- The JWST MIRI isochrones were missing one.
- The installation failed on a stray comma.

## [0.4.0] - 2024-10-17

### Added

- The documentation, built with Sphinx and hosted on Read the Docs, with the example
  notebook.
- `n_jobs`, to run several instances of Galaxia in parallel, each generating an
  independent catalog from the same input and parameters, and `verbose`, to control what
  the pipeline prints.
- A `partitionid` column in the input and the output: the catalog is written in one HDF5
  file per partition. The conversion from EBF to HDF5, the post-processing and the
  flushing of the output run in several threads, one per partition.
- A check that the system has what building the backend needs (git, GNU make and g++, with
  their minimum versions), before the build is attempted.
- The metadata of the module hold the full license text, and `.zenodo.json` holds the
  metadata of the Zenodo record.

### Changed

- The version is now derived from the git tags by versioneer.
- The modules `constants`, `defaults` and `templates` are private, so that the
  documentation leaves them out: `_constants`, `_defaults` and `_templates`.
- The querying of the SVO Filter Profile Service is improved, and each parent particle
  kernel has only one star at its center.

## Before 0.4.0

py-Galaxia-ananke was in development long before 0.4.0: the repository dates from June
2022, and the versions below are pre-releases, tagged while the package was taking shape
and not documented as they were made: their GitHub releases hold little more than a link to
the comparison. What follows is a rough history, put together afterwards from the commit
history: the dates are those of the tagged commits, and the descriptions are far from
exhaustive. None of these versions, nor 0.4.0 and 0.5.0, was uploaded to PyPI.

### Alpha versions, 2022-2023

- `v0.0.0-alpha`, 2022-12-14: the first implementation of the wrapper, from June 2022,
  with a `setup.py` that deals with the Galaxia git submodule, a testing notebook (from
  November), the adaptation to the new `photoCateg` parameter of Galaxia, the loading of
  new isochrones on the fly with more control over which isochrone to use (December), and
  the reorganization of the input, survey and output code into the `Input`, `Survey` and
  `Output` classes.
- `v0.1.0-alpha`, 2023-03-21: the output opened as a vaex DataFrame, the location of the
  cache chosen at installation to fit the Python environment, pre-generated Galaxia input
  files accepted by `Input`, the missing right ascension, declination, longitude and
  latitude added to the output, faster conversion from EBF to HDF5, the LSST photometric
  system corrected, a consolidated default color-magnitude box, vaex among the
  dependencies, and `git submodule update --init` forced at installation.
- `v0.1.1-alpha`, 2023-03-21: attempts at fixing an issue with `/dev/tty`, and the version
  number corrected.
- `v0.1.2-alpha`, 2023-04-07: the velocities among the exported keys, more of the output
  exposed as attributes of the vaex DataFrame, and `cmd_magnames` accepted as a
  dictionary.
- `v0.2.0-alpha` and `v0.2.1-alpha`, 2023-04-07, and `v0.2.2-alpha`, 2023-04-12: a change
  of name, from `Galaxia` to `Galaxia_ananke`, and the repository URL that goes with it
  (py-Galaxia became py-Galaxia-ananke), then the renaming of the git submodule (`galaxia`
  became `galaxia_ananke`, then `galaxia-ananke`) and a typo patched in it.
- `v0.2.3-alpha`, 2023-06-04: work on the installation: an error when `install-sh` has no
  permission to execute and the permission set on it, an error when `sys.prefix` cannot be
  written to, a function that compares sets of keys, and looser requirements on the input
  particles.
- `v0.2.4-alpha`, 2023-07-02: the Euclid (Padova) and the JWST (NIRCam and MIRI)
  isochrones, an attempt at Roman ones, `__setitem__` for `Output`, and class methods for
  what py-ananke needs (setting the isochrones, making the export keys).

### Beta versions, 2023-2024

- `v0.3.0-beta`, 2023-11-23: the first GitHub Actions workflow (August 2023), with flake8
  linting and the first tests, run with pytest; typing annotations; the post-processing of
  the output, and a method that converts ICRS coordinates to galactic ones; Gaia DR2
  isochrones, and GAIA as the default photometric system; dummy input generators; a new
  determination of the smoothing lengths, with a scale factor `k_factor` and a flag for the
  former one; checks of the lengths of the particle arrays and a default `parentid`; the
  tags of the file templates held in a singleton; a requirement of Python 3.7.12 to 3.10;
  and better documentation and error catching at installation.
- `v0.3.0-beta2`, 2023-11-30: a minimum version of vaex (4.17) among the dependencies.
- `v0.3.1-beta`, 2024-04-17: the versions of the dependencies constrained, the Python
  requirement aligned with numpy's (3.8 to 3.10), `PhotoSystem` objects in the
  `photometry` module that interface with the isochrones, a new naming convention for the
  photometric systems and their filters (which the release notes warn about) made safe for
  Python names, the SVO Filter Profile Service interface with band names that match it and
  the effective wavelength of each band, the solar position and velocity taken from
  astropy's Galactocentric frame by default, the constants split into templates and
  defaults, and `setuptools` in place of `distutils`.
- `v0.4.0-beta`, 2024-08-16, and `v0.4.0-beta2`, 2024-08-21: the first Sphinx documentation,
  on Read the Docs.
- `v0.4.0b3`, 2024-09-01: versioneer added, and the documentation formatting worked on.

Nothing was tagged between `v0.3.1-beta` (April 2024) and `v0.4.0-beta` (August 2024), but
the work went on: the commits of those months add the partition id to the input and the
output, run several instances of Galaxia in parallel, convert and post-process in several
threads, improve the metadata, and prepare the documentation, which 0.4.0 describes.

[Unreleased]: https://github.com/athob/py-Galaxia-ananke/compare/v0.8.0...develop
[0.8.0]: https://github.com/athob/py-Galaxia-ananke/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/athob/py-Galaxia-ananke/compare/v0.6.1...v0.7.0
[0.6.1]: https://github.com/athob/py-Galaxia-ananke/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/athob/py-Galaxia-ananke/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/athob/py-Galaxia-ananke/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/athob/py-Galaxia-ananke/compare/v0.3.1-beta...v0.4.0
