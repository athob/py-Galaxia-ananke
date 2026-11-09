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
- A check of the workflow files in CI, with actionlint, whose download is verified against a
  checksum.
- A README that gives the installation requirements, the size of the C++ backend among
  them, and says how to cite the package.
- A script, `.github/scripts/release_meta.py`, that the CD workflow runs on the pushed tag,
  and tests of it.
- A job of CI, `package`, that builds the sdist and checks it with `twine check --strict` and
  `check_dist.py`, which `stable` waits for like for the others.
- A script, `.github/scripts/check_dist.py`, that checks the sdist before it is published: it
  has to hold what installing needs, which is the Python sources, `.gitmodules` and the version
  that versioneer baked in, with a revision id the installer can look up, and it must not hold
  the backend.
- The CD workflow can be run by hand: it builds the sdist and keeps it as an artifact,
  and publishes nothing.
- The CD workflow creates the GitHub release of each release once it is on PyPI, with its
  section of this changelog as the notes and the sdist attached. A release made by hand
  beforehand is completed, not refused.
- A releasing guide, `RELEASING.md`, which also says why the package is published as an sdist
  only, and what the pinned backend commit has to satisfy for an sdist to install.
- The requirements of the tools that CI and CD install, in `.github/requirements/`, with a
  script, `update.sh`, that locks them with hashes (`uv pip compile`) and writes the checksum
  of the actionlint download.
- Tests of the repository itself, in `src/tests/test_repository.py`, that keep the project
  files, the workflows and their habits, the locks of the CI requirements, the supported
  Pythons, the license and the citation files, the git submodule, the changelog and the
  contributing guide in agreement. They do not need Galaxia.
- Dependabot also watches the locked requirements of the CI and CD tools, with its `uv`
  ecosystem.
- The test coverage is measured by branch too, without the footer of the modules that no test
  can run, and it has a floor, 95% for now, which the Python 3.12 run of CI fails under.
- Tests of the utilities that the package is built on: the state that tracks the steps of a
  survey and lets a run with caching skip the ones already done, `execute`, the checks of the
  dictionaries of particles and of the lengths of their arrays, the symbolic links, the version
  of a command, the loading of a source file, and the metadata of the classes and the compact
  form of the logs of the backend. The known bug of `execute`, which cuts a command on its
  spaces, is an expected failure that fails the run on the day it is fixed.
- Tests of the templates that write the parameter file and the commands of the backend, with
  the parameter file in the format that the backend reads; of the constants, which name what the
  files of a survey are written with and where the cache is; of the metadata of the package,
  including its date when versioneer has none; of the defaults, including the position and the
  velocity of the Sun that astropy gives; of the hash that names the cache of an input; and of
  what importing the package does to vaex. The known bug of a prefix that cannot be written to,
  which is not refused, is an expected failure.
- Tests of the tools that build the package: what the system must have, the lookup on GitHub
  of the commit that a version pins, the cache of an earlier installation with the custom
  photometric systems in it, the build of a toy backend with the real make, and the fetch of the
  backend with the real git on repositories laid out as on GitHub, including the sources that are
  not a git repository, as an sdist unpacks them. The weaknesses that they found are expected
  failures: a step of the build, a clone, a checkout or a git that fails goes unnoticed, the
  backend is built from a name relative to the working directory, and the custom systems cannot be
  moved between two file systems.
- Tests of `Input`: how it is made, what it refuses and what it completes, its hash (pinned,
  since it names the cache of an input), the names of the files that follow from it, the files
  that it writes and how it prepares them for the backend. The weaknesses that they found are
  expected failures: the explanation of the change of the mass keys fails on an import, a wrong
  number of arguments reaches a bare `raise`, the default number of neighbours is the name of its
  tag, an input that was just made cannot say its hash, a name with a dot loses its hash in the
  names of its files, the kernels and the particles cannot be given as lists, and the kernels are
  not checked to be as many as the particles.
- Tests of `Survey`, without running the backend: its parameters and the colors of its diagram,
  the photometric systems after the first, the generations, the name that follows from its hash,
  the parameter file, the commands that run the backend and append a photometric system, running
  with and without caching, and what making a survey does in what order. The weaknesses that they
  found are expected failures: a survey that is made again with other parameters keeps the name
  of the first, a parameter that does not exist is not refused, no generation fails with a message
  that says nothing of the generations, and `n_gens` and `fsample` of a survey that was not made
  fail on `None`.
- Tests of the files and layouts of the isochrones (`IsochroneFile`, `Formatting`, `Isochrone`):
  the name of the files and the header that is read, the four layouts that are known, the tables
  of a star for each age with the radius that is made from its mass and gravity, the descriptor
  of an isochrone made in Python, and the points and values that are made of the tables. The
  weaknesses that they found are expected failures: a header with white space at its ends (a tab
  first, spaces or a tab last) has an empty name, an isochrone that is refused leaves its
  directory behind and cannot be made again under its name, a system with a single magnitude gives
  its column and not a tuple of columns, and an oldest age that does not end with a white dwarf
  fails on an index.
- Tests of the photometric systems, without the network: the registry that is read from the
  installation, the dictionary with slashes that it is made of, the isochrones that are added
  and the clearing of them, what a system says about itself and the check of the colors of a
  diagram, the instruments that a system is asked about, and the effective wavelengths and zero
  points from the saved answers of the SVO filter profile service (Gaia, and LSST). The
  weaknesses that they found are expected failures: the protection of the padova category compares
  by identity, a name with a slash is set in a copy, an isochrone that is added is not a
  photometric system, the filters of LSST are eight for six magnitudes, the tables of an
  isochrone that is made in Python fail on duplicate column names, and the particles that share
  their nearest metallicity fail on an index.
- Tests of what `Output` is made of that needs no survey: the longitude within bounds, the
  partitions shared between processes, the keys of a catalogue, the conversions of the coordinates
  (checked against physics: a star at 1 kpc that moves across the sky at 4.74 km/s moves by 1
  mas/yr; the center and the pole of the Galaxy; formulas that do not use astropy; a round trip),
  the helpers that flush what a pipeline made into the HDF5 file of a catalogue, and the patch
  that the package makes to vaex to open the groups of its files and lists of files.
- Tests of `Output` on synthetic output of the backend (EBF files laid out as Galaxia writes
  them): what it says of itself, the files of the backend that it looks for and clears, the stars
  of each partition and the runs of them that are copied, the names of the HDF5 files, what those
  hold (the stars, the groups, the attributes, the logs of the backend), the redefinition of the
  partitions by a rule of the user, and the parser of the logs of the backend.
- More tests of `Output` on synthetic output of the backend: the post-processing of the
  catalogue end to end (the new columns checked against the pipelines applied to the table of
  the stars, what is flushed into the files, what is said when it runs), the group of the parent
  particles with the statistics of their stars worked out by hand, the use of an output like the
  vaex table that it holds, the pipelines of the user run on the file of each partition with
  their options, and the flushing of the columns of the loaded catalogue. The length of an
  output, which is not the number of its stars, and the pipelines of a catalogue whose partitions
  are not flushed by themselves, which are given the numbers of the partitions and not their
  tables, are expected failures.
- Tests of the two functions that the processes of `Output` run, called in the process of the
  tests (the measure of the coverage does not follow a pool of processes): the HDF5 file of a
  partition, checked against the one that the processes wrote, what it says when asked to, the file
  of the stars that are in no partition, and the redefinition of the partitions with a rule. A
  catalogue without any star is post-processed in a probe.
- Integration tests that run a small survey with the compiled backend (skipped where it is not
  built), once for all the tests: the files, the keys and the types of the catalogue, the limits
  that were asked for, the relations between the columns, the consistency of what is said of the
  parent particles, and that the files of the backend are laid out as the synthetic ones that the
  other tests of `Output` use.
- Tests of `make_survey_from_particles`: it runs a survey end to end, the same input or the same
  files make the same survey, and a call with the particles and no kernels, which fails on a
  variable that was not made instead of saying what the calls are, is an expected failure.
- The classifiers of the package list the Python versions it supports, 3.9 to 3.12, so that
  PyPI shows them.
- CodeQL code scanning of the Python source, on `develop`, on its pull requests and every
  week, which reports to Security > Code scanning and which `stable` does not wait for.
- OpenSSF Scorecard, which grades the supply-chain habits of the repository every week from
  `main`, reports to Security > Code scanning and which `stable` does not wait for.
- Dependabot also watches the pinned commit of the C++ backend, with its `gitsubmodule`
  ecosystem: it proposes the newer ones as pull requests that change the pin alone.
- A contributing guide, `CONTRIBUTING.md`: the branches, how to set up a copy with the C++
  backend that is a git submodule, how the package is built, what CI does, and how versions
  and releases work.

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
- The CD workflow refuses a tag that is not a plain PEP 440 version, and for a release it
  checks that `CHANGELOG.md` has the dated section of that version before anything is
  built or published.
- The CD workflow builds the sdist directly on the runner, with Python 3.12, instead of
  inside a manylinux2014 container that a package which compiles nothing when it is built
  has no use for, and lists what it built in the summary of the run.
- The install job of CI fetches the C++ backend in its checkout, shallow (only the commit
  this repository pins), instead of leaving the installer to clone it whole, and runs
  Python 3.9 on `ubuntu-22.04`.
- CI installs the tools it uses, ruff and pytest, from the hash-locked requirements: pip
  refuses anything that is not what the lock says. The Pythons past their end of life, 3.9
  and 3.10, use a lock of their own for the test tools. So does the CD workflow with the
  tools that build the sdist.

### Fixed

- The metadata of the package states its license in a `License` field, where only a
  classifier did.
- Building the package from files that have no git information, such as an archive of a commit
  that no tag names or a copy of the files made without `.git`, failed for want of a date, and so
  did importing a package built from them. The date is that of now when versioneer has none, and
  the version `0+unknown`: releases, which have both, are not affected.
- Two regular expressions and a docstring held a backslash that Python reads as an invalid
  escape sequence, and warns about (a `SyntaxWarning` on Python 3.12, a
  `DeprecationWarning` before): they are raw strings now, with the same value.
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

- **Breaking:** `Input` takes the mass of the particles under a new key, `massinit`, which is
  their initial mass, in place of the key `mass`, which it no longer accepts. A script that
  passes only `mass` fails: the error meant to explain the change fails itself, on a wrong
  import. The key of the initial mass in the output catalog is renamed too, in the backend.
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
