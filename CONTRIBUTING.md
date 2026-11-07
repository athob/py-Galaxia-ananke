# Contributing to py-Galaxia-ananke

Thanks for helping out! This page is the short version of "how work gets into
py-Galaxia-ananke": the branches, how to set up a development copy, what the
automatic checks do, and how versions and releases work. Everyone taking part is
expected to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

py-Galaxia-ananke is the wrapper [py-ananke](https://github.com/athob/py-ananke)
uses to run Galaxia, so a change to a function's, a class's or a method's signature
or behaviour can break py-ananke. If you plan something like that, open an issue
first so it can be discussed.

## The branches

| Branch    | What it is                                                         | Who writes to it                              |
|-----------|--------------------------------------------------------------------|-----------------------------------------------|
| `main`    | The latest **release**.                                            | The maintainer, when releasing.               |
| `stable`  | The newest development state that passed the automatic checks.     | **Only CI** -- never commit to it by hand.    |
| `develop` | Where development happens; every other branch starts from it.      | The maintainer, and merged pull requests.     |

Everything else (`my-feature`, `fix-something`, ...) is a short-lived branch
created from `develop`. When the checks pass on a push to `develop`, CI moves
`stable` forward to that commit by itself.

## Making a change

1. Branch from `develop`: `git switch develop && git pull && git switch -c my-change`.
2. Make the change, with a test for it in `src/tests`.
3. Run the checks locally (next section).
4. Open a pull request into `develop`. Keep a pull request about one thing.

## Setting up

```bash
git clone https://github.com/athob/py-Galaxia-ananke
cd py-Galaxia-ananke
git switch develop
python -m venv .venv && source .venv/bin/activate      # or a conda environment
pip install .
pip install pytest ruff
```

`pip install .` compiles the C++ backend, [Galaxia-ananke](https://github.com/athob/Galaxia-ananke),
and installs it with the package. It is a git submodule of this repository, and it is large:
about 5.8 GB of files, most of them data, and about a minute to build. So installing needs
network access, git, GNU make and g++ (the versions are in the README), several GB of free
disk space, and Windows is not supported. The tests import the package straight from `src/`
(they use relative imports), so once the install has worked your edits are picked up without
installing again. The version is read from the git tags (see "Versioning" below), so work from
a clone that has them: `git fetch --tags` if you cloned a fork.

The URL of the submodule in `.gitmodules` is relative, `../../athob/Galaxia-ananke.git`, with
the owner spelled out: git resolves it against the URL of the repository you cloned, so a
fork gets the same backend, and CI runs in a fork. Keep it that way, because the installer
reads it too, to clone the backend when the sources are not a git repository, which is what
an sdist is.

A change to the C++ code itself is made in the Galaxia-ananke repository, which has its own
issues and pull requests. This repository then pins the new commit of the submodule, and its
changelog tells the change under the version that pins it.

Before you push, run what CI runs:

```bash
ruff check .      # must pass
pytest src
```

To have ruff, and actionlint if you have it, run on every commit by itself, install the optional
[pre-commit](https://pre-commit.com/) hook: `pip install pre-commit && pre-commit install`.

## How the package is built

`setup.py` only describes the package: reading it builds nothing. What builds the backend is in
`src/_build_utils.py` (`make_cmdclass`): a `build_ext` command, which setuptools runs when the
wheel is built, checks that the system has git, GNU make and g++, makes sure the submodule is
there (with `git submodule update`, and, when the sources are not a git repository, by cloning
the backend at the commit this version pins, which it asks the GitHub API for), then runs
`configure`, `make` and `make install` for it, and copies its data. All of it goes to a cache
directory, `.cache/galaxia_ananke` under `sys.prefix` (or under the folder
`ANANKE_SYSTEM_PREFIX` names), with the log of each step in its `log` folder. Installing again
rebuilds that directory, and keeps the custom photometric catalogue (`py_custom`) that the
user made in it.

The backend is therefore not in the package, and not in the sdist either: the sdist holds the
Python sources, `.gitmodules`, and the version that versioneer wrote in it. That is also why only
an sdist is published, and why a release needs the commit of the backend it pins to be on
GitHub: [`RELEASING.md`](RELEASING.md) says it all.

## The automatic checks (CI)

The **CI** workflow ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) runs on
every push to a branch other than `main` and `stable`, on pull requests from forks
into `develop`, every Monday (to notice new releases of the dependencies, or trouble
fetching the backend, before anyone else does) and on demand from the Actions tab. A run
takes about four minutes. Its jobs:

- **lint**: [ruff](https://docs.astral.sh/ruff/). Only the serious findings fail it
  (syntax errors and undefined names, selected in `pyproject.toml`); the style
  warnings are listed for information.
- **workflows**: [actionlint](https://github.com/rhysd/actionlint) on the workflow
  files. You can run it locally too, and the pre-commit hook does so when it is
  installed.
- **install**: for each Python version from 3.9 to 3.12, fetches the backend (shallow: the
  one commit this repository pins), installs the package and runs `pytest src`.
- **package**: builds the sdist the way a release would, then runs `twine check --strict`
  and `.github/scripts/check_dist.py`: the sdist has to hold what installing needs
  (`.gitmodules` among it), with a version and a revision id that the installer can use, and
  not the backend.
- **promote**: on a push to `develop`, once all the jobs above have passed, moves `stable`
  forward to that commit.

Apart from CI, two workflows report to Security > Code scanning: **CodeQL** (GitHub's code
scanning of the Python source; on `develop`, on pull requests and every Wednesday) and
**Scorecard** (the [OpenSSF](https://securityscorecards.dev) grading of the repository's
supply-chain habits; from `main`, every Thursday). Their findings are meant to be looked at,
and `stable` does not wait for them.

The tools the workflows install (ruff, pytest, build and twine, and what the CD workflow builds with) are
locked with hashes in
[`.github/requirements/`](.github/requirements): each `requirements-X.in` says what is
wanted, and `requirements-X.txt` pins every package, transitive ones included, to a version
and to the hashes of its files, which `pip install --require-hashes` then insists on.
Dependabot is set up to propose updates to the locks, by running again the `uv pip compile`
command a lock records in its header (which is why `dependabot.yml` says `uv`, not `pip`:
pip-tools would resolve for one Python only). As in py-EnBiD-ananke, GitHub may refuse to
open those pull requests (HTTP 400, see
[its issue #12](https://github.com/athob/py-EnBiD-ananke/issues/12)), and until that is
understood the locks are refreshed by hand, as described below. Dependabot takes one
version of each package as the current one, and checks that the newest can be pinned for
every Python the lock covers, so it cannot maintain a lock that holds a package once per
Python. The test tools are therefore locked twice. `requirements-test.txt` covers the Pythons
that are still supported (3.11 and up), and is the one Dependabot is set up to maintain.
`requirements-frozen/requirements-test.txt` is a complete lock of every Python from the
oldest the package supports (3.9) up, used by the ones at the end of their life (3.9 and
3.10). It is next to `requirements/` and not in it, because Dependabot fetches the
subfolders of its folder too. Because it is complete, retiring a Python to it never touches
it: raise the `--python-version` of the test lock in `update.sh`, add that Python to the
matrix of `ci.yml` (the tests check that both agree), and run `update.sh`. Never copy the
live lock into `requirements-frozen/`, which would lose the entries of the older Pythons (the
tests catch that too). `update.sh --upgrade` refreshes the frozen lock when that is worth it,
and uv gives each Python the newest version it supports.
To refresh the locks, or to change one by hand (edit its `.in` file first), run
`bash .github/requirements/update.sh` (add `--upgrade` to move everything to the newest
versions; it needs [uv](https://docs.astral.sh/uv/)), then commit the `.txt` files it wrote.
If a test then says that a lock pins a package more than once, that package has a newer
version that only newer Pythons can take: bound it in the `.in` file, with a comment that says
until when, because the ignores of `dependabot.yml` are not read by `update.sh`.
The same script writes the checksum the actionlint download is checked against: run it
again whenever `ACTIONLINT_VERSION` changes in `ci.yml`.

Dependabot raises the security advisories of the packages in a lock as alerts of its own,
and each new advisory on an old pin does the same. The frozen lock holds the pins of the
Pythons that are past their end of life, and a fix often needs a newer Python than they
are, so it cannot always take it. Before dismissing an alert, check that the pin for the
Pythons that are still supported, in `requirements-test.txt` and in the entries of the
frozen lock for them, is a patched one: if it is not, the alert is real, and the lock needs
updating. Dismiss the others as a tolerable risk, so that Dependabot does not try to open a
security update for a lock that cannot take one (the number is in the alert's URL):

```bash
gh api -X PATCH repos/athob/py-Galaxia-ananke/dependabot/alerts/NUMBER -f state=dismissed \
  -f dismissed_reason=tolerable_risk \
  -f dismissed_comment="Test-time only. The fix needs a newer Python than the lock's end-of-life ones: see CONTRIBUTING.md."
```

Three alerts have been dismissed this way, all against the pins that the frozen lock holds for
Python 3.9: the releases that fix them need Python 3.10 or later, and the lock pins patched
versions for 3.10 and later, as `requirements-test.txt` does.

- **pytest**: [GHSA-6w46-j5rx-g56g](https://github.com/advisories/GHSA-6w46-j5rx-g56g)
  (CVE-2025-71176, [PYSEC-2026-1845](https://osv.dev/PYSEC-2026-1845)) affects every pytest
  before 9.0.3. It only matters when the tests run on a machine shared with other users (it
  is about pytest's temporary directories in `/tmp`), and it goes away with the support of
  Python 3.9.
- **setuptools**: [GHSA-h35f-9h28-mq5c](https://github.com/advisories/GHSA-h35f-9h28-mq5c)
  (CVE-2026-59890: on macOS, the exclusions of a `MANIFEST.in` can be bypassed when an sdist
  is built; fixed in 83.0.0) and
  [GHSA-xmj3-9gf7-h52w](https://github.com/advisories/GHSA-xmj3-9gf7-h52w) (on Windows, a path
  can escape while a legacy wheel is installed; fixed in 84.1.0). The tests that use this lock
  run on Linux only, and the release is built on Python 3.12 (`cd.yml`), where the build
  environment resolves a setuptools that has the fixes.

If CI fails on one Python version only, it is usually a real incompatibility
rather than bad luck.

## Conventions

- **Supported Pythons** are the ones in `python_requires` in `setup.py`
  (currently 3.9 up to 3.12). CI covers all of them, and the classifiers in
  `src/galaxia_ananke/__metadata__.py` name them for PyPI: the tests check that the three agree.
  The ceiling is the one of vaex.
- **vaex is pinned to one version** in `setup.py`, because `src/galaxia_ananke/_patch_vaex.py`
  patches the way vaex reads HDF5 files: moving to another version of vaex means checking that
  the patch still works.
- **Do not edit `versioneer.py` or `src/galaxia_ananke/_version.py` by hand.**
  They are generated by [versioneer](https://github.com/python-versioneer/python-versioneer),
  configured in `setup.cfg`.
- **New behaviour comes with a test** in `src/tests`. A survey run takes under a minute when its
  input is small and its colour-magnitude box (`cmd_box`) cuts the faint stars, which are the
  numerous ones: keep the tests that run the real backend to that.
- **The repository checks itself**: `src/tests/test_repository.py` verifies, among other
  things, that the DOIs and the license wording agree between the README, `CITATION.cff`,
  `.zenodo.json` and the metadata, that the actions are pinned and the workflows set their
  permissions, that what they install is locked with hashes, that the submodule is declared
  the way the installer reads it, and that the changelog has its format. If you change one of
  those, it tells you what else needs to follow.
- **A change a user would notice gets a line in `CHANGELOG.md`**, under
  `[Unreleased]`, in the matching group (Added, Changed, Deprecated, Removed,
  Fixed or Security). The release moves those lines under the new version.

## Versioning

The package version is not written down anywhere in the source: versioneer
derives it from the git tags, which start with `v` (`v0.8.0`). Version numbers
follow [Semantic Versioning](https://semver.org/). While the major number is 0, up to
`1.0.0`, a minor version may change the interface, as 0.6.0, 0.7.0 and 0.8.0 did; from `1.0.0`
on, a change that breaks it (a function, a class, a method or a property of the documentation,
its arguments or its results) is a new major version, and the changelog says so. Pre-releases
are written the Python way ([PEP 440](https://peps.python.org/pep-0440/)): `v0.9.0a1`,
`v0.9.0b1`, `v0.9.0rc0`, `v0.9.0.dev1`.

## Releasing (maintainer)

Pushing a tag that starts with `v` runs the **CD** workflow
([`.github/workflows/cd.yml`](.github/workflows/cd.yml)), which publishes to TestPyPI
and PyPI. What the tags may look like, and the steps, are in
[`RELEASING.md`](RELEASING.md).

## Questions and bug reports

Security problems are the exception: report them privately, as
[`SECURITY.md`](SECURITY.md) describes. For anything else, open an issue; the
templates ask for what is needed to reproduce it.

By contributing you agree that your contribution is licensed under the
project's license, the GNU Affero General Public License, version 3 or (at your option)
any later version (see `LICENSE` and `COPYRIGHT`).
