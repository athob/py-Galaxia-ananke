# Releasing

How a version of py-Galaxia-ananke gets released. It is a short checklist: the
[CD workflow](.github/workflows/cd.yml) does the publishing, and starts when a tag is
pushed.

## Which tag

A tag is `v` followed by a plain [PEP 440](https://peps.python.org/pep-0440/) version.
What happens depends on its kind:

| Tag                                   | Kind        | Published to                                     |
|---------------------------------------|-------------|--------------------------------------------------|
| `v0.9.0`, `v0.9.0.post1`              | release     | TestPyPI, then PyPI, then a GitHub release       |
| `v0.9.0a1`, `v0.9.0b1`, `v0.9.0rc0`   | pre-release | TestPyPI only, once CI has passed on the commit  |
| `v0.9.0.dev0`, `v0.9.0rc0.dev1`       | development | TestPyPI only, straight away                     |

Anything else (`v0.9.0-test`, `v0.9.0dev1`, ...) is refused by the first job of the
workflow, so that nothing that is not a clean release can reach PyPI. Development tags
are for trying the pipeline out. TestPyPI refuses a file name it has already seen, even
after deleting it, so every try needs a new number (`.dev1`, `.dev2`, ...).

## What is published

py-Galaxia-ananke is published as an sdist, and as nothing else. Installing it compiles
the C++ backend, Galaxia-ananke, on the machine of the user, and puts it in a cache
directory under `sys.prefix`, outside the package. A wheel would therefore not carry the
backend and could not work, and a wheel specific to a platform is refused by PyPI unless it
is a manylinux one. So the build is only Python files, and takes seconds.

The sdist holds the Python sources, `.gitmodules`, and the `_version.py` in which
versioneer wrote the version and the full revision id of the tagged commit. It does not
hold the backend, which is the git submodule `Galaxia-ananke` of this repository, about
5.8 GB of files. When pip builds the package, the installer first asks git for the
submodule. An sdist is not a git repository, so that fails, and the installer then asks
the GitHub API which commit of Galaxia-ananke the tagged revision pins, clones the
repository and checks that commit out.

Two things follow for a release:

- The tagged commit has to be on GitHub, because the installer asks GitHub about that
  revision. Pushing the tag does it.
- The commit of Galaxia-ananke that the release pins has to be on GitHub too, in
  `athob/Galaxia-ananke`: push it before the tag. `git submodule status` shows which commit
  it is.

## Trying the build

To try a change to the build without publishing anything, run the CD workflow by hand:
from the Actions tab, or with `gh workflow run cd.yml --ref BRANCH`. It builds the sdist,
lists it in the summary of the run and keeps it as the `python-package-distributions`
artifact; it publishes nothing, on any branch or tag. GitHub only offers this for a
workflow whose trigger is on the default branch, so `main` has to have it first.

## Steps

1. Everything for the release is on `develop`, and the CI run of that commit is green
   (`stable` then points at it). The commit of Galaxia-ananke it pins is on GitHub.
2. In [`CHANGELOG.md`](CHANGELOG.md), turn `## [Unreleased]` into
   `## [0.9.0] - YYYY-MM-DD`, add a new empty `## [Unreleased]` above it, and update the
   links at the bottom of the file (`[Unreleased]` compares the new tag with `develop`;
   `[0.9.0]` compares the previous tag with the new one). Commit it.
3. Optionally try the pipeline first with a pre-release, `v0.9.0rc0`: it goes to TestPyPI
   after CI, and the release itself then only differs by going to PyPI. To see that the
   sdist installs, do it in a clean environment, on a machine that has git, g++ and make,
   network access and the disk space the backend needs:

   ```bash
   pip install --index-url https://test.pypi.org/simple/ \
               --extra-index-url https://pypi.org/simple/ galaxia-ananke==0.9.0rc0
   ```

4. Tag the commit and push the tag:

   ```bash
   git tag v0.9.0
   git push origin v0.9.0
   ```

5. Move `main` to the released commit: `git switch main && git merge --ff-only v0.9.0
   && git push`.

The workflow then runs these jobs, and stops at the first one that fails:

1. **Classify the tag**: the tag has to be a version as above, and for a release
   `CHANGELOG.md` has to have its dated section. This is the reason to write the
   changelog first: nothing is built or published without it.
2. **Build**: for a pre-release, waits for the CI run of the tagged commit to succeed;
   then builds the sdist on the runner.
3. **TestPyPI**, then, for a release only, **PyPI**.
4. **GitHub release**, for a release only: its notes are the changelog section, and the
   sdist is attached. Zenodo archives the published release and gives that version its
   DOI.

Publishing uses PyPI's trusted publishing, through the `release` environment: no API
token is stored anywhere. The trusted publisher registered on both PyPI and TestPyPI
names the workflow file, `cd.yml`, and that environment.

## Afterwards

- Check the new version on [PyPI](https://pypi.org/project/galaxia-ananke/), the
  [release](https://github.com/athob/py-Galaxia-ananke/releases) and its
  [Zenodo record](https://doi.org/10.5281/zenodo.13942945) (that DOI resolves to the
  latest version), and that the documentation built on
  [Read the Docs](https://py-galaxia-ananke.readthedocs.io).
- Update what depends on the new version (py-ananke, for one).

## If it goes wrong

- A tag refused by the first job, or a job that failed before anything was published:
  fix the cause, delete the tag (`git tag -d v0.9.0 && git push --delete origin v0.9.0`)
  and tag again.
- A job that failed after the file reached TestPyPI or PyPI: do not tag again, as those
  services never accept the same file twice. Use **Re-run failed jobs** on the run in the
  Actions tab, which only repeats what did not succeed. If the file itself was wrong,
  release a new version instead (`v0.9.0.post1`, or the next one).
- A GitHub release created by hand before the workflow got there is completed by it (its
  notes and files are replaced), not refused.
- An installation that fails with `GitHub API rate limit exceeded`: the installer asks the
  API without a token, which allows few requests per hour for an address, so a shared
  network can run out. Waiting is the remedy.
