#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2026  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later).
#
# The full copyright notice, including terms governing use, modification,
# and redistribution, is contained in the files LICENSE and COPYRIGHT,
# which can be found at the root of the source code distribution tree:
# - LICENSE <https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE>
# - COPYRIGHT <https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT>
#
"""
Tests of the repository itself rather than of the code: the project files are there, the
workflows keep the habits they were given, and what has to agree with something else does
(the locks of the CI requirements and the workflows that install them, the matrix and the
lock each of its Pythons uses). Standard library only, and nothing here needs Galaxia.
"""
import configparser
import json
import pathlib
import re
import warnings

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
WORKFLOWS = sorted((ROOT / ".github" / "workflows").glob("*.yml"))

REQUIREMENTS = ROOT / ".github" / "requirements"
FROZEN = ROOT / ".github" / "requirements-frozen"
LOCKS = sorted(REQUIREMENTS.glob("requirements-*.txt")) + sorted(FROZEN.glob("requirements-*.txt"))
SOURCES = sorted(REQUIREMENTS.glob("requirements-*.in")) + sorted(FROZEN.glob("requirements-*.in"))


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def prose(path):
    """A document as running text: the way its lines break is not part of what it says."""
    return " ".join(read(path).split())


def header_of(lock):
    """The first lines of a lock, where uv records the command that wrote it."""
    return "\n".join(lock.read_text(encoding="utf-8").splitlines()[:3])


# -- the project files and the workflows ------------------------------------------

@pytest.mark.parametrize("path", [
    "README.md", "LICENSE", "COPYRIGHT", "CHANGELOG.md", "CONTRIBUTING.md", "RELEASING.md", "CODE_OF_CONDUCT.md",
    "SECURITY.md", "CITATION.cff", ".zenodo.json", ".editorconfig", ".gitattributes", ".gitmodules",
    ".pre-commit-config.yaml", ".markdownlint.yaml", ".readthedocs.yaml",
    ".github/dependabot.yml", ".github/PULL_REQUEST_TEMPLATE.md",
    ".github/ISSUE_TEMPLATE/bug_report.yml", ".github/ISSUE_TEMPLATE/feature_request.yml",
    ".github/ISSUE_TEMPLATE/config.yml",
])
def test_project_file_exists(path):
    assert (ROOT / path).is_file(), f"{path} is missing"


def test_gitattributes_keeps_versioneer_export_subst():
    # versioneer needs this to find the version in a `git archive`
    assert "src/galaxia_ananke/_version.py export-subst" in read(".gitattributes")


def test_there_are_workflows():
    assert {p.name for p in WORKFLOWS} >= {"ci.yml", "cd.yml"}


@pytest.mark.parametrize("workflow", WORKFLOWS, ids=lambda p: p.name)
def test_actions_are_pinned_to_a_commit(workflow):
    unpinned = []
    for line in workflow.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^\s*(?:-\s+)?uses:\s+(?P<action>[^@\s]+)@(?P<ref>\S+)(?P<rest>.*)$", line)
        if match and not match.group("action").startswith("./"):
            if not re.fullmatch(r"[0-9a-f]{40}", match.group("ref")) or "#" not in match.group("rest"):
                unpinned.append(line.strip())
    assert not unpinned, f"{workflow.name}: pin to a commit SHA, with the version in a comment: {unpinned}"


@pytest.mark.parametrize("workflow", WORKFLOWS, ids=lambda p: p.name)
def test_workflows_set_their_permissions(workflow):
    text = workflow.read_text(encoding="utf-8")
    assert re.search(r"^permissions:", text, flags=re.MULTILINE), \
        f"{workflow.name} has no top-level `permissions:`: the token would get the repository default"


def test_promote_waits_for_every_other_job_of_ci():
    # a job added to ci.yml and forgotten in promote's needs would not hold stable back
    ci = read(".github/workflows/ci.yml")
    jobs = re.findall(r"^  ([a-z][\w-]*):\s*$", ci.split("\njobs:\n", 1)[1], flags=re.MULTILINE)
    needs = re.search(r"^  promote:.*?^    needs:\s*\[([^\]]*)\]", ci, flags=re.MULTILINE | re.DOTALL).group(1)
    assert set(jobs) - {"promote"} == {name.strip() for name in needs.split(",")}


def test_cd_publishes_only_on_a_push():
    # the workflow can also be run by hand, to try a build out: that must never publish anything
    cd = read(".github/workflows/cd.yml")
    jobs = dict(re.findall(r"^  ([a-z][\w-]*):\s*\n(.*?)(?=^  [a-z][\w-]*:\s*\n|\Z)",
                           cd.split("\njobs:\n", 1)[1], flags=re.MULTILINE | re.DOTALL))
    for name in ("publish-to-testpypi", "publish-to-pypi", "github-release"):
        condition = re.search(r"^    if:(.*)$", jobs[name], flags=re.MULTILINE)
        assert condition and "github.event_name == 'push'" in condition.group(1), \
            f"{name} has to run only on a push: a run by hand must not publish"


def test_no_job_of_ci_is_allowed_to_fail():
    # a job that may fail is red without holding stable back: if one is ever meant to, say so here
    assert "continue-on-error" not in read(".github/workflows/ci.yml")


def test_cd_waits_for_a_workflow_that_exists():
    match = re.search(r'WORKFLOW_FILE="([^"]+)"', read(".github/workflows/cd.yml"))
    assert match, "the CD workflow no longer says which workflow it waits for"
    assert (ROOT / ".github" / "workflows" / match.group(1)).is_file()


# -- the submodule, and what CI and CD do with it ---------------------------------

FORK_GUARD = "github.event_name != 'pull_request' || github.event.pull_request.head.repo.full_name != github.repository"


def jobs_of(workflow):
    """The jobs of a workflow, by name, as their text."""
    text = read(f".github/workflows/{workflow}")
    return dict(re.findall(r"^  ([a-z][\w-]*):\s*\n(.*?)(?=^  [a-z][\w-]*:\s*\n|\Z)",
                           text.split("\njobs:\n", 1)[1], flags=re.MULTILINE | re.DOTALL))


def needs_of(job):
    """The jobs a job needs, written as `needs: job` or as a list."""
    return re.findall(r"^    - (\S+)\s*$", job, flags=re.MULTILINE) + re.findall(r"^    needs:\s*(\w[\w-]*)\s*$", job, flags=re.MULTILINE)


def test_the_submodule_is_declared_the_way_the_installer_reads_it():
    # When the sources are not a git repository (an sdist), the installer reads .gitmodules itself: it
    # joins the url, which has to be relative to the one of this repository, to the url of this repository.
    parser = configparser.ConfigParser()
    parser.read(ROOT / ".gitmodules", encoding="utf-8")
    name = re.search(r"^GALAXIA_SUBMODULE_NAME\s*=\s*['\"]([^'\"]+)['\"]", read("src/galaxia_ananke/_constants.py"),
                     flags=re.MULTILINE).group(1)
    assert parser.sections() == [f'submodule "{name}"'], parser.sections()
    section = parser[parser.sections()[0]]
    assert section["path"] == name
    assert re.fullmatch(r"\.\./\.\./[\w.\-]+/Galaxia-ananke\.git", section["url"]), \
        f"the url of the submodule has to stay relative, like ../../owner/Galaxia-ananke.git: {section['url']}"


def test_the_manifest_ships_the_declaration_of_the_submodule_and_not_the_submodule():
    manifest = read("MANIFEST.in").splitlines()
    assert "include .gitmodules" in manifest, "an sdist needs .gitmodules: the installer reads the url of the backend in it"
    name = re.search(r"^GALAXIA_SUBMODULE_NAME\s*=\s*['\"]([^'\"]+)['\"]", read("src/galaxia_ananke/_constants.py"),
                     flags=re.MULTILINE).group(1)
    assert not [line for line in manifest if name in line and not line.startswith("prune")], \
        "the backend is fetched when the package is installed: it is not part of the sdist"


def test_only_the_jobs_of_ci_that_install_the_package_fetch_the_submodule():
    # the commit the repository pins weighs about 1 GB in git's store: only a job that builds the backend pays for it
    for name, job in jobs_of("ci.yml").items():
        installs = re.search(r"^\s*pip install \.\s*$", job, flags=re.MULTILINE) is not None
        fetches = re.search(r"^\s*submodules:\s*recursive\s*$", job, flags=re.MULTILINE) is not None
        assert fetches == installs, (f"the {name} job of ci.yml " + (
            "installs the package without fetching the submodule" if installs else "fetches the submodule without building it"))


def test_every_job_of_ci_but_promote_skips_the_pull_requests_of_this_repository():
    # a pull request from a branch of this repository is already covered by the push to that branch
    for name, job in jobs_of("ci.yml").items():
        if name != "promote":
            assert f"    if: {FORK_GUARD}\n" in job, f"the {name} job of ci.yml runs twice for a pull request of this repository"


def test_cd_can_be_run_by_hand_and_does_not_classify_a_tag_then():
    assert re.search(r"^  workflow_dispatch:", read(".github/workflows/cd.yml"), flags=re.MULTILINE)
    assert "if: github.event_name == 'push'" in jobs_of("cd.yml")["meta"], "a run by hand is on a branch: there is no tag to classify"


def test_the_build_of_cd_gets_the_whole_history():
    # versioneer reads the git tags to name the version
    assert re.search(r"^\s*fetch-depth:\s*0\s*$", jobs_of("cd.yml")["build"], flags=re.MULTILINE)


def test_cd_builds_the_sdist_and_nothing_else():
    # the backend is not built into a wheel (RELEASING.md says why): the build is the Python files only
    build = jobs_of("cd.yml")["build"]
    assert "python -m build --sdist" in build
    for leftover in ("docker", "manylinux", "auditwheel", "--wheel", "yum "):
        assert leftover not in build, f"the build of cd.yml still has {leftover!r}, which an sdist does not need"
    assert build.index("actions/setup-python") < build.index("python -m build"), "Python has to be set up before the build"


def test_pypi_is_for_releases_and_the_github_release_comes_after_it():
    jobs = jobs_of("cd.yml")
    for name in ("publish-to-pypi", "github-release"):
        assert "needs.meta.outputs.kind == 'release'" in jobs[name], f"{name} is not reserved for releases"
    assert "publish-to-pypi" in needs_of(jobs["github-release"]), "the GitHub release has to wait for the upload to PyPI"
    wait = re.search(r"- name: Wait for CI workflow to succeed\n\s*if:(.*)", jobs["build"])
    assert wait and "needs.meta.outputs.kind == 'prerelease'" in wait.group(1), "only a pre-release waits for CI"


def test_a_job_of_cd_that_reads_the_kind_of_the_tag_needs_the_job_that_works_it_out():
    for name, job in jobs_of("cd.yml").items():
        if "needs.meta." in job:
            assert "meta" in needs_of(job), f"{name} reads needs.meta without needing it"


# -- the supported Pythons --------------------------------------------------------

def supported_pythons():
    """The first and the last minor version of Python 3 that python_requires in setup.py allows."""
    bounds = re.search(r"python_requires\s*=\s*['\"]>=3\.(\d+),\s*<3\.(\d+)", read("setup.py"))
    return int(bounds.group(1)), int(bounds.group(2)) - 1


def test_the_classifiers_list_every_supported_python():
    # the Pythons the package declares, the ones CI tests, and the ones PyPI is told about agree
    metadata = read("src/galaxia_ananke/__metadata__.py")
    classified = set(re.findall(r'"Programming Language :: Python :: (3\.\d+)"', metadata))
    first, last = supported_pythons()
    supported = {f"3.{minor}" for minor in range(first, last + 1)}
    ci = read(".github/workflows/ci.yml")
    tested = set(re.findall(r'"(3\.\d+)"', re.search(r"python-version:\s*\[([^\]]*)\]", ci).group(1)))
    assert classified == supported == tested, (sorted(classified), sorted(supported), sorted(tested))
    assert '"Programming Language :: Python :: 3 :: Only"' in metadata


def test_the_documents_state_the_supported_python_range():
    # the range in setup.py is the one the README gives, in its introduction and in its requirements
    first, last = supported_pythons()
    readme = read("README.md")
    assert f"compatible with Python versions 3.{first} to 3.{last}" in readme
    assert f"Python 3.{first} to 3.{last} (" in readme


def test_ruff_checks_for_the_oldest_python_supported():
    # syntax of a newer Python would be accepted by ruff, and fail for the users of the oldest one
    first, _ = supported_pythons()
    assert re.search(r'^target-version\s*=\s*"py3%d"\s*$' % first, read("pyproject.toml"), flags=re.MULTILINE)


def test_cd_builds_with_a_supported_python():
    first, last = supported_pythons()
    python = int(re.search(r'python-version:\s*"3\.(\d+)"', read(".github/workflows/cd.yml")).group(1))
    assert first <= python <= last


# -- the locks of the CI and CD requirements --------------------------------------

def test_pip_installs_in_workflows_are_pinned_by_hash():
    # installing this repository itself (`pip install .`) cannot be pinned by hash; any other
    # pip install has to be of a lock (see .github/requirements/update.sh)
    unpinned = []
    for workflow in WORKFLOWS:
        for number, line in enumerate(workflow.read_text(encoding="utf-8").splitlines(), 1):
            code = line.split("#", 1)[0].strip()
            if re.search(r"\bpip install\b", code) and "--require-hashes" not in code \
                    and not re.search(r"\bpip install \.$", code):
                unpinned.append(f"{workflow.name}:{number}: {code}")
    assert not unpinned, f"pip install without --require-hashes: {unpinned}"


@pytest.mark.parametrize("workflow", WORKFLOWS, ids=lambda p: p.name)
def test_workflows_install_locks_that_exist(workflow):
    # the names are found wherever they are written, in a command or in a matrix, relative to .github/
    text = workflow.read_text(encoding="utf-8")
    for lock in set(re.findall(r"requirements(?:-frozen)?/requirements-[\w.\-]+\.txt", text)):
        assert (ROOT / ".github" / lock).is_file(), f"{workflow.name} installs {lock}, which does not exist"


def test_the_cd_build_python_is_the_one_its_lock_is_compiled_for():
    # the build tools are locked for one Python: the runner has to build with that one
    python = re.search(r'python-version:\s*"3\.(\d+)"', read(".github/workflows/cd.yml"))
    assert python, "cd.yml does not say which Python it builds with"
    assert re.search(rf"--python-version\s+3\.{python.group(1)}\b", header_of(REQUIREMENTS / "requirements-cd-build.txt")), \
        f"cd.yml builds with Python 3.{python.group(1)}: run .github/requirements/update.sh"


def test_each_python_of_the_matrix_has_the_right_test_lock():
    # the Pythons below the one requirements-test.txt is compiled for use the frozen lock, the others do not
    ci = read(".github/workflows/ci.yml")
    versions = re.findall(r'"(3\.\d+)"', re.search(r"python-version:\s*\[([^\]]*)\]", ci).group(1))
    frozen = dict(re.findall(r'-\s*python-version:\s*"(3\.\d+)"\s*\n\s*test-lock:\s*(\S+)', ci))
    lowest = int(re.search(r"--python-version\s+3\.(\d+)", header_of(REQUIREMENTS / "requirements-test.txt")).group(1))
    for version in versions:
        if int(version.split(".")[1]) < lowest:
            assert frozen.get(version) == "requirements-frozen/requirements-test.txt", \
                f"Python {version} needs the frozen lock"
        else:
            assert version not in frozen, f"Python {version} is covered by requirements-test.txt"


def test_the_frozen_lock_is_a_complete_lock():
    # Compiled with the oldest Python the package supports, it holds an answer for every Python
    # from there up: that is what lets a Python be retired to it without touching it, and what a
    # copy of the live lock (compiled for a newer Python) would lose.
    oldest = re.search(r"python_requires\s*=\s*['\"]>=3\.(\d+)", read("setup.py")).group(1)
    assert re.search(rf"--python-version\s+3\.{oldest}\b", header_of(FROZEN / "requirements-test.txt")), \
        f"requirements-frozen/requirements-test.txt has to be compiled for Python 3.{oldest}, the oldest the package supports"


def requirement_names(path):
    """The names of the packages a requirements file asks for, the files it includes with -r too."""
    names = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line.startswith("-r"):
            names |= requirement_names((path.parent / line[2:].strip()).resolve())
        elif line and not line.startswith("-"):
            names.add(re.sub(r"[-_.]+", "-", re.match(r"[A-Za-z0-9][A-Za-z0-9._-]*", line).group(0)).lower())
    return names


@pytest.mark.parametrize("source", SOURCES, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_every_requirement_is_in_its_lock(source):
    # a requirement added to a .in file, and no run of update.sh after it, would leave the tools without it
    lock = source.with_suffix(".txt").read_text(encoding="utf-8")
    pinned = {re.sub(r"[-_.]+", "-", name).lower()
              for name in re.findall(r"^([A-Za-z0-9][A-Za-z0-9._-]*)(?:\[[^\]]*\])?==", lock, flags=re.MULTILINE)}
    missing = sorted(requirement_names(source) - pinned)
    assert not missing, f"{source.name} asks for {missing}, which {source.with_suffix('.txt').name} does not pin: run update.sh"


def test_the_live_and_the_frozen_test_requirements_ask_for_the_same_tools():
    def tools(path):
        lines = path.read_text(encoding="utf-8").splitlines()
        return sorted(line.strip() for line in lines if line.strip() and not line.startswith("#"))
    assert tools(REQUIREMENTS / "requirements-test.in") == tools(FROZEN / "requirements-test.in")


def test_there_are_ci_requirements():
    assert SOURCES


@pytest.mark.parametrize("source", SOURCES, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_ci_requirements_are_locked_with_hashes(source):
    lock = source.with_suffix(".txt")
    assert lock.is_file(), f"{lock.name} is missing: run .github/requirements/update.sh"
    pins = [line for line in lock.read_text(encoding="utf-8").splitlines()
            if re.match(r"^[A-Za-z0-9_.\-]+(\[[^\]]*\])?==", line)]
    assert pins, f"{lock.name} pins nothing"
    # every pin is followed by the hashes of its files, which is what the trailing backslash says
    assert all(line.rstrip().endswith("\\") for line in pins), f"a pin without hashes in {lock.name}"
    assert "--hash=sha256:" in lock.read_text(encoding="utf-8")


@pytest.mark.parametrize("lock", LOCKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_locks_record_the_command_dependabot_runs_again(lock):
    header = header_of(lock)
    assert f"uv pip compile {lock.stem}.in" in header, f"{lock.name} does not record its uv pip compile command"
    for option in ("--universal", "--generate-hashes", "--python-version"):
        assert option in header, f"{lock.name} was not compiled with {option}: Dependabot would not keep it"


def test_update_sh_locks_every_input_with_the_python_its_lock_records():
    # a .in file nobody locks, or a lock compiled with another Python than update.sh asks for,
    # is a script and its files that have drifted apart
    script = (REQUIREMENTS / "update.sh").read_text(encoding="utf-8")
    wanted = {}
    for name, python, folder in re.findall(r"^lock (\S+) (3\.\d+)(?: (\S+))?\s*(?:#.*)?$", script, flags=re.MULTILINE):
        wanted[((REQUIREMENTS / folder).resolve() if folder else REQUIREMENTS) / f"requirements-{name}.in"] = python
    assert {p.resolve() for p in SOURCES} == {p.resolve() for p in wanted}, \
        "update.sh does not lock exactly the requirements-*.in files there are"
    for source, python in wanted.items():
        lock = source.with_suffix(".txt")
        assert re.search(rf"--python-version\s+{re.escape(python)}\b", header_of(lock)), \
            f"update.sh locks {source.name} for Python {python}, but {lock.parent.name}/{lock.name} was compiled otherwise"


def test_dependabot_updates_the_locks_with_uv():
    # The pip ecosystem would regenerate a lock with pip-compile for its own Python only, and drop
    # what the other Pythons need. The uv ecosystem runs the lock's own `uv pip compile` command again.
    text = read(".github/dependabot.yml")
    entries = [e for e in re.split(r"^\s*- package-ecosystem:\s*", text, flags=re.MULTILINE)[1:]
               if re.search(r"^\s*directory:\s*/\.github/requirements\s*$", e, flags=re.MULTILINE)]
    assert len(entries) == 1, "dependabot.yml needs exactly one entry for /.github/requirements"
    assert entries[0].split()[0].strip("\"'") == "uv", "the locks must be updated by the uv ecosystem, not pip"


def test_dependabot_proposes_its_updates_against_develop():
    # without target-branch it would target main, where nothing is merged between releases
    entries = re.split(r"^\s*- package-ecosystem:\s*", read(".github/dependabot.yml"), flags=re.MULTILINE)[1:]
    assert entries
    for entry in entries:
        assert re.search(r"^\s*target-branch:\s*develop\s*$", entry, flags=re.MULTILINE), \
            f"the {entry.split()[0]} entry of dependabot.yml does not target develop"


@pytest.mark.parametrize("lock", sorted(REQUIREMENTS.glob("requirements-*.txt")), ids=lambda p: p.name)
def test_locks_dependabot_maintains_pin_each_package_once(lock):
    # Dependabot takes one version as the current one and checks the newest can be pinned for
    # every Python the lock covers: a package locked once per Python is reported as an error.
    # Such a lock belongs in requirements-frozen/, which it does not scan.
    pins = [line for line in lock.read_text(encoding="utf-8").splitlines()
            if re.match(r"^[A-Za-z0-9_.\-]+(\[[^\]]*\])?==", line)]
    names = [re.match(r"^([A-Za-z0-9_.\-]+)", line).group(1).lower() for line in pins]
    twice = sorted({name for name in names if names.count(name) > 1})
    assert not twice, f"{lock.name} pins {twice} more than once: it cannot be maintained by Dependabot"


def test_the_folder_dependabot_scans_has_no_subfolders():
    # Dependabot fetches the subfolders of the folder it is given, not only its files. A lock it
    # must not maintain lives next to the folder, as requirements-frozen/ does.
    subfolders = [p.name for p in REQUIREMENTS.iterdir() if p.is_dir() and p.name != "__pycache__"]
    assert not subfolders, f".github/requirements/ has subfolders {subfolders}: Dependabot would scan them"
    assert not re.search(r"directory:\s*/\.github/requirements-frozen", read(".github/dependabot.yml"))


def test_the_actionlint_download_has_its_checksum():
    checksum = (REQUIREMENTS / "actionlint.sha256").read_text(encoding="utf-8").strip()
    assert re.fullmatch(r"[0-9a-f]{64}", checksum)


# -- what has to agree with something else ----------------------------------------

ZENODO_DOI = "10.5281/zenodo.13942945"   # all versions of this package
JOSS_DOI = "10.21105/joss.06234"         # the py-ananke article


def test_the_dois_agree():
    metadata = re.search(r'^__doi__ = "([^"]+)"', read("src/galaxia_ananke/__metadata__.py"), flags=re.MULTILINE)
    assert metadata and metadata.group(1) == JOSS_DOI
    cff = read("CITATION.cff")
    assert re.search(rf"^doi: {re.escape(ZENODO_DOI)}$", cff, flags=re.MULTILINE)
    assert re.search(rf"^  doi: {re.escape(JOSS_DOI)}$", cff, flags=re.MULTILINE)
    readme = read("README.md")
    assert ZENODO_DOI in readme and JOSS_DOI in readme
    related = json.loads(read(".zenodo.json"))["related_identifiers"]
    assert any(item["identifier"] == JOSS_DOI for item in related)


def test_the_license_is_said_the_same_way_everywhere():
    assert json.loads(read(".zenodo.json"))["license"] == "AGPL-3.0-or-later"
    assert re.search(r"^license: AGPL-3.0-or-later$", read("CITATION.cff"), flags=re.MULTILINE)
    assert "(AGPLv3+)" in read("src/galaxia_ananke/__metadata__.py")
    assert "GNU AFFERO GENERAL PUBLIC LICENSE" in read("LICENSE")


def test_the_license_shipped_in_the_package_is_the_one_of_the_repository():
    assert read("LICENSE") == read("src/galaxia_ananke/__license__/LICENSE")


def test_zenodo_and_citation_files_name_the_same_software():
    title = json.loads(read(".zenodo.json"))["title"]
    assert re.search(rf'^title: "{re.escape(title)}"$', read("CITATION.cff"), flags=re.MULTILINE)


def test_source_files_carry_the_license_header():
    # empty files, Sphinx's conf.py, and what versioneer generates are the exceptions; so is
    # a symlink that a checkout without symlink support holds as a line of text
    exempt = {"conf.py", "_version.py", "versioneer.py"}
    missing = []
    candidates = list((ROOT / "src").rglob("*.py")) + list((ROOT / ".github").rglob("*.py")) + [ROOT / "setup.py"]
    for path in candidates:
        if path.is_symlink() or path.name in exempt:
            continue
        text = path.read_text(encoding="utf-8")
        if len(text.splitlines()) >= 3 and "(AGPL-3.0-or-later)" not in text:
            missing.append(str(path.relative_to(ROOT)))
    assert not missing, f"no AGPL-3.0-or-later header in: {missing}"


def test_the_readme_says_what_the_build_requires():
    # the versions of the tools the installation checks for, and the Pythons, are those the code enforces
    build = read("src/_build_utils.py")
    make = re.search(r'MIN_MAKE_VERSION\s*=\s*Version\("([\d.]+)"\)', build).group(1)
    gpp = re.search(r'MIN_GPP_VERSION\s*=\s*Version\("([\d.]+)"\)', build).group(1)
    readme = read("README.md")
    assert f"GNU make ({make} or later)" in readme, f"make {make} in the README"
    assert f"g++ ({gpp} or later)" in readme, f"g++ {gpp} in the README"
    assert "Windows is not supported" in readme
    assert "pip install galaxia_ananke" in readme                    # from PyPI first


def test_python_sources_compile_without_warnings():
    # An invalid escape sequence ("\s" in a plain string) is only a warning today: a DeprecationWarning, then
    # a SyntaxWarning, and the one a future Python turns into an error. What versioneer generates is exempt.
    exempt = {"versioneer.py", "_version.py"}
    sources = (list((ROOT / "src").rglob("*.py")) + list((ROOT / ".github").rglob("*.py"))
               + [ROOT / "setup.py", ROOT / "docs" / "source" / "conf.py"])
    problems = []
    for path in sources:
        if path.is_symlink() or path.name in exempt or not path.is_file():
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            try:
                compile(path.read_text(encoding="utf-8"), str(path), "exec")
            except Exception as error:
                problems.append(f"{path.relative_to(ROOT)}: {error}")
    assert not problems, problems


# -- the changelog -----------------------------------------------------------------

def test_changelog_format():
    text = read("CHANGELOG.md")
    headings = re.findall(r"^## \[([^\]]+)\](.*)$", text, flags=re.MULTILINE)
    names = [name for name, _ in headings]
    assert names and names[0] == "Unreleased", "the first section has to be [Unreleased]"
    for name, rest in headings[1:]:
        assert re.fullmatch(r" - \d{4}-\d{2}-\d{2}", rest), f"[{name}] needs '## [{name}] - YYYY-MM-DD'"
    links = set(re.findall(r"^\[([^\]]+)\]:\s", text, flags=re.MULTILINE))
    assert set(names) <= links, f"no link at the bottom for: {sorted(set(names) - links)}"


def test_changelog_links_compare_consecutive_versions():
    # at a release `[Unreleased]` becomes the new version and the links move with it: this is the step to get wrong
    text = read("CHANGELOG.md")
    versions = re.findall(r"^## \[(\d[^\]]*)\]", text, flags=re.MULTILINE)
    links = dict(re.findall(r"^\[([^\]]+)\]:\s*(\S+)\s*$", text, flags=re.MULTILINE))
    base = "https://github.com/athob/py-Galaxia-ananke/compare/"
    assert versions, "the changelog has no released version"
    assert links["Unreleased"] == f"{base}v{versions[0]}...develop"
    for newer, older in zip(versions, versions[1:]):
        assert links[newer] == f"{base}v{older}...v{newer}", f"[{newer}] has to compare v{older} with v{newer}"
    # the oldest version compares with the pre-release its development started from
    assert re.fullmatch(re.escape(base) + r"v[\w.\-]+\.\.\.v" + re.escape(versions[-1]), links[versions[-1]])


def test_changelog_tells_what_came_before_0_4_0():
    # the package did not start with 0.4.0, and the changelog must not read as if it did
    text = read("CHANGELOG.md")
    assert re.search(r"^## Before 0\.4\.0$", text, flags=re.MULTILINE)
    for tag in ("v0.0.0-alpha", "v0.3.1-beta", "v0.4.0b3"):
        assert tag in text


# -- the contributing guide -------------------------------------------------------

def test_contributing_describes_every_job_of_ci():
    # a job added to ci.yml and left out of the guide would be a check nobody knows to run
    contributing = read("CONTRIBUTING.md")
    for name in jobs_of("ci.yml"):
        assert f"- **{name}**" in contributing, f"CONTRIBUTING.md does not describe the {name} job of ci.yml"


def test_contributing_points_only_to_files_that_exist():
    unresolved = [link for link in re.findall(r"\]\(([^)#]+)(?:#[^)]*)?\)", read("CONTRIBUTING.md"))
                  if not link.startswith("http") and not (ROOT / link).exists()]
    assert not unresolved, f"CONTRIBUTING.md points to files that do not exist: {unresolved}"


def test_contributing_is_about_this_repository():
    # the guide began as the one of the sibling project: what was not adapted would still talk about it
    contributing = prose("CONTRIBUTING.md")
    # the name of the sibling project is allowed where it is named as such, even if a line breaks inside it
    assert "enbid" not in re.sub(r"py-EnBiD-\s*ananke", "", contributing).lower(), "CONTRIBUTING.md still talks about EnBiD"
    for command in re.findall(r"gh api[^\n]*?(repos/\S+)", contributing):
        assert command.startswith("repos/athob/py-Galaxia-ananke/"), f"a gh command of CONTRIBUTING.md names {command}"


def test_contributing_says_which_pythons_use_which_test_lock():
    ci = read(".github/workflows/ci.yml")
    frozen = [version for version, _ in re.findall(r'-\s*python-version:\s*"(3\.\d+)"\s*\n\s*test-lock:\s*(\S+)', ci)]
    assert frozen, "no Python uses the frozen lock"
    listed = ", ".join(frozen[:-1]) + (" and " if len(frozen) > 1 else "") + frozen[-1]
    lowest = int(re.search(r"--python-version\s+3\.(\d+)", header_of(REQUIREMENTS / "requirements-test.txt")).group(1))
    oldest = supported_pythons()[0]
    contributing = prose("CONTRIBUTING.md")
    assert f"(3.{lowest} and up)" in contributing, "the Pythons of requirements-test.txt"
    assert f"oldest the package supports (3.{oldest})" in contributing
    assert f"at the end of their life ({listed})" in contributing, f"the Pythons of the frozen lock are {listed}"


def test_contributing_states_the_supported_python_range():
    first, last = supported_pythons()
    contributing = prose("CONTRIBUTING.md")
    assert f"(currently 3.{first} up to 3.{last})" in contributing
    assert f"each Python version from 3.{first} to 3.{last}," in contributing


def test_contributing_says_vaex_is_pinned_if_it_is():
    # the day the pin becomes a range, the guide has to stop saying so
    pinned = re.search(r"['\"]vaex==[\d.]+['\"]", read("setup.py")) is not None
    assert ("vaex is pinned to one version" in prose("CONTRIBUTING.md")) == pinned
