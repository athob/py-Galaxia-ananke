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
Tests of the repository itself rather than of the code: what has to agree with something
else does (the locks of the CI requirements and the workflows that install them, the
matrix and the lock each of its Pythons uses). Standard library only, and nothing here
needs Galaxia.
"""
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
WORKFLOWS = sorted((ROOT / ".github" / "workflows").glob("*.yml"))

REQUIREMENTS = ROOT / ".github" / "requirements"
FROZEN = ROOT / ".github" / "requirements-frozen"
LOCKS = sorted(REQUIREMENTS.glob("requirements-*.txt")) + sorted(FROZEN.glob("requirements-*.txt"))
SOURCES = sorted(REQUIREMENTS.glob("requirements-*.in")) + sorted(FROZEN.glob("requirements-*.in"))


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def header_of(lock):
    """The first lines of a lock, where uv records the command that wrote it."""
    return "\n".join(lock.read_text(encoding="utf-8").splitlines()[:3])


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


def test_the_actionlint_download_has_its_checksum():
    checksum = (REQUIREMENTS / "actionlint.sha256").read_text(encoding="utf-8").strip()
    assert re.fullmatch(r"[0-9a-f]{64}", checksum)
