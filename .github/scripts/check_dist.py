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
Check the sdist that `python -m build --sdist` made, before anyone publishes it.

    python .github/scripts/check_dist.py dist

py-Galaxia-ananke is published as an sdist and as nothing else (RELEASING.md says why): a wheel
would not carry the C++ backend, which installing builds into a cache directory. What the
installer needs from the sdist is its Python sources, `.gitmodules`, which gives the url of the
backend, and the `_version.py` that versioneer wrote, whose full revision id is what the
installer asks GitHub for, to find the commit of the backend this version pins. An sdist whose
revision id is unknown (one built from files without git information) cannot be installed. The
backend is not part of the sdist, and must not become part of it: it is about 5.8 GB.
"""
import json
import re
import sys
import tarfile
from pathlib import Path

PACKAGE = "galaxia_ananke"
SUBMODULE = "Galaxia-ananke"

SDIST_REQUIRED = [  # relative to the directory the sdist unpacks into
    "setup.py",
    "setup.cfg",
    "pyproject.toml",
    "MANIFEST.in",
    "README.md",
    "LICENSE",
    ".gitmodules",
    "versioneer.py",
    "src/_build_utils.py",
    f"src/{PACKAGE}/__init__.py",
    f"src/{PACKAGE}/_version.py",
    f"src/{PACKAGE}/__license__/LICENSE",
]


def normalized(version):
    return re.sub(r"[-_.]+", ".", version.lower())


def baked_version(text):
    """What versioneer wrote in the sdist's _version.py: the JSON between its VERSION_JSON markers."""
    match = re.search(r"version_json = '''\n(.*?)'''  # END VERSION_JSON", text, flags=re.DOTALL)
    return json.loads(match.group(1)) if match else None


def check_sdist(path):
    problems = []
    with tarfile.open(path) as sdist:
        members = [n for n in sdist.getnames() if "/" in n]
        # the first path component is the directory it unpacks into
        root = sdist.getnames()[0].split("/", 1)[0]
        names = {n.split("/", 1)[1] for n in members}
        version_file = f"src/{PACKAGE}/_version.py"
        text = sdist.extractfile(f"{root}/{version_file}").read().decode("utf-8") if version_file in names else ""
    for name in SDIST_REQUIRED:
        if name not in names:
            problems.append(f"{path.name}: missing {name}")
    for name in sorted(n for n in names if n == SUBMODULE or n.startswith(f"{SUBMODULE}/")):
        problems.append(f"{path.name}: carries {name}, and the backend is fetched when the package is installed")
    if text:
        baked = baked_version(text)
        named = re.fullmatch(rf"{PACKAGE}-(.+)\.tar\.gz", path.name)
        if baked is None:
            problems.append(f"{path.name}: {version_file} has no version baked in by versioneer")
        else:
            if baked.get("error") is not None:
                problems.append(f"{path.name}: versioneer reports an error: {baked['error']}")
            if not re.fullmatch(r"[0-9a-f]{40}", str(baked.get("full-revisionid"))):
                problems.append(f"{path.name}: the revision id baked in is {baked.get('full-revisionid')!r}, not a commit: "
                                f"the installer could not find the backend this version pins")
            if baked.get("dirty"):
                problems.append(f"{path.name}: it was built from files that differ from the commit (dirty)")
            if baked.get("date") is None:
                problems.append(f"{path.name}: no date baked in")
            if named and normalized(baked.get("version", "")) != normalized(named.group(1)):
                problems.append(f"{path.name}: the version baked in is {baked.get('version')}, not the one of its name")
    return problems


def main(argv):
    dist = Path(argv[1] if len(argv) > 1 else "dist")
    wheels, sdists = sorted(dist.glob("*.whl")), sorted(dist.glob("*.tar.gz"))
    problems = []
    if wheels or len(sdists) != 1:
        problems.append(f"expected one sdist and no wheel in {dist}/, found {[p.name for p in sdists]} "
                        f"and {[p.name for p in wheels]}")
    for sdist in sdists:
        problems += check_sdist(sdist)
    for problem in problems:
        print(f"::error::{problem}")
    if not problems:
        print(f"OK: {sdists[0].name} has what it should.")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
