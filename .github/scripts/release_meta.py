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
What the CD workflow needs to know about a tag, and about the changelog.

    release_meta.py check TAG [--changelog CHANGELOG.md] [--github-output FILE]
    release_meta.py notes TAG [--changelog CHANGELOG.md] [--output FILE]

``check`` classifies the tag, and for a release makes sure the changelog has the
section of that version, dated, so that nothing is published before it is written.
``notes`` prints the release notes, which are that section.

After the leading ``v`` a tag has to be a plain PEP 440 version:

    release      v0.3.4  v0.3.4.post1                  TestPyPI, then PyPI
    pre-release  v0.3.4a1  v0.3.4b1  v0.3.4rc0         TestPyPI, after CI has passed
    development  v0.3.4.dev0  v0.3.4rc0.dev1           TestPyPI, straight away

Anything else (``v0.3.4-test``, ``v0.3.4dev1``, ...) is refused, so that nothing
that is not a clean release can reach PyPI. Only the standard library is used.
"""
import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

TAG_RE = re.compile(r"^v(?P<version>[0-9]+(\.[0-9]+)*((a|b|rc)[0-9]+)?(\.post[0-9]+)?(\.dev[0-9]+)?)$")
PRE_RE = re.compile(r"^[0-9]+(\.[0-9]+)*(a|b|rc)[0-9]+")

SECTION_RE = re.compile(r"^## \[(?P<name>[^\]]+)\](?P<rest>.*)$")
LINK_REF_RE = re.compile(r"^\[[^\]]+\]:\s")
DATE_RE = re.compile(r"^ - [0-9]{4}-[0-9]{2}-[0-9]{2}$")


class ReleaseError(Exception):
    pass


def classify(tag: str) -> Tuple[str, str]:
    """Return (kind, version) for a tag: kind is release, prerelease or dev."""
    match = TAG_RE.match(tag)
    if not match:
        raise ReleaseError(
            f"The tag '{tag}' is not a version the CD workflow knows how to publish. "
            "Use v<N>[.<N>...] with an optional a<N>, b<N>, rc<N>, .post<N> or .dev<N> "
            "suffix, e.g. v0.3.4rc0 or v0.3.4.dev0.")
    version = match.group("version")
    if ".dev" in version:
        return "dev", version
    if PRE_RE.match(version):
        return "prerelease", version
    return "release", version


def changelog_section(text: str, name: str) -> Optional[Tuple[str, str]]:
    """The (heading remainder, body) of the ``## [name]`` section, or None.

    The remainder is what follows the bracketed name on the heading line (the date);
    the body is trimmed, and stops before the next second-level heading (whether or not
    it is a version: the history before 0.2.0 has one of its own) or the link references.
    """
    lines: List[str] = text.splitlines()
    for i, line in enumerate(lines):
        match = SECTION_RE.match(line)
        if match and match.group("name") == name:
            body: List[str] = []
            for following in lines[i + 1:]:
                if following.startswith("## ") or LINK_REF_RE.match(following):
                    break
                body.append(following)
            return match.group("rest"), "\n".join(body).strip()
    return None


def release_section(version: str, changelog_text: str) -> str:
    """The body of the dated changelog section of a release, or ReleaseError."""
    found = changelog_section(changelog_text, version)
    if found is None:
        raise ReleaseError(
            f"CHANGELOG.md has no '## [{version}]' section. Rename '## [Unreleased]' to "
            f"'## [{version}] - YYYY-MM-DD' (and add a fresh '## [Unreleased]') before tagging.")
    rest, body = found
    if not DATE_RE.match(rest):
        raise ReleaseError(f"The CHANGELOG.md heading of {version} must read "
                           f"'## [{version}] - YYYY-MM-DD', not '## [{version}]{rest}'.")
    if not body:
        raise ReleaseError(f"The CHANGELOG.md section of {version} is empty.")
    return body


def read_changelog(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ReleaseError(f"cannot read {path}: {exc}")


def cmd_check(args: argparse.Namespace) -> None:
    kind, version = classify(args.tag)
    if kind == "release":
        release_section(version, read_changelog(args.changelog))
    print(f"{args.tag}: version {version}, {kind}")
    if args.github_output:
        with open(args.github_output, "a", encoding="utf-8") as handle:
            handle.write(f"kind={kind}\nversion={version}\n")


def cmd_notes(args: argparse.Namespace) -> None:
    kind, version = classify(args.tag)
    if kind != "release":
        raise ReleaseError(f"{args.tag} is a {kind} version: only releases have release notes.")
    notes = release_section(version, read_changelog(args.changelog)) + "\n"
    if args.output:
        Path(args.output).write_text(notes, encoding="utf-8")
    else:
        sys.stdout.write(notes)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("check", "notes"):
        p = sub.add_parser(name)
        p.add_argument("tag")
        p.add_argument("--changelog", type=Path, default=Path("CHANGELOG.md"))
    sub.choices["check"].add_argument("--github-output")
    sub.choices["notes"].add_argument("--output")
    args = parser.parse_args(argv)
    try:
        {"check": cmd_check, "notes": cmd_notes}[args.command](args)
    except ReleaseError as exc:
        print(f"::error::{exc}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
