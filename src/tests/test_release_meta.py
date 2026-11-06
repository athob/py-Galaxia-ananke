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
Tests of .github/scripts/release_meta.py, which decides where a tag is published to
and what its release notes are.
"""
import importlib.util
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("release_meta", ROOT / ".github" / "scripts" / "release_meta.py")
release_meta = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(release_meta)

CHANGELOG = """\
# Changelog

## [Unreleased]

- not yet

## [1.1.0] - 2026-02-01

### Added

- the new thing

## [1.0.0] - 2026-01-01

- first

[Unreleased]: https://example.org/compare/v1.1.0...HEAD
[1.1.0]: https://example.org/compare/v1.0.0...v1.1.0
[1.0.0]: https://example.org/releases/tag/v1.0.0
"""


@pytest.mark.parametrize("tag, kind, version", [
    ("v0.3.4", "release", "0.3.4"),
    ("v1", "release", "1"),
    ("v0.3.4.post1", "release", "0.3.4.post1"),
    ("v0.3.4a1", "prerelease", "0.3.4a1"),
    ("v0.3.4b2", "prerelease", "0.3.4b2"),
    ("v0.3.4rc0", "prerelease", "0.3.4rc0"),
    ("v0.3.4.dev0", "dev", "0.3.4.dev0"),
    ("v0.3.4rc0.dev1", "dev", "0.3.4rc0.dev1"),
])
def test_classify(tag, kind, version):
    assert release_meta.classify(tag) == (kind, version)


@pytest.mark.parametrize("tag", [
    "v0.3.4-test", "v0.3.2rc0-test", "v0.3.4dev1", "v0.3.4.dev", "vfoo", "v0.3.4.foo", "0.3.4", "v", "v0.3.4;rm",
])
def test_classify_refuses(tag):
    with pytest.raises(release_meta.ReleaseError):
        release_meta.classify(tag)


def test_changelog_section_stops_at_the_next_section_and_at_the_links():
    rest, body = release_meta.changelog_section(CHANGELOG, "1.1.0")
    assert rest == " - 2026-02-01"
    assert body == "### Added\n\n- the new thing"
    assert release_meta.changelog_section(CHANGELOG, "1.0.0")[1] == "- first"
    assert release_meta.changelog_section(CHANGELOG, "9.9.9") is None


def test_changelog_section_stops_at_any_second_level_heading():
    text = "## [1.0.0] - 2026-01-01\n\n- a\n\n## Before 1.0.0\n\nnot part of it\n\n[1.0.0]: https://example.org\n"
    assert release_meta.changelog_section(text, "1.0.0")[1] == "- a"


def test_release_section():
    assert release_meta.release_section("1.1.0", CHANGELOG).startswith("### Added")


def test_release_section_refuses_a_missing_undated_or_empty_section():
    with pytest.raises(release_meta.ReleaseError):
        release_meta.release_section("2.0.0", CHANGELOG)
    with pytest.raises(release_meta.ReleaseError):
        release_meta.release_section("2.0.0", "## [2.0.0]\n\n- x\n")
    with pytest.raises(release_meta.ReleaseError):
        release_meta.release_section("2.0.0", "## [2.0.0] - 2026-03-01\n\n## [1.0.0] - 2026-01-01\n\n- x\n")


def test_every_released_version_of_the_real_changelog_can_be_released():
    # the notes of every version in CHANGELOG.md can be extracted: dated, not empty
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    versions = [line.split("]")[0][4:] for line in text.splitlines()
                if line.startswith("## [") and not line.startswith("## [Unreleased]")]
    assert versions
    for version in versions:
        assert release_meta.release_section(version, text)


def test_the_notes_of_0_4_0_stop_before_the_history_that_precedes_it():
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    notes = release_meta.release_section("0.4.0", text)
    assert "versioneer" in notes
    assert "alpha" not in notes and "Before 0.4.0" not in notes
