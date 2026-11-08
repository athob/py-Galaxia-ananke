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
Tests of the metadata of the package: the version and the date that versioneer gives (and the date of now when
it gives none, which is what building from files without git information needs), what is printed about the project,
and the classifiers. The metadata is imported afresh, in a package of its own, with the answer of versioneer that
each test sets.
"""
import datetime
import importlib
import re

import pytest

from .utils import fresh_package

UTC = datetime.timezone.utc


def metadata_for(monkeypatch, version="1.2.3", date="2026-11-07T03:16:56+0000"):
    """The metadata as a fresh import makes it when versioneer answers with this version and this date."""
    with fresh_package() as package:
        versioneer = importlib.import_module(f"{package}._version")
        monkeypatch.setattr(versioneer, "get_versions", lambda: {"version": version, "date": date})
        yield importlib.import_module(f"{package}.__metadata__")


@pytest.fixture
def make(monkeypatch):
    generators = []

    def build(**kwargs):
        generator = metadata_for(monkeypatch, **kwargs)
        generators.append(generator)
        return next(generator)

    yield build
    for generator in generators:
        generator.close()


def one_line(text):
    return re.sub(r"\s+", " ", text)


# -- the version and the date -------------------------------------------------------

def test_the_version_is_the_one_versioneer_gives(make):
    assert make(version="0.8.0+69.gfbdc9b2").__version__ == "0.8.0+69.gfbdc9b2"


def test_the_date_is_the_one_versioneer_gives_with_its_offset(make):
    assert make(date="2026-11-07T03:16:56+0000").__date__ == datetime.datetime(2026, 11, 7, 3, 16, 56, tzinfo=UTC)
    shifted = make(date="2026-11-07T03:16:56+0530").__date__
    assert shifted.utcoffset() == datetime.timedelta(hours=5, minutes=30)
    assert shifted == datetime.datetime(2026, 11, 6, 21, 46, 56, tzinfo=UTC)


def test_without_a_date_it_is_the_date_of_now_in_utc(make):
    before = datetime.datetime.now(UTC)
    metadata = make(version="0+unknown", date=None)
    after = datetime.datetime.now(UTC)
    assert metadata.__version__ == "0+unknown"
    assert metadata.__date__.tzinfo == UTC
    assert before <= metadata.__date__ <= after


def test_the_copyright_runs_from_the_year_of_the_project_to_the_year_of_the_date(make):
    assert make(date="2031-02-03T04:05:06+0000").__copyright_short__ == "Copyright (C) 2022-2031  Adrien Thob"


def test_the_copyright_of_a_package_without_a_date_ends_this_year(make):
    year = datetime.datetime.now(UTC).year
    assert make(version="0+unknown", date=None).__copyright_short__ == f"Copyright (C) 2022-{year}  Adrien Thob"


# -- what is printed about the project ----------------------------------------------

def test_the_copyright_text_names_the_project_its_version_and_its_license(make):
    metadata = make(version="9.8.7")
    text = one_line(metadata.__copyright__)
    assert f"{metadata.__project__} v9.8.7" in text
    assert "GNU Affero General Public License" in text
    assert one_line(metadata.__copyright_short__) in text


def test_the_credits_thank_the_contributors_and_the_supporters_together(make):
    text = one_line(make().__credits__)
    assert "Robyn Sanderson & Andrew Eden" in text
    assert "Farnik Nikakhtar, Nondh Panithanpaisal, Nicolas Garavito-Camargo & Sanjib Sharma" in text


def test_the_citation_gives_the_doi_the_reference_and_the_bibtex_entry(make):
    metadata = make()
    text = one_line(metadata.__citing__)
    assert f"please use DOI {metadata.__doi__} as a reference" in text      # the citation holds the DOI too: this is the sentence
    assert one_line(metadata.__citation__) in text
    assert metadata.__bibtex__.splitlines()[0] in metadata.__citing__
    assert "{__citation__}" not in metadata.__citing__ and "{__bibtex__}" not in metadata.__citing__


def test_the_readme_is_the_copyright_and_the_credits(make):
    metadata = make()
    assert metadata.__copyright__ in metadata.__readme__
    assert metadata.__credits__ in metadata.__readme__


@pytest.mark.parametrize("name, text", [
    ("copyright", "__copyright__"),
    ("credits", "__credits__"),
    ("license", "__license_full__"),
    ("readme", "__readme__"),
])
def test_the_printers_are_those_of_the_builtins_and_print_their_text(make, name, text):
    # copyright(), credits(), license() and readme() print the texts, and are shown as the builtins of Python are
    metadata = make()
    printer = getattr(metadata, name)
    assert type(printer).__name__ == "_Printer"
    assert printer._Printer__name == name
    assert printer._Printer__data == getattr(metadata, text)


def test_everything_in_all_is_defined(make):
    metadata = make()
    assert [name for name in metadata.__all__ if not hasattr(metadata, name)] == []


# -- the classifiers ----------------------------------------------------------------

def test_the_license_of_the_classifiers_is_the_license_of_the_package(make):
    metadata = make()
    assert metadata.__license_classifier__ in metadata.__classifiers__
    assert metadata.__license_classifier__ == f"License :: OSI Approved :: {metadata.__license__}"


def test_the_status_is_among_the_classifiers(make):
    metadata = make()
    assert metadata.__status_classifier__ in metadata.__classifiers__


if __name__ == '__main__':
    pass
