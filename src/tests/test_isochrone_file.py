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
Tests of IsochroneFile, which reads the files of the isochrones of a photometric system (a commented header with the names of
the columns, then a row for each star) and writes them for the isochrones that are made in Python.
"""
import types

import pandas as pd
import pytest
from astropy.table import Table

from ..galaxia_ananke.photometry.IsochroneFile import IsochroneFile
from .utils import write_isochrone_file

ISOCHRONE = types.SimpleNamespace(_age='Age', _mini='M_ini')            # what IsochroneFile asks of the isochrone that it belongs to
HEADER = ['Age', 'M_ini', 'M_act', 'G']


# -- how an isochrone file is made --------------------------------------------------

def test_an_isochrone_file_is_a_path_that_belongs_to_an_isochrone(tmp_path):
    iso_file = IsochroneFile(tmp_path / 'output_0.019000.dat', isochrone=ISOCHRONE)
    assert iso_file.path == tmp_path / 'output_0.019000.dat'
    assert iso_file.filename == 'output_0.019000.dat'
    assert iso_file.isochrone is ISOCHRONE
    assert repr(iso_file) == 'IsochroneFile(filename=output_0.019000.dat)'


def test_an_isochrone_file_can_be_given_as_a_string(tmp_path):
    assert IsochroneFile(str(tmp_path / 'output_0.019000.dat'), isochrone=ISOCHRONE).path == tmp_path / 'output_0.019000.dat'


def test_the_name_of_the_files_is_output_then_the_metallicity():
    assert IsochroneFile._file_format == 'output_{}.dat'


@pytest.mark.parametrize("arguments, error", [
    ((), "Keyword argument 'isochrone' missing"),
    (('a', 'b', 'c'), "Too many arguments \\(3 given\\)"),
])
def test_an_isochrone_file_that_is_wrongly_made_says_what_is_wrong(arguments, error):
    with pytest.raises(TypeError, match=error):
        IsochroneFile(*arguments)


def test_an_isochrone_file_needs_at_least_a_path():
    with pytest.raises(TypeError, match="requires at least one argument"):
        IsochroneFile(isochrone=ISOCHRONE)


# -- the metallicity ----------------------------------------------------------------

@pytest.mark.parametrize("filename, metallicity", [
    ('output_0.019000.dat', 0.019),
    ('output_0.000100.dat', 0.0001),
    ('output_0.0152.dat', 0.0152),
    ('output_0.030000.dat', 0.03),
])
def test_the_metallicity_is_in_the_name_of_the_file(filename, metallicity, tmp_path):
    assert IsochroneFile(tmp_path / filename, isochrone=ISOCHRONE).metallicity == metallicity


# -- writing a table ----------------------------------------------------------------

ROWS = {'Age': [2.0, 1.0, 1.0], 'M_ini': [1.0, 2.0, 1.0], 'M_act': [0.3, 0.2, 0.1], 'G': [5.0, 6.0, 7.0]}


@pytest.mark.parametrize("kind", ["pandas", "astropy"])
def test_a_table_is_written_in_the_order_of_the_ages_and_then_the_masses(kind, tmp_path):
    data = pd.DataFrame(ROWS) if kind == "pandas" else Table(ROWS)
    path = tmp_path / 'output_0.019000.dat'
    IsochroneFile(path, data, isochrone=ISOCHRONE)
    read = IsochroneFile(path, isochrone=ISOCHRONE)
    assert read.column_names == HEADER
    assert read.data['Age'].tolist() == [1.0, 1.0, 2.0]
    assert read.data['M_ini'].tolist() == [1.0, 2.0, 1.0]
    assert read.data['M_act'].tolist() == [0.1, 0.2, 0.3]                       # the rows kept together
    assert read.data['G'].tolist() == [7.0, 6.0, 5.0]


def test_the_header_of_a_written_file_is_a_comment_with_the_names(tmp_path):
    path = tmp_path / 'output_0.019000.dat'
    IsochroneFile(path, pd.DataFrame(ROWS), isochrone=ISOCHRONE)
    assert path.read_text().splitlines()[0] == '# Age M_ini M_act G'


@pytest.mark.parametrize("data", [{'Age': [1.0]}, [[1.0, 2.0]], 'text', None])
def test_data_that_is_not_a_table_is_refused(data, tmp_path):
    with pytest.raises(ValueError, match="should either be a pandas DataFrame or an astropy Table"):
        IsochroneFile(tmp_path / 'output_0.019000.dat', data, isochrone=ISOCHRONE)


# -- reading a file -----------------------------------------------------------------

@pytest.mark.parametrize("header", [
    "# a b c",
    "#a\tb\tc",
    "#   a    b   c   ",
    "# first line\n# a b c",                                # the last line of comments is the one with the names
    "# a\t\tb c",
], ids=["spaces", "tabs", "many spaces", "two comments", "mixed"])
def test_the_names_of_the_columns_are_in_the_last_comment_of_the_header(header, tmp_path):
    path = tmp_path / 'output_0.019000.dat'
    path.write_text(f"{header}\n1 2 3\n4 5 6\n")
    assert IsochroneFile(path, isochrone=ISOCHRONE).column_names == ['a', 'b', 'c']


@pytest.mark.xfail(raises=AssertionError, reason="the spaces around the names are stripped before the tabs are made spaces, so a "
                   "tab right after the # leaves an empty name in front of the others, and every column is one off")
def test_a_header_that_starts_with_a_tab_has_no_empty_name(tmp_path):
    path = tmp_path / 'output_0.019000.dat'
    path.write_text("#\ta\tb\tc\n1 2 3\n")
    assert IsochroneFile(path, isochrone=ISOCHRONE).column_names == ['a', 'b', 'c']


def test_the_data_of_a_file_has_the_names_of_its_header(tmp_path):
    path = write_isochrone_file(tmp_path / 'output_0.019000.dat', ['a', 'b', 'c'], [[1, 2.5, 3], [4, 5.5, 6]])
    data = IsochroneFile(path, isochrone=ISOCHRONE).data
    assert data.colnames == ['a', 'b', 'c']
    assert data['b'].tolist() == [2.5, 5.5] and data['a'].tolist() == [1, 4]


if __name__ == '__main__':
    pass
