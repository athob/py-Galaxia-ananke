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
Tests of Isochrone, which is the isochrones of one photometric system: a directory with a descriptor (the name of the system, where
the magnitudes start in the files, and their names) and a file per metallicity. They cover how one is made, from a directory or from
data in Python, what the descriptor says, which layout the files have, and the tables that are made from them. The isochrones that
are made in Python go to a temporary directory and not to the installation, and the others are small ones that are written here.
"""
import pathlib
import re
import sys

import numpy as np
import pytest

from ..galaxia_ananke._constants import CUSTOM_PHOTOCAT
from ..galaxia_ananke.photometry.Formatting import (default_formatting, oldpadova_fomatting, oldpadova_fomatting_withlogage,
                                                    padova_formatting)
from ..galaxia_ananke.photometry.Isochrone import Isochrone
from ..galaxia_ananke.photometry.IsochroneFile import IsochroneFile
from .utils import custom_isochrone_data as custom_data, write_isochrone_file

ISOCHRONE_MODULE = sys.modules[Isochrone.__module__]
HEADER = ['Age', 'M_ini', 'M_act', 'Lum', 'T_eff', 'Grav', 'Label', 'G', 'BP']          # the layout of default_formatting, with a column that it does not take


@pytest.fixture
def isochrones_path(tmp_path, monkeypatch):
    """Where the isochrones that are made in Python go: a temporary directory, in Isochrone, which took it from the constants."""
    monkeypatch.setattr(ISOCHRONE_MODULE, 'ISOCHRONES_PATH', tmp_path / 'Isochrones')
    return tmp_path / 'Isochrones'


def isochrone_directory(base, name, stars_by_metallicity, descriptor='X 9 7 2 G BP', category='mine'):
    """An isochrone as a directory of the installation has it, with the layout of default_formatting."""
    directory = base / category / name
    directory.mkdir(parents=True)
    if descriptor is not None:
        (directory / Isochrone._file_descriptor).write_text(descriptor + "\n\n")
    for metallicity, stars in stars_by_metallicity.items():
        write_isochrone_file(directory / f"output_{metallicity:.6f}.dat", HEADER, stars)
    return directory


# -- what is needed -----------------------------------------------------------------

def test_what_an_isochrone_needs():
    assert Isochrone._required_keys == ['Age', 'M_ini', 'M_act', 'Lum', 'T_eff', 'Grav']
    assert Isochrone._file_descriptor == 'IsoFileDescriptor.txt'
    needed = Isochrone._required_metallicities
    assert len(needed) == 34 and min(needed) == 0.0001 and max(needed) == 0.03
    assert len({format(metallicity, '.6f') for metallicity in needed}) == 34         # each one has a file name of its own


# -- how an isochrone is made -------------------------------------------------------

def test_an_isochrone_made_from_a_directory_is_that_directory(tmp_path):
    iso = Isochrone(tmp_path / 'padova' / 'LSST')
    assert iso.path == tmp_path / 'padova' / 'LSST' and iso.category == 'padova' and iso.name == 'LSST'
    assert repr(iso) == 'Isochrone(category=padova, name=LSST)'
    assert Isochrone(str(tmp_path / 'padova' / 'LSST')).path == tmp_path / 'padova' / 'LSST'


def test_an_isochrone_needs_at_least_an_argument():
    with pytest.raises(TypeError, match="requires at least one argument"):
        Isochrone()


def test_an_isochrone_does_not_take_more_than_two_arguments():
    with pytest.raises(TypeError, match=r"Too many arguments \(3 given\)"):
        Isochrone('a', 'b', 'c')


def test_an_isochrone_made_from_data_is_made_in_the_custom_category(isochrones_path):
    iso = Isochrone('mine', custom_data())
    assert iso.path == isochrones_path / CUSTOM_PHOTOCAT / 'mine'
    assert (iso.category, iso.name) == (CUSTOM_PHOTOCAT, 'mine')


def test_an_isochrone_made_from_data_has_a_file_for_each_metallicity_and_a_descriptor(isochrones_path):
    iso = Isochrone('mine', custom_data())
    assert sorted(path.name for path in iso.path.glob('output_*.dat')) == sorted(
        f"output_{metallicity:.6f}.dat" for metallicity in Isochrone._required_metallicities)
    assert iso.has_file_descriptor and iso.file_descriptor_path == iso.path / 'IsoFileDescriptor.txt'


def test_the_descriptor_gives_the_name_the_columns_and_the_names_of_the_magnitudes_in_order(isochrones_path):
    iso = Isochrone('mine', custom_data())
    assert iso.file_descriptor_path.read_text() == "Python_mine 8 6 2 BP G\n\n"           # the magnitudes after the six others, sorted


def test_the_files_have_the_columns_in_the_order_that_the_descriptor_says(isochrones_path):
    iso = Isochrone('mine', custom_data())
    assert iso.isochrone_files[0].column_names == ['Age', 'M_ini', 'M_act', 'Lum', 'T_eff', 'Grav', 'BP', 'G']


def test_an_isochrone_that_exists_is_not_made_again(isochrones_path):
    first = Isochrone('mine', custom_data())
    with pytest.raises(FileExistsError, match=f"Isochrone 'mine' already exists at '{re.escape(str(first.path))}'"):
        Isochrone('mine', custom_data())


def test_an_isochrone_that_exists_can_be_made_again_over_the_old_one(isochrones_path):
    first = Isochrone('mine', custom_data())
    (first.path / 'left_by_the_first.txt').write_text('')
    again = Isochrone('mine', custom_data(), overwrite=True)
    assert again.path == first.path and not (again.path / 'left_by_the_first.txt').exists() and again.has_file_descriptor


# -- what is refused ----------------------------------------------------------------

def refused(isochrones_path, data, error):
    with pytest.raises(ValueError, match=error):
        Isochrone('mine', data)
    assert not (isochrones_path / CUSTOM_PHOTOCAT / 'mine').exists()                 # nothing is left behind
    old = Isochrone('mine', custom_data())                                           # and the one that is to be replaced is not taken away
    before = sorted((path.name, path.read_bytes()) for path in old.path.iterdir())
    with pytest.raises(ValueError, match=error):
        Isochrone('mine', data, overwrite=True)
    assert sorted((path.name, path.read_bytes()) for path in old.path.iterdir()) == before


def test_data_that_misses_a_metallicity_is_refused(isochrones_path):
    data = custom_data()
    del data[0.0001]
    refused(isochrones_path, data, "Given isochrone data covers wrong set of metallicities: misses {0.0001}")


def test_data_with_a_metallicity_that_is_not_needed_is_refused(isochrones_path):
    data = custom_data()
    data[0.5] = data[0.03]
    refused(isochrones_path, data, "misincludes {0.5}")


def test_data_whose_tables_do_not_have_the_same_columns_is_refused(isochrones_path):
    data = custom_data()
    data[0.03] = data[0.03].assign(extra=1.0)
    refused(isochrones_path, data, "Given isochrone dataframes have unequal headers")


def test_data_without_a_required_column_is_refused(isochrones_path):
    data = {metallicity: table.drop(columns='Lum') for metallicity, table in custom_data().items()}
    refused(isochrones_path, data, "incomplete headers: missing {'Lum'}")


def test_data_without_magnitudes_is_refused(isochrones_path):
    data = {metallicity: table.drop(columns=['G', 'BP']) for metallicity, table in custom_data().items()}
    refused(isochrones_path, data, "have no magnitude columns")


def test_a_name_that_was_refused_can_be_used_again_without_overwrite(isochrones_path):
    data = custom_data()
    incomplete = {metallicity: table for metallicity, table in data.items() if metallicity != 0.0001}
    with pytest.raises(ValueError):
        Isochrone('mine', incomplete)
    assert Isochrone('mine', data).has_file_descriptor


def test_what_cannot_be_written_leaves_nothing_behind_either(isochrones_path, monkeypatch):
    def cannot_write(*args, **kwargs):
        raise OSError("no room left")
    monkeypatch.setattr(IsochroneFile, '__init__', cannot_write)
    with pytest.raises(OSError, match="no room left"):
        Isochrone('mine', custom_data())
    assert not (isochrones_path / CUSTOM_PHOTOCAT / 'mine').exists()


# -- what the descriptor says -------------------------------------------------------

def test_the_descriptor_is_read_as_the_name_the_columns_the_start_of_the_magnitudes_and_their_names(isochrones_path):
    iso = Isochrone('mine', custom_data())
    assert iso.file_descriptor_content == ['Python_mine', '8', '6', '2', 'BP', 'G']
    assert (iso._fd_startid, iso._fd_nmags, iso.mag_names) == (6, 2, ['BP', 'G'])
    assert iso.magnitude_itemgetter([f"column{i}" for i in range(8)]) == ('column6', 'column7')


def test_a_directory_without_a_descriptor_is_not_a_photometric_system(tmp_path):
    (tmp_path / 'padova' / 'Nothing').mkdir(parents=True)
    iso = Isochrone(tmp_path / 'padova' / 'Nothing')
    assert not iso.has_file_descriptor
    with pytest.raises(NotImplementedError, match="doesn't have an IsoFileDescriptor.txt file"):
        iso.file_descriptor_content


def test_the_magnitudes_of_a_system_with_a_single_one_are_a_tuple(tmp_path):
    directory = isochrone_directory(tmp_path, 'Single', {0.0152: []}, descriptor='Single 8 7 1 G')
    assert Isochrone(directory).magnitude_itemgetter([f"column{i}" for i in range(8)]) == ('column7',)


def test_the_isochrone_files_are_those_of_the_directory_in_the_order_of_their_names(tmp_path):
    stars = [(9.0, 1.0, 0.99, 0.0, 3.76, 4.4, 0, 5.0, 5.5)]
    iso = Isochrone(isochrone_directory(tmp_path, 'X', {0.152: stars, 0.0152: stars, 0.00152: stars}))
    assert [iso_file.filename for iso_file in iso.isochrone_files] == ['output_0.001520.dat', 'output_0.015200.dat', 'output_0.152000.dat']
    assert all(iso_file.isochrone is iso for iso_file in iso.isochrone_files)


# -- the layout of the files --------------------------------------------------------

@pytest.mark.parametrize("category, name, expected", [
    ('padova', 'WFIRST+HST__WFC3', oldpadova_fomatting_withlogage),
    *[('padova', name, oldpadova_fomatting) for name in ('WFIRST', 'LSST', 'GAIA__0', 'GAIA__DR2', 'HST__WFC3')],
    *[('padova', name, padova_formatting) for name in ('CTIO__DECam', 'LSST_DP0', 'Roman', 'Euclid', 'JWST__NIRCam', 'JWST__MIRI')],
    ('padova', 'Another', default_formatting),
    ('py_custom', 'LSST', default_formatting),                      # the names are those of padova only
    ('elsewhere', 'GAIA__DR2', default_formatting),
])
def test_the_layout_of_the_files_is_that_of_the_category_and_the_name(category, name, expected):
    assert Isochrone(pathlib.Path('nowhere') / category / name).formatting is expected


# What the backend reads of each system (IsochroneBase::readfile of Galaxia-ananke, src/IsochroneBase.cpp): the column of the age and
# whether it is linear (the backend takes the log of it), then those of the initial mass, the current mass, the luminosity, the
# temperature and the gravity. The backend recognises the systems made in Python by the descriptor that starts with Python (the package
# by their category), the ones of this list by their name, and goes through the same branch as the systems made in Python for the others.
BACKEND_READS = [
    ('py_custom', 'Test', (0, False, 1, 2, 3, 4, 5)),
    *[('padova', name, (0, False, 1, 2, 3, 4, 5)) for name in ('APO__SDSS', 'DCMC', 'Generic__Johnson_UBVRIJHK', 'Generic__Stroemgren')],
    ('padova', 'WFIRST+HST__WFC3', (1, False, 2, 3, 4, 5, 6)),
    *[('padova', name, (1, True, 2, 3, 4, 5, 6)) for name in ('WFIRST', 'LSST', 'GAIA__0', 'GAIA__DR2', 'HST__WFC3')],
    *[('padova', name, (2, False, 3, 5, 6, 7, 8)) for name in ('CTIO__DECam', 'LSST_DP0', 'Roman', 'Euclid', 'JWST__NIRCam', 'JWST__MIRI')],
]


@pytest.mark.parametrize("category, name, reads", BACKEND_READS, ids=[f"{category}/{name}" for category, name, _ in BACKEND_READS])
def test_the_package_lays_each_system_out_as_the_backend_reads_it(category, name, reads):
    formatting = Isochrone(pathlib.Path('nowhere') / category / name).formatting
    assert (formatting.i_age, formatting.lin_age, formatting.i_mini, formatting.i_mact, formatting.i_lum, formatting.i_teff,
            formatting.i_grav) == reads


# -- the tables ---------------------------------------------------------------------

STARS = [(9.0, 1.0, 0.99, 0.0, 3.76, 4.4, 0, 5.0, 5.5),
         (9.0, 2.0, 1.9, 1.0, 3.9, 4.0, 0, 3.0, 3.2),
         (10.0, 1.0, 0.98, 0.1, 3.7, 4.3, 0, 5.2, 5.7),
         (10.0, 2.0, 0.6, -9.999, 4.4, 8.1, 9, 16.0, 16.5)]


def test_the_tables_are_for_each_metallicity_in_log_and_then_each_age(tmp_path):
    iso = Isochrone(isochrone_directory(tmp_path, 'X', {0.0152: STARS, 0.152: STARS}))
    assert list(iso.qtables_dictionary) == pytest.approx([0.0, 1.0])
    assert [list(ages) for ages in iso.qtables_dictionary.values()] == [[9.0, 10.0], [9.0, 10.0]]
    assert [list(ages) for ages in iso.qtables_unique_ages_dictionary.values()] == [[9.0, 10.0], [9.0, 10.0]]


def test_the_tables_of_an_isochrone_made_in_python_hold_its_magnitudes(isochrones_path):
    iso = Isochrone('mine', custom_data())
    table = next(iter(iso.qtables_dictionary.values()))[9.0]
    assert np.asarray(table['G'].value).tolist() == [5.0, 3.0] and np.asarray(table['BP'].value).tolist() == [5.5, 3.2]


# the last star of the first age is a white dwarf: it is carried to the last age, with its magnitudes
WHITE_DWARFS = [(9.0, 1.0, 0.99, 0.0, 3.76, 4.4, 0, 5.0, 5.5), (9.0, 2.0, 1.9, 1.0, 3.9, 4.0, 0, 3.0, 3.2),
                (9.0, 3.0, 0.6, -9.999, 4.5, 8.0, 9, 15.0, 15.5),
                (10.0, 1.0, 0.98, 0.1, 3.7, 4.3, 0, 5.2, 5.7), (10.0, 2.0, 0.6, -9.999, 4.4, 8.1, 9, 16.0, 16.5)]
# the first ages end with stars that are more massive than any of the next age: they are dead there, with 999 for magnitudes, and
# the last of them is carried to the last age
BLACK_HOLES = [(9.0, 1.0, 0.99, 0.0, 3.76, 4.4, 0, 5.0, 5.5), (9.0, 2.0, 1.9, 1.0, 3.9, 4.0, 0, 3.0, 3.2),
               (9.0, 3.0, 2.8, 1.5, 4.0, 3.8, 3, 2.0, 2.2), (9.0, 4.0, 3.7, 1.6, 4.1, 3.7, 3, 1.0, 1.2),
               (10.0, 1.0, 0.98, 0.2, 3.7, 4.3, 0, 5.2, 5.7), (10.0, 2.0, 1.9, 0.9, 3.8, 4.1, 0, 3.1, 3.3),
               (10.0, 3.0, 2.7, 1.2, 3.9, 3.9, 3, 2.1, 2.3),
               (11.0, 1.0, 0.97, 0.3, 3.7, 4.2, 0, 5.4, 5.9), (11.0, 2.0, 0.6, -9.999, 4.4, 8.1, 9, 16.0, 16.5)]


def consolidated(tmp_path, stars):
    return Isochrone(isochrone_directory(tmp_path, 'X', {0.0152: stars})).consolidated_points_and_values


def test_the_points_and_values_of_isochrones_that_end_with_white_dwarfs(tmp_path):
    points, values = consolidated(tmp_path, WHITE_DWARFS)
    np.testing.assert_allclose(points, [[0.0, 9.0, 1.0], [0.0, 9.0, 2.0], [0.0, 9.0, 3.0], [0.0, 10.0, 3.0],
                                        [0.0, 10.0, 1.0], [0.0, 10.0, 2.0]])
    np.testing.assert_allclose(values, [[5.0, 5.5], [3.0, 3.2], [15.0, 15.5], [15.0, 15.5], [5.2, 5.7], [16.0, 16.5]])


def test_the_points_and_values_of_ages_that_end_with_stars_that_are_not_white_dwarfs(tmp_path):
    points, values = consolidated(tmp_path, BLACK_HOLES)
    np.testing.assert_allclose(points, [[0.0, 9.0, 1.0], [0.0, 9.0, 2.0], [0.0, 9.0, 3.0], [0.0, 9.0, 4.0],
                                        [0.0, 10.0, 4.0], [0.0, 11.0, 4.0],
                                        [0.0, 10.0, 1.0], [0.0, 10.0, 2.0], [0.0, 10.0, 3.0],
                                        [0.0, 11.0, 3.0], [0.0, 11.0, 3.0],
                                        [0.0, 11.0, 1.0], [0.0, 11.0, 2.0]])
    np.testing.assert_allclose(values, [[5.0, 5.5], [3.0, 3.2], [2.0, 2.2], [1.0, 1.2], [999.0, 999.0], [999.0, 999.0],
                                        [5.2, 5.7], [3.1, 3.3], [2.1, 2.3], [999.0, 999.0], [999.0, 999.0],
                                        [5.4, 5.9], [16.0, 16.5]])


def test_an_oldest_age_that_does_not_end_with_a_white_dwarf_is_refused_with_what_the_isochrone_needs(tmp_path):
    # Completing the (mass, age) space takes the stars that are dead at an age for white dwarfs below a mass and for black holes above it,
    # which needs the oldest age to be one where the most massive star of the track has died as a white dwarf.
    with pytest.raises(ValueError, match="oldest age of the isochrone .* at metallicity .* does not end with a white dwarf"):
        consolidated(tmp_path, [(9.0, 1.0, 0.99, 0.0, 3.76, 4.4, 0, 5.0, 5.5), (9.0, 2.0, 1.9, 1.0, 3.9, 4.0, 0, 3.0, 3.2)])


def test_points_and_values_have_the_dimensions_of_the_space_and_of_the_magnitudes(tmp_path):
    points, values = consolidated(tmp_path, WHITE_DWARFS)
    assert points.shape == (6, 3) and values.shape == (6, 2)


if __name__ == '__main__':
    pass
