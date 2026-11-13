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
Tests of Formatting, which says how the files of the isochrones of a photometric system are laid out (which column is the age, the
mass, ...) and turns what they hold into one table per age with units. The four layouts that the package knows are pinned, and the
tables are checked on small isochrone files whose numbers are worked out here, outside astropy.
"""
import dataclasses
from operator import itemgetter
import types

import numpy as np
import pytest
from astropy import units

from ..galaxia_ananke.photometry.Formatting import (cmd_2_1_formatting, python_formatting, cmd_3_0_formatting,
                                                    cmd_2_8_formatting, cmd_3_6_formatting)
from ..galaxia_ananke.photometry.IsochroneFile import IsochroneFile
from .utils import write_isochrone_file

G, M_SUN, R_SUN = 6.6743e-11, 1.988409870698051e30, 6.957e8             # SI, as astropy has them (to a part in 10^4 or better)


def radius_in_solar_radii(mass, logg):
    """What the radius of a star is, in solar radii, from its mass in solar masses and its gravity: log10 of cm/s^2."""
    return ((G * mass * M_SUN) / 10 ** (logg - 2)) ** 0.5 / R_SUN


def values(table, name):
    return np.asarray(table[name].value).tolist()


# -- the four layouts ---------------------------------------------------------------

@pytest.mark.parametrize("formatting, layout", [
    # CMD2-label: (python_formatting, (0, 1, 2, 3, 4, 5, 6, False)),                    # what is made in Python, and others
    # CMD2-label: (cmd_3_6_formatting, (2, 3, 5, 6, 7, 8, 9, False)),                      # CTIO__DECam, LSST_DP0, Roman, Euclid, JWST
    # CMD2-label: (cmd_3_0_formatting, (1, 2, 3, 4, 5, 6, 7, True)),                     # WFIRST, LSST, GAIA__0, GAIA__DR2: linear ages
    # CMD2-label: (cmd_2_8_formatting, (1, 2, 3, 4, 5, 6, 7, False)),         # WFIRST+HST__WFC3
    (python_formatting, (0, 1, 2, 3, 4, 5, False)),                       # what is made in Python, and others
    (cmd_3_6_formatting, (2, 3, 5, 6, 7, 8, False)),                         # CTIO__DECam, LSST_DP0, Roman, Euclid, JWST
    (cmd_3_0_formatting, (1, 2, 3, 4, 5, 6, True)),                        # WFIRST, LSST, GAIA__0, GAIA__DR2: linear ages
    (cmd_2_8_formatting, (1, 2, 3, 4, 5, 6, False)),            # WFIRST+HST__WFC3
], ids=["python", "cmd_3_6", "cmd_3_0", "cmd_2_8"])
def test_the_layouts_of_the_isochrone_files_that_the_package_knows(formatting, layout):
    assert dataclasses.astuple(formatting) == layout + (0.0152,)


def test_the_names_of_the_columns_of_a_table_of_isochrones():
    assert (python_formatting._zini, python_formatting._age, python_formatting._mini, python_formatting._mact,
            python_formatting._radi, python_formatting._lum, python_formatting._teff,
            python_formatting._grav) == ('Z_ini', 'Age', 'M_ini', 'M_act', 'Radius', 'Lum', 'T_eff', 'Grav')
    # CMD2-label: assert (python_formatting._zini, python_formatting._age, python_formatting._mini, python_formatting._mact,
    # CMD2-label:         python_formatting._radi, python_formatting._lum, python_formatting._teff, python_formatting._grav,
    # CMD2-label:         python_formatting._label) == ('Z_ini', 'Age', 'M_ini', 'M_act', 'Radius', 'Lum', 'T_eff', 'Grav', 'Label')


def test_every_column_of_a_table_of_isochrones_has_a_unit():
    mapper = python_formatting.units_mapper
    assert set(mapper) == {'Z_ini', 'Age', 'M_ini', 'M_act', 'Radius', 'Lum', 'T_eff', 'Grav'}
    # CMD2-label: assert set(mapper) == {'Z_ini', 'Age', 'M_ini', 'M_act', 'Radius', 'Lum', 'T_eff', 'Grav', 'Label'}
    assert mapper['M_ini'] == mapper['M_act'] == units.Msun and mapper['Radius'] == units.Rsun


def test_importing_the_formatting_defines_the_solar_metallicity_unit_in_astropy():
    # How it is today: astropy.units gets an attribute for the solar metallicity (and, with the label, one for the integers)
    assert units.Zsun is not None
    # CMD2-label: assert units.integer.dtype.kind == 'i'


def test_the_columns_are_mapped_to_the_names_by_their_position():
    assert python_formatting.format_mapper == {0: 'Age', 1: 'M_ini', 2: 'M_act', 3: 'Lum', 4: 'T_eff', 5: 'Grav'}
    assert cmd_3_6_formatting.format_mapper == {2: 'Age', 3: 'M_ini', 5: 'M_act', 6: 'Lum', 7: 'T_eff', 8: 'Grav'}
    # CMD2-label: assert python_formatting.format_mapper == {0: 'Age', 1: 'M_ini', 2: 'M_act', 3: 'Lum', 4: 'T_eff', 5: 'Grav', 6: 'Label'}
    # CMD2-label: assert cmd_3_6_formatting.format_mapper == {2: 'Age', 3: 'M_ini', 5: 'M_act', 6: 'Lum', 7: 'T_eff', 8: 'Grav', 9: 'Label'}


def test_the_columns_that_are_not_magnitudes_are_picked_in_the_order_of_the_mapper():
    names = [f"column{i}" for i in range(12)]
    assert cmd_3_6_formatting.nonmagnitude_itemgetter(names) == ('column2', 'column3', 'column5', 'column6', 'column7', 'column8')
    # CMD2-label: assert cmd_3_6_formatting.nonmagnitude_itemgetter(names) == ('column2', 'column3', 'column5', 'column6', 'column7', 'column8', 'column9')


# -- the metallicity ----------------------------------------------------------------

@pytest.mark.parametrize("metallicity, expected", [(0.0152, 0.0), (0.152, 1.0), (0.00152, -1.0)])
def test_the_metallicity_is_the_log_of_its_ratio_to_that_of_the_sun(metallicity, expected):
    assert python_formatting.metallicity_converter(metallicity) == pytest.approx(expected)


def test_the_systems_of_cmd_2_1_and_2_2_have_the_columns_of_the_python_layout_and_the_solar_metallicity_of_marigo_et_al():
    assert dataclasses.astuple(cmd_2_1_formatting) == (0, 1, 2, 3, 4, 5, False, 0.019)
    assert cmd_2_1_formatting == dataclasses.replace(python_formatting, feh_sun=0.019)
    # CMD2-label: assert dataclasses.astuple(cmd_2_1_formatting) == (0, 1, 2, 3, 4, 5, 6, False, 0.019)


@pytest.mark.parametrize("metallicity, expected", [(0.019, 0.0), (0.19, 1.0), (0.0019, -1.0)])
def test_the_metallicity_of_cmd_2_1_and_2_2_is_relative_to_that_of_the_sun_of_marigo_et_al(metallicity, expected):
    assert cmd_2_1_formatting.metallicity_converter(metallicity) == pytest.approx(expected)


def test_the_metallicity_of_several_isochrones_can_be_converted_together():
    assert python_formatting.metallicity_converter(np.array([0.0152, 0.152])).tolist() == pytest.approx([0.0, 1.0])


# -- a table for each age -----------------------------------------------------------

# age, initial mass, current mass, log luminosity, log temperature, log gravity, label, and two magnitudes
STARS = [(9.0, 1.0, 0.99, 0.0, 3.76, 4.4, 0, 5.0, 5.5),
         (9.0, 0.5, 0.5, -1.0, 3.5, 4.7, 0, 8.0, 8.5),
         (9.0, 2.0, 1.9, 1.0, 3.9, 4.0, 1, 3.0, 3.2),
         (10.0, 0.8, 0.79, -0.3, 3.7, 4.5, 0, 6.0, 6.4),
         (10.0, 0.4, 0.4, -1.5, 3.4, 4.8, 0, 9.0, 9.4)]
NAMES = ['age', 'mini', 'mact', 'lum', 'teff', 'grav', 'label', 'G', 'BP']


def isochrone_file(tmp_path, metallicity=0.0152, stars=STARS, first_magnitude=7, extra_columns=(), names=NAMES):
    path = write_isochrone_file(tmp_path / f"output_{metallicity}.dat", list(extra_columns) + names,
                                [[0.0152] * len(extra_columns) + list(star) for star in stars])
    isochrone = types.SimpleNamespace(magnitude_itemgetter=itemgetter(first_magnitude, first_magnitude + 1), mag_names=['G', 'BP'])
    return IsochroneFile(path, isochrone=isochrone)


def test_there_is_a_table_for_each_age(tmp_path):
    tables = python_formatting.qtable_per_age_from_isochronefile(isochrone_file(tmp_path))
    assert list(tables) == [9.0, 10.0]
    assert [len(tables[age]) for age in tables] == [3, 2]


def test_a_table_has_the_columns_in_the_same_order_with_the_metallicity_first_and_the_radius_added(tmp_path):
    table = python_formatting.qtable_per_age_from_isochronefile(isochrone_file(tmp_path))[9.0]
    assert table.colnames == ['Z_ini', 'Age', 'M_ini', 'M_act', 'Lum', 'Radius', 'T_eff', 'Grav', 'G', 'BP']
    # CMD2-label: assert table.colnames == ['Z_ini', 'Age', 'M_ini', 'M_act', 'Lum', 'Radius', 'T_eff', 'Grav', 'Label', 'G', 'BP']


def test_a_column_that_the_layout_does_not_take_is_left_out_of_the_tables(tmp_path):
    # the file has a column between the gravity and the magnitudes (the label of CMD 3.x, the bolometric magnitude of CMD 2.x)
    table = python_formatting.qtable_per_age_from_isochronefile(isochrone_file(tmp_path))[9.0]
    assert 'label' not in table.colnames and len(table.colnames) == 10


def test_a_table_is_made_of_files_whose_magnitudes_follow_the_gravity_as_in_the_ones_made_in_python(tmp_path):
    stars, names = [star[:6] + star[7:] for star in STARS], NAMES[:6] + NAMES[7:]                # no column between them
    table = python_formatting.qtable_per_age_from_isochronefile(isochrone_file(tmp_path, stars=stars, names=names, first_magnitude=6))[9.0]
    assert table.colnames == ['Z_ini', 'Age', 'M_ini', 'M_act', 'Lum', 'Radius', 'T_eff', 'Grav', 'G', 'BP']
    assert values(table, 'G') == [8.0, 5.0, 3.0] and values(table, 'BP') == [8.5, 5.5, 3.2]


def test_the_stars_of_an_age_are_in_the_order_of_their_initial_mass(tmp_path):
    table = python_formatting.qtable_per_age_from_isochronefile(isochrone_file(tmp_path))[9.0]
    assert values(table, 'M_ini') == [0.5, 1.0, 2.0]
    assert values(table, 'M_act') == [0.5, 0.99, 1.9]                          # the rows stay together
    assert values(table, 'Lum') == [-1.0, 0.0, 1.0] and values(table, 'T_eff') == [3.5, 3.76, 3.9]
    assert values(table, 'Grav') == [4.7, 4.4, 4.0]
    # CMD2-label: assert values(table, 'Label') == [0, 0, 1]
    assert values(table, 'G') == [8.0, 5.0, 3.0] and values(table, 'BP') == [8.5, 5.5, 3.2]
    assert values(table, 'Age') == [9.0, 9.0, 9.0]


def test_every_star_has_the_metallicity_of_the_file_in_log_relative_to_the_sun(tmp_path):
    tables = python_formatting.qtable_per_age_from_isochronefile(isochrone_file(tmp_path, metallicity=0.152))
    assert values(tables[9.0], 'Z_ini') == pytest.approx([1.0, 1.0, 1.0])


def test_every_star_of_cmd_2_1_and_2_2_has_the_metallicity_of_the_file_in_log_relative_to_the_sun_of_marigo_et_al(tmp_path):
    tables = cmd_2_1_formatting.qtable_per_age_from_isochronefile(isochrone_file(tmp_path, metallicity=0.19))
    assert values(tables[9.0], 'Z_ini') == pytest.approx([1.0, 1.0, 1.0])


def test_the_radius_of_a_star_is_made_from_its_current_mass_and_its_gravity(tmp_path):
    table = python_formatting.qtable_per_age_from_isochronefile(isochrone_file(tmp_path))[9.0]
    expected = [radius_in_solar_radii(mass, logg) for mass, logg in ((0.5, 4.7), (0.99, 4.4), (1.9, 4.0))]
    assert values(table, 'Radius') == pytest.approx(expected, rel=1e-3)
    assert table['Radius'].unit == units.Rsun


def test_the_sun_has_the_radius_of_the_sun():
    assert radius_in_solar_radii(1.0, 4.438) == pytest.approx(1.0, abs=2e-4)          # the check of the check


def test_a_layout_with_the_age_in_years_gives_tables_for_the_log_of_the_ages(tmp_path):
    stars = [(1e9, 1.0, 0.99, 0.0, 3.76, 4.4, 0, 5.0, 5.5), (2e9, 0.8, 0.79, -0.3, 3.7, 4.5, 0, 6.0, 6.4)]
    iso_file = isochrone_file(tmp_path, stars=stars, first_magnitude=8, extra_columns=['Z'])
    tables = cmd_3_0_formatting.qtable_per_age_from_isochronefile(iso_file)
    assert list(tables) == pytest.approx([9.0, np.log10(2e9)])
    assert values(tables[list(tables)[0]], 'Age') == pytest.approx([9.0])


def test_the_same_layout_with_the_age_in_log_keeps_the_ages_as_they_are(tmp_path):
    stars = [(9.0, 1.0, 0.99, 0.0, 3.76, 4.4, 0, 5.0, 5.5), (9.5, 0.8, 0.79, -0.3, 3.7, 4.5, 0, 6.0, 6.4)]
    iso_file = isochrone_file(tmp_path, stars=stars, first_magnitude=8, extra_columns=['Z'])
    assert list(cmd_2_8_formatting.qtable_per_age_from_isochronefile(iso_file)) == [9.0, 9.5]


if __name__ == '__main__':
    pass
