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
Tests of PhotoSystem, which is what a survey uses of a photometric system: its names, the check of the colors of a diagram, the
keys of its magnitudes in the output, its effective wavelengths and zero points (from the saved answers of the SVO filter profile
service, not from the service), and the isochrone tracks nearest to given ages and metallicities. Its isochrone is a stand-in with
the little that a PhotoSystem asks of one.
"""
import types

import numpy as np
import pytest
from astropy import units
from astropy.table import Table

from ..galaxia_ananke.photometry.PhotoSystem import PhotoSystem
from .utils import replay_svo


def isochrone(category='padova', name='GAIA__DR2', mag_names=('G', 'Gbp', 'Grp'), **more):
    return types.SimpleNamespace(category=category, name=name, mag_names=list(mag_names), **more)


def system(**kwargs):
    return PhotoSystem(isochrone(**kwargs))


# -- what a photometric system is called --------------------------------------------

def test_a_photometric_system_has_the_category_and_the_name_of_its_isochrone():
    stand_in = isochrone()
    photo_system = PhotoSystem(stand_in)
    assert (photo_system.category, photo_system.name) == ('padova', 'GAIA__DR2')
    assert photo_system.categ_and_name == ['padova', 'GAIA__DR2'] and photo_system.key == 'padova/GAIA__DR2'
    assert photo_system._isochrone is stand_in
    assert repr(photo_system) == 'PhotoSystem(category=padova, name=GAIA__DR2)'


def test_the_magnitudes_are_those_of_the_isochrone():
    assert system(mag_names=['G', 'Gbp', 'Grp']).mag_names == ['G', 'Gbp', 'Grp']


def test_the_keys_of_the_magnitudes_in_the_output_are_the_name_of_the_system_and_the_magnitude_in_lower_case():
    photo_system = system()
    assert photo_system.to_export_keys == ['gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp']
    assert photo_system.mag_names_to_export_keys_mapping == {'G': 'gaia__dr2_g', 'Gbp': 'gaia__dr2_gbp', 'Grp': 'gaia__dr2_grp'}


# -- the colors of a diagram --------------------------------------------------------

def test_what_the_colors_of_a_dictionary_must_have():
    assert PhotoSystem._required_cmd_magnames_dictkeys == {'magnitude', 'color_minuend', 'color_subtrahend'}


@pytest.mark.parametrize("cmd", ['G,Gbp-Grp', 'Gbp,G-Grp', 'G Gbp Grp'])
def test_colors_that_the_system_has_are_given_back_as_they_were_given(cmd):
    assert system().check_cmd_magnames(cmd) == cmd                    # not rewritten: 'G Gbp Grp' stays what it is


def test_the_colors_of_a_dictionary_are_written_as_the_magnitude_and_the_color():
    cmd = {'magnitude': 'Grp', 'color_minuend': 'Gbp', 'color_subtrahend': 'G'}
    assert system().check_cmd_magnames(cmd) == 'Grp,Gbp-G'


def test_colors_that_the_system_does_not_have_are_refused_with_those_it_has():
    with pytest.raises(ValueError, match=r"CMD magnitude names 'G,Gbp-X' don't match isochrone names: choose among \['G', 'Gbp', 'Grp'\]"):
        system().check_cmd_magnames('G,Gbp-X')


def test_the_colors_of_a_dictionary_that_the_system_does_not_have_are_refused():
    with pytest.raises(ValueError, match="don't match isochrone names"):
        system().check_cmd_magnames({'magnitude': 'G', 'color_minuend': 'Gbp', 'color_subtrahend': 'X'})


@pytest.mark.parametrize("cmd, error", [
    ({'magnitude': 'G', 'color_minuend': 'Gbp'}, "Given cmd_magnames dict covers wrong set of keys: misses {'color_subtrahend'}"),
    ({'magnitude': 'G', 'color_minuend': 'Gbp', 'color_subtrahend': 'Grp', 'other': 'G'}, "misincludes {'other'}"),
])
def test_a_dictionary_of_colors_that_does_not_have_the_right_keys_is_refused(cmd, error):
    with pytest.raises(ValueError, match=error):
        system().check_cmd_magnames(cmd)


@pytest.mark.parametrize("cmd", [['G', 'Gbp', 'Grp'], 5, None])
def test_colors_that_are_neither_a_string_nor_a_dictionary_are_refused(cmd):
    with pytest.raises(ValueError, match="cmd_magnames should be either a string or a dict: given"):
        system().check_cmd_magnames(cmd)


# -- the wavelengths and the zero points --------------------------------------------

GAIA_WAVELENGTHS = [0.58295799337521, 0.50176937808197, 0.7587213835951]            # G, Gbp, Grp: the saved table in angstrom, in micron
GAIA_ZERO_POINTS = [3243.848262579, 3438.2895767386, 2566.0776344001]                # in Jy


def test_the_effective_wavelengths_are_those_of_the_filters_in_micron(monkeypatch):
    replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    wavelengths = system().effective_wavelengths
    assert wavelengths.unit == units.micron
    assert wavelengths.value.tolist() == pytest.approx(GAIA_WAVELENGTHS, rel=1e-9)


def test_the_zero_points_are_those_of_the_filters_in_jansky(monkeypatch):
    replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    zero_points = system().zeropoints
    assert zero_points.unit == units.Jy
    assert zero_points.value.tolist() == pytest.approx(GAIA_ZERO_POINTS, rel=1e-9)


def test_the_wavelengths_and_the_zero_points_are_in_the_order_of_the_magnitudes(monkeypatch):
    replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    photo_system = system(mag_names=['Grp', 'G', 'Gbp'])
    assert photo_system.effective_wavelengths.value.tolist() == pytest.approx([GAIA_WAVELENGTHS[2], GAIA_WAVELENGTHS[0], GAIA_WAVELENGTHS[1]], rel=1e-9)
    assert photo_system.zeropoints.value.tolist() == pytest.approx([GAIA_ZERO_POINTS[2], GAIA_ZERO_POINTS[0], GAIA_ZERO_POINTS[1]], rel=1e-9)


def test_the_service_is_asked_once_for_the_wavelengths_and_the_zero_points(monkeypatch):
    asked = replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    photo_system = system()
    photo_system.effective_wavelengths, photo_system.zeropoints, photo_system.effective_wavelengths
    assert asked == ['GAIA']


def two_instruments_table():
    """What the service says of a facility TEST with two instruments, A and B, that each have a filter x."""
    return Table({'filterID': ['TEST/A.x', 'TEST/B.x'], 'WavelengthEff': [1000.0, 3000.0] * units.AA, 'ZeroPoint': [3000.0, 1000.0] * units.Jy})


def test_the_wavelengths_of_a_system_of_several_instruments_are_the_mean_of_those_of_its_instruments(monkeypatch):
    asked = replay_svo(monkeypatch, TEST=two_instruments_table())
    wavelengths = system(name='TEST__A+TEST__B', mag_names=['x']).effective_wavelengths
    assert wavelengths.unit == units.micron and wavelengths.value.tolist() == pytest.approx([0.2])     # (1000 + 3000) / 2 angstrom
    assert asked == ['TEST', 'TEST']                                                                     # one question for each instrument


def test_the_zero_points_of_a_system_of_several_instruments_are_the_mean_of_those_of_its_instruments(monkeypatch):
    replay_svo(monkeypatch, TEST=two_instruments_table())
    zero_points = system(name='TEST__A+TEST__B', mag_names=['x']).zeropoints
    assert zero_points.unit == units.Jy and zero_points.value.tolist() == pytest.approx([2000.0])


def test_a_system_that_is_a_facility_only_takes_the_filters_of_all_the_instruments_of_the_facility(monkeypatch):
    table = Table({'filterID': ['TEST/A.x', 'TEST/B.y'], 'WavelengthEff': [1000.0, 3000.0] * units.AA, 'ZeroPoint': [3000.0, 1000.0] * units.Jy})
    asked = replay_svo(monkeypatch, TEST=table)
    photo_system = system(name='TEST', mag_names=['x', 'y'])
    assert photo_system.effective_wavelengths.value.tolist() == pytest.approx([0.1, 0.3])           # in micron: 1000 and 3000 angstrom
    assert photo_system.zeropoints.value.tolist() == pytest.approx([3000.0, 1000.0])
    assert asked == ['TEST']


@pytest.mark.parametrize("name", ['LSST', 'LSST_DP0'])
def test_there_is_a_wavelength_for_each_magnitude_of_lsst(monkeypatch, name):
    replay_svo(monkeypatch, LSST='svo_raw_LSST.ecsv')
    assert len(system(name=name, mag_names=['u', 'g', 'r', 'i', 'z', 'y']).effective_wavelengths) == 6


# -- the tracks nearest to given ages and metallicities -----------------------------

LOG_METALLICITIES = {-1.0: [8.0, 9.0, 10.0], 0.0: [8.0, 9.0, 10.0], 0.5: [8.5, 9.5]}      # the log ages that each metallicity has


def tracks_isochrone():
    """An isochrone whose tables are the names of the metallicity and the age that they are for."""
    return dict(qtables_unique_ages_dictionary={z: list(ages) for z, ages in LOG_METALLICITIES.items()},
                qtables_dictionary={z: {age: f"{z}/{age}" for age in ages} for z, ages in LOG_METALLICITIES.items()})


def test_each_particle_gets_the_table_of_the_nearest_metallicity_and_then_the_nearest_age():
    photo_system = PhotoSystem(isochrone(**tracks_isochrone()))
    tables = photo_system.get_nearest_isochrone_track_for_age_and_metallicity(np.array([9.1, 8.2, 9.8, 8.6]), np.array([0.1, -0.9, 0.4, 0.45]))
    # the nearest metallicities are 0.0, -1.0, 0.5 and 0.5; the nearest ages for them are 9.0, 8.0, 9.5 and 8.5
    assert tables == ['0.0/9.0', '-1.0/8.0', '0.5/9.5', '0.5/8.5']


def test_the_nearest_is_the_one_across_the_middle_between_two_values():
    photo_system = PhotoSystem(isochrone(**tracks_isochrone()))
    first, second = photo_system.get_nearest_isochrone_track_for_age_and_metallicity(np.array([9.0, 9.0]), np.array([-0.4, -0.6]))
    assert (first, second) == ('0.0/9.0', '-1.0/9.0')                  # -0.5 is the middle between -1 and 0


def test_particles_that_all_have_the_same_nearest_metallicity_each_get_a_table():
    photo_system = PhotoSystem(isochrone(**tracks_isochrone()))
    tables = photo_system.get_nearest_isochrone_track_for_age_and_metallicity(np.array([9.1, 8.2]), np.array([0.1, 0.05]))
    assert tables == ['0.0/9.0', '0.0/8.0']


def test_a_single_particle_gets_its_table():
    photo_system = PhotoSystem(isochrone(**tracks_isochrone()))
    assert photo_system.get_nearest_isochrone_track_for_age_and_metallicity(np.array([9.1]), np.array([0.1])) == ['0.0/9.0']


if __name__ == '__main__':
    pass
