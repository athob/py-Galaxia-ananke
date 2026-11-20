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
from astropy.table import MaskedColumn, Table

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


def test_the_plus_of_a_name_that_joins_facilities_stays_in_the_keys_of_the_magnitudes():
    photo_system = system(name='2MASS+Spitzer+WISE', mag_names=['J', 'Ks', '24mu', 'W1'])
    assert photo_system.to_export_keys == ['2mass+spitzer+wise_j', '2mass+spitzer+wise_ks', '2mass+spitzer+wise_24mu', '2mass+spitzer+wise_w1']
    assert photo_system.mag_names_to_export_keys_mapping['Ks'] == '2mass+spitzer+wise_ks'


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


# -- a facility that the service has no filter of -----------------------------------

NO_FILTER = 'No filter found for requested Facilty (and Instrument)'


def test_the_wavelengths_and_the_zero_points_of_a_facility_that_the_service_has_no_filter_of_are_nan(monkeypatch):
    replay_svo(monkeypatch, Generic=IndexError(NO_FILTER))
    photo_system = system(name='Generic__Johnson_UBVRIJHK', mag_names=['U', 'B'])
    with pytest.warns(UserWarning, match="no filter for the facility 'Generic'") as warned:
        wavelengths, zero_points = photo_system.effective_wavelengths, photo_system.zeropoints
    assert wavelengths.unit == units.micron and np.isnan(wavelengths.value).all() and len(wavelengths) == 2
    assert np.isnan(zero_points.value).all() and len(zero_points) == 2
    assert len(warned) == 1                                                    # once, not for each of the two quantities


def test_the_instruments_that_the_service_has_filters_of_keep_their_values_beside_one_that_it_has_not(monkeypatch):
    replay_svo(monkeypatch, TEST=two_instruments_table(), NOTHING=IndexError(NO_FILTER))
    photo_system = system(name='TEST__A+NOTHING', mag_names=['x'])
    with pytest.warns(UserWarning, match="no filter for the facility 'NOTHING'") as warned:
        wavelengths, zero_points = photo_system.effective_wavelengths, photo_system.zeropoints
    # the 1000 angstrom of the instrument A: the other has no value
    assert wavelengths.value.tolist() == pytest.approx([0.1])
    alone = system(name='TEST__A', mag_names=['x'])
    assert not np.isnan(zero_points.value).any()
    assert zero_points.value.tolist() == pytest.approx(alone.zeropoints.value.tolist())
    assert len(warned) == 1


def test_the_offsets_of_a_facility_without_filters_are_nan_with_a_warning_for_it_and_one_for_the_filters(monkeypatch):
    replay_svo(monkeypatch, Generic=IndexError(NO_FILTER))
    photo_system = system(name='Generic__Johnson_UBVRIJHK', mag_names=['U', 'B'])
    with pytest.warns(UserWarning) as warned:
        offsets = photo_system.vega_to_ab_offsets
    messages = [str(record.message) for record in warned]
    assert any("no filter for the facility 'Generic'" in message for message in messages)
    assert any("no zero point of the Vega system for the filters ['U', 'B']" in message for message in messages)
    assert list(offsets) == ['generic__johnson_ubvrijhk_u', 'generic__johnson_ubvrijhk_b']
    assert np.isnan(list(offsets.values())).all()


# -- the magnitude system, and the offsets between Vega and AB ----------------------

GAIA_VEGA_TO_AB = [0.1224, 0.0592, 0.3769]          # G, Gbp, Grp: 2.5 log10(3631 Jy / the zero points above), in magnitude
# of the filters of the saved table, with their total throughputs
LSST_VEGA_TO_AB = {'lsst_u': 0.6269, 'lsst_g': -0.1025, 'lsst_r': 0.1495}


def vega_instruments_table(systems=('Vega', 'Vega')):
    """What the service says of a facility TEST: instruments A and B that have a filter x, and the system of its zero point."""
    return Table({'filterID': ['TEST/A.x', 'TEST/B.x'], 'MagSys': list(systems), 'WavelengthEff': [1000.0, 3000.0] * units.AA,
                  'ZeroPoint': [3000.0, 1000.0] * units.Jy})


def test_the_magnitude_system_is_that_of_the_isochrone():
    assert system(magnitude_system='AB').magnitude_system == 'AB'
    assert system(magnitude_system=None).magnitude_system is None


def test_the_offsets_to_ab_are_those_of_the_vega_zero_points_by_the_keys_of_the_output(monkeypatch):
    replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    offsets = system().vega_to_ab_offsets
    assert list(offsets) == ['gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp']
    assert list(offsets.values()) == pytest.approx(GAIA_VEGA_TO_AB, abs=5e-5)
    by_the_formula = [2.5 * np.log10(3631 / zero_point) for zero_point in GAIA_ZERO_POINTS]
    assert list(offsets.values()) == pytest.approx(by_the_formula, rel=1e-9)


def test_the_offsets_to_ab_of_lsst_are_those_of_its_total_throughputs(monkeypatch):
    replay_svo(monkeypatch, LSST='svo_raw_LSST.ecsv')
    offsets = system(name='LSST', mag_names=['u', 'g', 'r']).vega_to_ab_offsets
    assert offsets == pytest.approx(LSST_VEGA_TO_AB, abs=5e-5)


def test_the_offsets_are_in_the_order_of_the_magnitudes(monkeypatch):
    replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    offsets = system(mag_names=['Grp', 'G', 'Gbp']).vega_to_ab_offsets
    assert list(offsets) == ['gaia__dr2_grp', 'gaia__dr2_g', 'gaia__dr2_gbp']
    assert list(offsets.values()) == pytest.approx([GAIA_VEGA_TO_AB[2], GAIA_VEGA_TO_AB[0], GAIA_VEGA_TO_AB[1]], abs=5e-5)


def test_the_service_is_asked_once_for_the_offsets_and_for_the_zero_points(monkeypatch):
    asked = replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    photo_system = system()
    photo_system.vega_to_ab_offsets, photo_system.zeropoints, photo_system.vega_to_ab_offsets
    assert asked == ['GAIA']


def test_the_offsets_of_a_system_of_several_instruments_use_the_mean_of_their_vega_zero_points(monkeypatch):
    replay_svo(monkeypatch, TEST=vega_instruments_table())
    offsets = system(name='TEST__A+TEST__B', mag_names=['x']).vega_to_ab_offsets
    assert offsets == pytest.approx({'test__a+test__b_x': 2.5 * np.log10(3631 / 2000)})            # (3000 + 1000) / 2 Jy


def test_a_zero_point_that_is_not_that_of_the_vega_system_is_left_out(monkeypatch):
    replay_svo(monkeypatch, TEST=vega_instruments_table(('AB', 'Vega')))
    offsets = system(name='TEST__A+TEST__B', mag_names=['x']).vega_to_ab_offsets
    assert offsets == pytest.approx({'test__a+test__b_x': 2.5 * np.log10(3631 / 1000)})            # that of B only


def test_a_filter_that_has_no_zero_point_of_the_vega_system_has_nan_and_a_warning(monkeypatch):
    replay_svo(monkeypatch, TEST=vega_instruments_table(('AB', 'AB')))
    photo_system = system(name='TEST__A+TEST__B', mag_names=['x'])
    with pytest.warns(UserWarning, match=r"no zero point of the Vega system for the filters \['x'\] of padova/TEST__A\+"):
        offsets = photo_system.vega_to_ab_offsets
    assert np.isnan(offsets['test__a+test__b_x'])


def test_a_filter_whose_row_is_masked_has_nan_and_the_others_keep_their_offsets(monkeypatch):
    table = Table({'filterID': ['TEST/A.x', 'TEST/A.y'], 'WavelengthEff': [1000.0, 3000.0] * units.AA,
                   'MagSys': MaskedColumn(['Vega', 'Vega'], mask=[False, True]),
                   'ZeroPoint': MaskedColumn([3000.0, 1000.0], mask=[False, True], unit=units.Jy)})
    replay_svo(monkeypatch, TEST=table)
    photo_system = system(name='TEST', mag_names=['x', 'y'])
    with pytest.warns(UserWarning, match=r"for the filters \['y'\]"):
        offsets = photo_system.vega_to_ab_offsets
    assert offsets['test_x'] == pytest.approx(2.5 * np.log10(3631 / 3000)) and np.isnan(offsets['test_y'])


def test_the_offsets_to_the_system_that_the_magnitudes_are_in_are_zero_and_the_service_is_not_asked(monkeypatch):
    asked = replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    nothing = {'gaia__dr2_g': 0.0, 'gaia__dr2_gbp': 0.0, 'gaia__dr2_grp': 0.0}
    assert system(magnitude_system='AB').offsets_to_magnitude_system('ab') == nothing
    assert system(magnitude_system='Vega').offsets_to_magnitude_system('Vega') == nothing
    assert asked == []


def test_the_offsets_from_vega_to_ab_are_the_offsets_to_ab(monkeypatch):
    replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    photo_system = system(magnitude_system='Vega')
    assert photo_system.offsets_to_magnitude_system('AB') == photo_system.vega_to_ab_offsets


def test_the_offsets_from_ab_to_vega_are_minus_the_offsets_to_ab(monkeypatch):
    replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    photo_system = system(magnitude_system='AB')
    to_vega = photo_system.offsets_to_magnitude_system('Vega')
    assert to_vega == {key: -offset for key, offset in photo_system.vega_to_ab_offsets.items()}
    assert list(to_vega.values()) == pytest.approx([-offset for offset in GAIA_VEGA_TO_AB], abs=5e-5)


def test_magnitudes_that_do_not_say_their_system_cannot_be_converted():
    with pytest.raises(ValueError, match="magnitude system of padova/GAIA__DR2 is not stated"):
        system(magnitude_system=None).offsets_to_magnitude_system('AB')


@pytest.mark.parametrize("native, target", [('ST', 'AB'), ('Vega', 'ST'), ('AB', 'ST')])
def test_only_the_conversions_between_vega_and_ab_are_there(native, target):
    with pytest.raises(NotImplementedError, match=f"not from {native} to {target}"):
        system(magnitude_system=native).offsets_to_magnitude_system(target)


@pytest.mark.parametrize("target", ['Johnson', None, ''])
def test_a_system_to_convert_to_that_is_not_one_is_refused(target):
    with pytest.raises(ValueError):
        system(magnitude_system='Vega').offsets_to_magnitude_system(target)


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
