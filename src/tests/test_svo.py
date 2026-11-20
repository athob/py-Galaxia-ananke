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
Tests of the drivers that ask the SVO filter profile service about the filters of a photometric system: how the name of the system
becomes the names of the instruments to ask about (InterfaceSvoFpsDriver), and the table of filters of an instrument, in the order
of the magnitudes of the system (SvoFpsDriver). The service is not used: its answers are tables that were saved (tests/data), the
Gaia one and the LSST one as astroquery returned them.
"""
import types
import warnings

import numpy as np
import pytest
from astropy.table import Table

from ..galaxia_ananke.photometry.InterfaceSvoFpsDriver import InterfaceSvoFpsDriver
from ..galaxia_ananke.photometry.SvoFpsDriver import SvoFpsDriver
from .utils import replay_svo


def interface(name, mag_names=('G',)):
    """The interface of a photometric system that has this name and these magnitudes, and nothing else."""
    return InterfaceSvoFpsDriver(types.SimpleNamespace(name=name, mag_names=list(mag_names)))


# -- from the name of a system to the names of its instruments ----------------------

def test_the_names_that_the_service_knows_differently():
    assert InterfaceSvoFpsDriver._special_cases == {
        'DCMC': 'DENIS+2MASS',
        'Euclid': 'Euclid__NISP0+Euclid__VIS',
        'GAIA__0': 'GAIA__GAIA0',
        'GAIA__DR2': 'GAIA__GAIA2r',
        'Generic__Johnson_UBVRIJHK': 'Generic__Johnson_UBVRIJHKL',
        'Generic__Stroemgren': 'Generic__Stromgren',
        'LSST': 'LSST__LSST',
        'LSST_DP0': 'LSST__LSST',
        'Pan_STARRS1': 'PAN-STARRS',
        'HST__WFC3': 'HST__WFC3_UVIS1+HST__WFC3_UVIS2+HST__WFC3_IR',
    }


def test_the_separators_of_the_names():
    assert InterfaceSvoFpsDriver._plus_in_isochrone_name == '+'                 # between the instruments of one system
    assert SvoFpsDriver._slash_in_isochrone_name == '__'                        # between the facility and the instrument


@pytest.mark.parametrize("name, instruments", [
    ('GAIA__DR2', ['GAIA__GAIA2r']),
    ('GAIA__0', ['GAIA__GAIA0']),
    ('Euclid', ['Euclid__NISP0', 'Euclid__VIS']),
    ('HST__WFC3', ['HST__WFC3_UVIS1', 'HST__WFC3_UVIS2', 'HST__WFC3_IR']),
    ('GAIA__DR2+HST__WFC3', ['GAIA__GAIA2r', 'HST__WFC3_UVIS1', 'HST__WFC3_UVIS2', 'HST__WFC3_IR']),   # each part has its own special case (a made-up name: no system joins them)
    ('DCMC', ['DENIS', '2MASS']),
    ('2MASS+Spitzer+WISE', ['2MASS', 'Spitzer', 'WISE']),                    # a name that joins facilities is asked as each of them
    ('LSST_DP0', ['LSST__LSST']),
    ('LSST', ['LSST__LSST']),
    ('Pan_STARRS1', ['PAN-STARRS']),
    ('Generic__Johnson_UBVRIJHK', ['Generic__Johnson_UBVRIJHKL']),
    ('Generic__Stroemgren', ['Generic__Stromgren']),
    ('CTIO__DECam', ['CTIO__DECam']),                                        # the names that are not special are the same
    ('JWST__NIRCam', ['JWST__NIRCam']),
])
def test_the_instruments_of_a_photometric_system(name, instruments):
    assert interface(name).list_of_instruments == instruments


def test_a_name_that_is_not_special_is_left_as_it_is():
    assert InterfaceSvoFpsDriver._special_cases_formatter('Something__Else') == 'Something__Else'
    assert InterfaceSvoFpsDriver._special_cases_formatter('GAIA__DR2') == 'GAIA__GAIA2r'


def test_there_is_a_driver_for_each_instrument_that_asks_for_it():
    iface = interface('GAIA__DR2+HST__WFC3')
    assert len(iface.svofpsdrivers) == 4 and all(isinstance(driver, SvoFpsDriver) for driver in iface.svofpsdrivers)
    assert [driver._instrument_name for driver in iface.svofpsdrivers] == iface.list_of_instruments
    assert all(driver._interface is iface for driver in iface.svofpsdrivers)


def test_the_interface_says_the_name_and_the_magnitudes_of_its_system():
    iface = interface('GAIA__DR2', ['G', 'Gbp', 'Grp'])
    assert iface.name == 'GAIA__DR2' and iface.mag_names == ['G', 'Gbp', 'Grp']


@pytest.mark.parametrize("instrument, facility, name", [
    ('GAIA__GAIA2r', 'GAIA', 'GAIA2r'),
    ('HST__WFC3_UVIS1', 'HST', 'WFC3_UVIS1'),
    ('Generic__Johnson_UBVRIJHKL', 'Generic', 'Johnson_UBVRIJHKL'),
    ('LSST__LSST', 'LSST', 'LSST'),
    ('PAN-STARRS', 'PAN-STARRS', None),                                       # a name without a separator is a facility only
])
def test_an_instrument_has_a_facility_and_perhaps_a_name_of_its_own(instrument, facility, name):
    driver = interface(instrument).svofpsdrivers[0]
    assert (driver.facility, driver.instrument) == (facility, name)


# -- the table of filters of an instrument ------------------------------------------

GAIA_EFFECTIVE_WAVELENGTHS = {'G': 5829.5799337521, 'Gbp': 5017.6937808197, 'Grp': 7587.213835951}        # in angstrom, in the saved table
GAIA_ZERO_POINTS = {'G': 3243.848262579, 'Gbp': 3438.2895767386, 'Grp': 2566.0776344001}                  # in Jy


def gaia_table(monkeypatch, mag_names=('G', 'Gbp', 'Grp')):
    asked = replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    return interface('GAIA__DR2', mag_names).svofpsdrivers[0].svo_filter_list, asked


def test_the_filters_of_an_instrument_are_those_of_it_alone_in_the_order_of_the_magnitudes(monkeypatch):
    table, _ = gaia_table(monkeypatch)
    assert [str(band) for band in table['Band']] == ['G', 'Gbp', 'Grp']
    assert [str(instrument) for instrument in table['Instrument']] == ['GAIA2r'] * 3          # not the other instruments of Gaia
    assert table['WavelengthEff'].tolist() == pytest.approx([GAIA_EFFECTIVE_WAVELENGTHS[band] for band in ('G', 'Gbp', 'Grp')], rel=1e-9)
    assert table['ZeroPoint'].tolist() == pytest.approx([GAIA_ZERO_POINTS[band] for band in ('G', 'Gbp', 'Grp')], rel=1e-9)


def test_the_order_of_the_filters_follows_that_of_the_magnitudes_of_the_system(monkeypatch):
    table, _ = gaia_table(monkeypatch, ['Grp', 'G', 'Gbp'])
    assert [str(band) for band in table['Band']] == ['Grp', 'G', 'Gbp']
    assert table['WavelengthEff'].tolist() == pytest.approx([GAIA_EFFECTIVE_WAVELENGTHS[band] for band in ('Grp', 'G', 'Gbp')], rel=1e-9)


def test_the_service_is_asked_about_the_facility_and_not_about_the_instrument(monkeypatch):
    _, asked = gaia_table(monkeypatch)
    assert asked == ['GAIA']                                          # not GAIA2r: the instruments are told apart in what comes back


def test_the_table_of_an_instrument_is_made_once(monkeypatch):
    asked = replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    driver = interface('GAIA__DR2', ['G', 'Gbp', 'Grp']).svofpsdrivers[0]
    assert driver.svo_filter_list is driver.svo_filter_list
    assert asked == ['GAIA']


def test_a_magnitude_that_is_a_filter_in_other_letters_takes_the_values_of_that_filter(monkeypatch):
    table, _ = gaia_table(monkeypatch, ['g', 'Gbp', 'Grp'])
    assert [str(band) for band in table['Band']] == ['g', 'Gbp', 'Grp']
    assert table['WavelengthEff'].tolist() == pytest.approx([GAIA_EFFECTIVE_WAVELENGTHS[band] for band in ('G', 'Gbp', 'Grp')], rel=1e-9)


@pytest.mark.parametrize("name", ['LSST', 'LSST_DP0'])
def test_there_is_a_filter_for_each_magnitude_of_lsst(monkeypatch, name):
    # the saved table of LSST has two filters named g and two named i (LSST.g and LSSTv19.g, LSST.i and LSSTv19.i): the instrument is asked for
    asked = replay_svo(monkeypatch, LSST='svo_raw_LSST.ecsv')
    table = interface(name, ['u', 'g', 'r', 'i', 'z', 'y']).svofpsdrivers[0].svo_filter_list
    assert asked == ['LSST']
    assert len(table) == 6
    assert set(table['Instrument']) == {'LSST'}


# -- a facility that the service has no filter of -----------------------------------

# astroquery's words (with their typo) when the service has no filter whose facility is the one that is asked: it is what
# it says of 'Generic' and 'DENIS', whose filters have no facility in the service
# (https://github.com/astropy/astroquery/issues/2899)
NO_FILTER = 'No filter found for requested Facilty (and Instrument)'


def test_a_facility_that_the_service_has_no_filter_of_gives_no_table_and_a_warning(monkeypatch):
    asked = replay_svo(monkeypatch, Generic=IndexError(NO_FILTER))
    driver = interface('Generic__Johnson_UBVRIJHK', ['U', 'B']).svofpsdrivers[0]
    with pytest.warns(UserWarning, match=r"no filter for the facility 'Generic' and the instrument 'Johnson_UBVRIJHKL', which "
                                         r"Generic__Johnson_UBVRIJHK needs.*astroquery/issues/2899"):
        assert driver.svo_filter_list is None
    assert asked == ['Generic']


def test_the_service_that_has_no_filter_of_a_facility_is_asked_and_the_warning_given_once(monkeypatch):
    asked = replay_svo(monkeypatch, Generic=IndexError(NO_FILTER))
    driver = interface('Generic__Johnson_UBVRIJHK', ['U']).svofpsdrivers[0]
    with pytest.warns(UserWarning):
        driver.svo_filter_list
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert driver.svo_filter_list is None
    assert asked == ['Generic']


def test_an_empty_table_is_the_service_having_no_filter_too(monkeypatch):
    replay_svo(monkeypatch, Generic=Table({'filterID': np.array([], dtype=str)}))              # what facility='' gives
    driver = interface('Generic__Johnson_UBVRIJHK', ['U']).svofpsdrivers[0]
    with pytest.warns(UserWarning, match="no filter for the facility 'Generic'"):
        assert driver.svo_filter_list is None


def test_an_instrument_that_the_facility_does_not_have_is_the_service_having_no_filter_of_it(monkeypatch):
    replay_svo(monkeypatch, GAIA='svo_raw_GAIA.ecsv')
    driver = interface('GAIA__NOSUCH', ['G']).svofpsdrivers[0]
    with pytest.warns(UserWarning, match="no filter for the facility 'GAIA' and the instrument 'NOSUCH'"):
        assert driver.svo_filter_list is None


def test_another_index_error_is_not_taken_for_the_service_having_no_filter(monkeypatch):
    replay_svo(monkeypatch, Generic=IndexError('list index out of range'))
    driver = interface('Generic__Johnson_UBVRIJHK', ['U']).svofpsdrivers[0]
    with pytest.raises(IndexError, match='out of range'):
        driver.svo_filter_list


def test_a_service_that_cannot_be_reached_is_not_taken_for_having_no_filter(monkeypatch):
    replay_svo(monkeypatch, Generic=ConnectionError('no network'))
    driver = interface('Generic__Johnson_UBVRIJHK', ['U']).svofpsdrivers[0]
    with pytest.raises(ConnectionError, match='no network'):
        driver.svo_filter_list


if __name__ == '__main__':
    pass
