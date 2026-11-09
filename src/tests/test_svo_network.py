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
Tests that use the network: the filter profile service of SVO, where the wavelengths and the zero points of the photometric systems
come from. They only run where GALAXIA_TEST_NETWORK is set (a survey does not need them: only these two properties of a system do).
They check that the table that the other tests replay from a file is still what the service answers, so that those tests do not
go on checking a service that is no longer there.
"""
import pytest
from astropy import units

from ..galaxia_ananke import photometry
from .test_photo_system import GAIA_WAVELENGTHS, GAIA_ZERO_POINTS      # what the saved table of the service gives
from .utils import needs_network

pytestmark = [pytest.mark.network, needs_network]


def gaia():
    return photometry.available_photo_systems['padova']['GAIA__DR2']


def test_the_effective_wavelengths_of_gaia_are_those_of_the_saved_table_of_the_service():
    wavelengths = gaia().effective_wavelengths
    assert wavelengths.unit == units.micron and wavelengths.value == pytest.approx(GAIA_WAVELENGTHS, rel=1e-6)


def test_the_zero_points_of_gaia_are_those_of_the_saved_table_of_the_service():
    zeropoints = gaia().zeropoints
    assert zeropoints.unit == units.Jy and zeropoints.value == pytest.approx(GAIA_ZERO_POINTS, rel=1e-6)


if __name__ == '__main__':
    pass
