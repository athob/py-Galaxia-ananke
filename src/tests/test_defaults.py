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
Tests of the defaults of the package: the names, the box of the color-magnitude diagram that cuts nothing, the units, the
position and the velocity of the Sun that astropy gives (which a survey takes as its viewpoint), and the parameters that a
survey starts from.
"""
import dataclasses
import pathlib

import pytest
from astropy import units

from ..galaxia_ananke._defaults import (DEF_UNIT, DEFAULT_CMD, DEFAULT_CMD_BOX, DEFAULT_PARFILE, DEFAULT_PSYS,
                                        DEFAULT_SIMNAME, DEFAULT_SURVEYNAME, DEFAULTS_FOR_PARFILE, GALAXIA_TMP,
                                        DefaultUnits, rSun, vSun)
from ..galaxia_ananke._templates import FTTAGS


def test_the_default_names():
    assert (DEFAULT_SIMNAME, DEFAULT_SURVEYNAME, DEFAULT_PARFILE) == ("sim", "survey", "survey_params")


def test_the_default_photometric_system_and_diagram():
    assert DEFAULT_PSYS == ["padova/GAIA__DR2"]
    assert DEFAULT_CMD == "G,Gbp-Grp"


def test_the_default_box_of_the_diagram_cuts_nothing_but_the_faintest_absolute_magnitudes():
    assert DEFAULT_CMD_BOX == {"app_mag": [-1000, 1000], "abs_mag": [-1000, 20], "color": [-1000, 1000]}


def test_the_temporary_directory_is_a_directory_that_exists():
    assert GALAXIA_TMP.is_absolute() and GALAXIA_TMP.is_dir()


# -- the units ----------------------------------------------------------------------

def test_the_default_units():
    assert DEF_UNIT.position == units.kpc
    assert DEF_UNIT.velocity == units.km / units.s
    assert DEF_UNIT.wavelength == units.micron
    assert DEF_UNIT.spectral == units.Jy
    assert DEF_UNIT.irradiance == units.Unit("erg/(cm2 s)")


def test_the_default_units_are_one_unchangeable_set():
    assert DefaultUnits() is DEF_UNIT
    with pytest.raises(dataclasses.FrozenInstanceError):
        DEF_UNIT.position = units.pc


# -- the Sun ------------------------------------------------------------------------

def test_the_viewpoint_of_a_survey_is_the_sun_of_astropy_as_it_is_today():
    # What astropy takes as the position and the velocity of the Sun in the Galaxy is the viewpoint of every survey:
    # if a version of astropy changes it, the surveys change with it, and this is where it would be seen.
    assert list(rSun) == pytest.approx([-8.121973366106948, -1.4729444519452402e-11, 0.020800000030482345], rel=1e-9, abs=1e-12)
    assert list(vSun) == pytest.approx([12.9, 245.6, 7.78], rel=1e-9)


# -- what a survey starts from ------------------------------------------------------

def test_the_defaults_and_what_the_survey_adds_are_all_the_tags_of_the_parameter_file():
    # what a survey writes by itself is its names, its output file and its photometric system
    added_by_the_survey = {FTTAGS.name, FTTAGS.pname, FTTAGS.output_file, FTTAGS.photo_categ, FTTAGS.photo_sys, FTTAGS.mag_color_names}
    every_tag = {field.name for field in dataclasses.fields(FTTAGS)}
    assert set(DEFAULTS_FOR_PARFILE) == every_tag - added_by_the_survey


def test_the_default_box_is_the_one_of_the_parameters():
    parameters = DEFAULTS_FOR_PARFILE
    assert [parameters[FTTAGS.app_mag_lim_lo], parameters[FTTAGS.app_mag_lim_hi]] == DEFAULT_CMD_BOX["app_mag"]
    assert [parameters[FTTAGS.abs_mag_lim_lo], parameters[FTTAGS.abs_mag_lim_hi]] == DEFAULT_CMD_BOX["abs_mag"]
    assert [parameters[FTTAGS.color_lim_lo], parameters[FTTAGS.color_lim_hi]] == DEFAULT_CMD_BOX["color"]


def test_the_viewpoint_of_the_parameters_is_the_sun():
    for tag, value in zip(FTTAGS.rSun, rSun):
        assert DEFAULTS_FOR_PARFILE[tag] == value
    for tag, value in zip(FTTAGS.vSun, vSun):
        assert DEFAULTS_FOR_PARFILE[tag] == value


def test_the_output_directory_is_the_temporary_one():
    assert DEFAULTS_FOR_PARFILE[FTTAGS.output_dir] == GALAXIA_TMP


@pytest.mark.parametrize("tag, value", [
    (FTTAGS.fsample, 1),            # a survey takes every star unless told otherwise
    (FTTAGS.rand_seed, 17052),
    (FTTAGS.nres, 64),              # the number of neighbours that Input's own default should be (it is the name of the tag, see Input)
    (FTTAGS.nstart, 0),
    (FTTAGS.r_max, 500),
    (FTTAGS.r_min, 0),
    (FTTAGS.geometry_opt, 0),
    (FTTAGS.pop_id, 10),
    (FTTAGS.star_type, 0),
    (FTTAGS.photo_error, 0),
    (FTTAGS.warp_flare_on, 0),
])
def test_the_defaults_a_survey_starts_from(tag, value):
    assert DEFAULTS_FOR_PARFILE[tag] == value


if __name__ == '__main__':
    pass
