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
Tests of the templates of the package, which write the files and the commands that the backend is given: the tags that
fill them, the parameter file (in the format that the backend reads, as the surveys of the maintainer wrote it), the
file that lists the input files, and the commands that run and that append to a survey.
"""
import dataclasses
import re

import pytest

from ..galaxia_ananke._constants import GALAXIA
from ..galaxia_ananke._templates import (APPEND_TEMPLATE, CTTAGS, FILENAME_TEMPLATE, FTTAGS, HDIMBLOCK_TEMPLATE,
                                         PARFILE_TEMPLATE, PARFILENAME_TEMPLATE, RUN_TEMPLATE,
                                         CommandTemplateTags, FileTemplateTags)


def placeholders(template):
    return re.findall(r"\$\{(\w+)\}", template.template)


def tags_of(dataclass_instance):
    return [field.name for field in dataclasses.fields(dataclass_instance)]


# -- the tags -----------------------------------------------------------------------

@pytest.mark.parametrize("tags", [FTTAGS, CTTAGS], ids=["file", "command"])
def test_every_tag_is_its_own_name(tags):
    # the templates are filled by keyword arguments: a tag that is not its name would never find its place
    for name in tags_of(tags):
        assert getattr(tags, name) == name


def test_the_tags_are_singletons():
    assert FileTemplateTags() is FTTAGS
    assert CommandTemplateTags() is CTTAGS


def test_the_tags_cannot_be_changed():
    with pytest.raises(dataclasses.FrozenInstanceError):
        FTTAGS.name = "another"
    with pytest.raises(dataclasses.FrozenInstanceError):
        CTTAGS.nfile = "another"


def test_the_position_and_velocity_of_the_sun_have_their_three_tags_in_order():
    assert FTTAGS.rSun == ["rSun0", "rSun1", "rSun2"]
    assert FTTAGS.vSun == ["vSun0", "vSun1", "vSun2"]


def test_the_tags_of_the_photometric_systems_after_the_first_are_numbered():
    assert FTTAGS.append_photo_categ(2) == "photo_categ_2"
    assert FTTAGS.append_photo_sys(2) == "photo_sys_2"
    assert FTTAGS.append_photo(3) == ["photo_categ_3", "photo_sys_3"]


# -- the parameter file -------------------------------------------------------------

SAMPLE_VALUES = {
    "output_file": "survey_cebfb6f",
    "output_dir": "/work/out",
    "photo_categ": "padova",
    "photo_sys": "GAIA__DR2",
    "mag_color_names": "G,Gbp-Grp",
    "app_mag_lim_lo": -1000, "app_mag_lim_hi": 15.0,
    "abs_mag_lim_lo": -1000, "abs_mag_lim_hi": 0.0,
    "color_lim_lo": -1000, "color_lim_hi": 1000,
    "geometry_opt": 0, "survey_area": 207.455, "fsample": 0.01, "pop_id": 10, "warp_flare_on": 0,
    "longitude": 76.273, "latitude": 13.4725, "star_type": 0, "photo_error": 0, "rand_seed": 17052,
    "r_max": 500, "r_min": 0, "nres": 64, "nstart": 0,
    "rSun0": -8.121973366106948, "rSun1": -1.4729444519452402e-11, "rSun2": 0.020800000030482345,
    "vSun0": 12.9, "vSun1": 245.6, "vSun2": 7.78,
}

# The parameter file of a survey of the maintainer, as the backend was given it (the tabs are the ones of the file)
EXPECTED_PARFILE = """outputFile\tsurvey_cebfb6f\t#don't fiddle
outputDir\t/work/out\t#where to output the survey
photoCateg\tpadova\t#name of folder where to select magnitude system
photoSys\tGAIA__DR2\t#magnitude system (see ananke-for-wings/GalaxiaData/Isochrones/padova/ for options)
magcolorNames\tG,Gbp-Grp\t#magnitude and color to use for selecting the CMD box
appMagLimits[0]\t-1000\t#upper and lower limits in apparent mag
appMagLimits[1]\t15.0
absMagLimits[0]\t-1000\t#upper and lower limits in absolute mag
absMagLimits[1]\t0.0
colorLimits[0]\t-1000\t#upper and lower limits in color defined on line 4
colorLimits[1]\t1000
geometryOption\t0\t#don't fiddle
surveyArea\t207.455\t#not used
fSample\t0.01\t#don't fiddle
popID\t10\t#don't fiddle
warpFlareOn\t0\t#not used
longitude\t76.273\t#not used
latitude\t13.4725\t#not used
starType\t0\t#don't fiddle
photoError\t0\t#not used
seed\t17052\t#change if you want a different random sample
r_max\t500\t#max distance from galactic center to include
r_min\t0\t#min distance from galactic center to include
nres\t64\t#nres
nstart\t0\t#integer at which to start numbering synthetic stars
rSun[0]\t-8.121973366106948\t#location of survey viewpoint relative to galactic center
rSun[1]\t-1.4729444519452402e-11
rSun[2]\t0.020800000030482345
vSun[0]\t12.9\t#velocity of survey viewpoint relative to galactic center (not used)
vSun[1]\t245.6
vSun[2]\t7.78
"""


def test_the_parameter_file_is_written_as_the_backend_reads_it():
    assert PARFILE_TEMPLATE.substitute(SAMPLE_VALUES) == EXPECTED_PARFILE


def test_the_parameter_file_uses_each_tag_once_but_the_names_of_the_input():
    used = placeholders(PARFILE_TEMPLATE)
    assert sorted(used) == sorted(set(tags_of(FTTAGS)) - {"name", "pname"})


def test_a_value_that_is_missing_is_an_error():
    for missing in SAMPLE_VALUES:
        values = {key: value for key, value in SAMPLE_VALUES.items() if key != missing}
        with pytest.raises(KeyError):
            PARFILE_TEMPLATE.substitute(values)


def test_the_sample_values_are_all_the_tags_of_the_parameter_file():
    assert sorted(SAMPLE_VALUES) == sorted(placeholders(PARFILE_TEMPLATE))


def test_the_parameter_file_is_named_after_the_input():
    assert PARFILENAME_TEMPLATE.substitute(name="sim_6c0ceae") == "sim_6c0ceae_params"


# -- the file that lists the input files --------------------------------------------

def test_the_list_of_input_files_names_the_directory_the_one_satellite_and_the_file():
    assert FILENAME_TEMPLATE.substitute(name="sim_6c0ceae", pname="sim_6c0ceae.ebf") == "nbody1/sim_6c0ceae/\n\t1\t1\nsim_6c0ceae.ebf\n"


# -- the commands -------------------------------------------------------------------

def test_the_run_command_is_the_backend_with_its_options_and_the_parameter_file():
    hdim_block = HDIMBLOCK_TEMPLATE.substitute(hdim=6)
    assert hdim_block == " --hdim=6"
    command = RUN_TEMPLATE.substitute(hdim_block=hdim_block, nfile="sim_6c0ceae", ngen=0, parfile="/work/in/survey_params")
    assert command == f"{GALAXIA} -r --hdim=6 --nfile=sim_6c0ceae --ngen=0 /work/in/survey_params"


def test_the_run_command_without_a_dimension_has_no_such_option():
    command = RUN_TEMPLATE.substitute(hdim_block="", nfile="sim", ngen=1, parfile="params")
    assert command == f"{GALAXIA} -r --nfile=sim --ngen=1 params"


def test_the_append_command_gives_the_system_and_the_file_to_add_to():
    command = APPEND_TEMPLATE.substitute(pcat="padova", psys="GAIA__DR2", filename="survey.0.ebf")
    assert command == f"{GALAXIA} -a --pcat=padova --psys=GAIA__DR2 survey.0.ebf"


def test_the_commands_use_every_command_tag_and_only_them():
    used = placeholders(RUN_TEMPLATE) + placeholders(HDIMBLOCK_TEMPLATE) + placeholders(APPEND_TEMPLATE)
    assert sorted(set(used)) == sorted(tags_of(CTTAGS))


def test_the_commands_are_cut_into_the_arguments_that_execute_gives_the_backend():
    # execute splits a command on its spaces: each option has to be one word
    command = RUN_TEMPLATE.substitute(hdim_block=HDIMBLOCK_TEMPLATE.substitute(hdim=3), nfile="n", ngen=0, parfile="p")
    assert command.split(" ") == [str(GALAXIA), "-r", "--hdim=3", "--nfile=n", "--ngen=0", "p"]


if __name__ == '__main__':
    pass
