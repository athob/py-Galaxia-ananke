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
Tests of Output as the survey that made it uses it, before and around the reading of the output of the backend: what it says
about itself (names, places, keys, parameters, metadata), the files of the backend that it looks for and clears, and what it
refuses before it has anything to read. The backend's files are synthetic (see synthetic.py).
"""
import numpy as np
import pytest

from ..galaxia_ananke._defaults import DEFAULTS_FOR_PARFILE, DEFAULT_CMD, rSun
from ..galaxia_ananke._templates import FTTAGS
from .synthetic import particles, prepared_survey, write_backend_output, write_ebf

EXPORT = {'px', 'py', 'pz', 'vx', 'vy', 'vz', 'rad', 'dmod', 'mtip', 'mact', 'minit', 'age', 'grav', 'feh', 'teff', 'lum',
          'parentid', 'partid', 'partitionid'}
POSTPROCESS = {'ra', 'dec', 'glon', 'glat', 'pi', 'mura', 'mudec', 'mul', 'mub', 'vr'}
GAIA_MAGNITUDES = {'gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp'}


@pytest.fixture
def survey(tmp_path, monkeypatch):
    return prepared_survey(tmp_path, monkeypatch)


# -- what an output says about itself -----------------------------------------------

def test_an_output_belongs_to_the_survey_that_made_it(survey):
    assert survey.output.survey is survey and survey.output.photosystems is survey.photosystems
    assert survey.output.caching is False and survey.output.verbose is False


def test_the_places_and_the_name_of_an_output(survey, tmp_path):
    output = survey.output
    name = f"{survey.surveyname_hash}.{survey.inputname_hash}"
    assert output.output_dir == tmp_path / 'out' and output.output_name == name
    assert output._file_base == tmp_path / 'out' / name
    assert output._ebf_glob_pattern == tmp_path / 'out' / f"{name}.*.ebf"
    assert output._hdf5_glob_pattern == tmp_path / 'out' / f"{name}.*.h5"


def test_what_an_output_says_of_the_parameters_of_its_survey(survey):
    output = survey.output
    assert output.parameters == survey.parameters
    assert output.parameter_mag_color_names == DEFAULT_CMD == 'G,Gbp-Grp' and output.parameter_magnitude_name == 'G'
    assert output.parameter_abs_mag_hi == DEFAULTS_FOR_PARFILE[FTTAGS.abs_mag_lim_hi]
    assert output.parameter_app_mag_hi == DEFAULTS_FOR_PARFILE[FTTAGS.app_mag_lim_hi]


def test_the_position_of_the_sun_is_a_sky_coordinate_in_the_galactic_frame(survey):
    sun = survey.output.rsun_skycoord
    assert sun.frame.name == 'galactic'
    assert sun.cartesian.xyz.to_value('kpc').tolist() == pytest.approx(list(rSun))


def test_the_keys_that_an_output_exports_and_the_keys_of_its_catalogue(survey):
    assert set(survey.output.export_keys) == EXPORT | GAIA_MAGNITUDES
    assert set(survey.output.catalogue_keys) == EXPORT | GAIA_MAGNITUDES | POSTPROCESS
    assert len(survey.output.export_keys) == len(EXPORT | GAIA_MAGNITUDES)


def test_the_keys_of_the_optional_properties_of_the_input_are_exported_and_the_id_is_the_id_of_the_satellite(tmp_path, monkeypatch):
    survey = prepared_survey(tmp_path, monkeypatch, input_particles=particles(id=np.array([1, 2, 3]), dform=np.zeros(3)))
    assert set(survey.output._make_input_optional_keys()) == {'satid', 'dform', 'parentid', 'partitionid'}
    assert {'satid', 'dform'} <= set(survey.output.export_keys) and 'id' not in survey.output.export_keys


def test_the_keys_of_the_optional_properties_of_the_input_are_in_the_catalogue_too(tmp_path, monkeypatch):
    survey = prepared_survey(tmp_path, monkeypatch, input_particles=particles(id=np.array([1, 2, 3]), dform=np.zeros(3)))
    assert {'satid', 'dform'} | POSTPROCESS <= set(survey.output.catalogue_keys) and 'id' not in survey.output.catalogue_keys
    assert set(survey.output.catalogue_keys) == set(survey.output.export_keys) | POSTPROCESS


def test_the_magnitudes_of_all_the_photometric_systems_are_exported(tmp_path, monkeypatch):
    survey = prepared_survey(tmp_path, monkeypatch, photo_sys=['padova/GAIA__DR2', 'padova/LSST'])
    assert set(survey.output.export_keys) == EXPORT | GAIA_MAGNITUDES | {'lsst_u', 'lsst_g'}


def test_all_the_metadata_of_an_output_is_that_of_its_parameters_its_input_and_its_survey(survey, tmp_path):
    metadata = survey.output.all_metadata
    assert len(metadata) == len(survey.parameters) + len(survey.input.metadata) + len(survey.metadata)
    assert metadata['SurveyParam_fsample'] == 1 and metadata[f"SurveyParam_{FTTAGS.output_dir}"] == str(tmp_path / 'out')
    assert (metadata['Input_name'], metadata['Input_hash']) == ('sim', survey.input.hash)
    assert (metadata['Survey_surveyname'], metadata['Survey_hash']) == ('survey', survey.hash)
    assert all(isinstance(value, (bool, int, float, str)) for value in metadata.values())


# -- the files of the backend --------------------------------------------------------

def test_an_output_without_files_of_the_backend_says_that_there_are_none(survey):
    assert list(survey.output._ebf_glob) == []
    with pytest.raises(RuntimeError, match="Galaxia didn't return any ebf output files"):
        survey.output._ebfs


def test_the_files_of_the_backend_are_found_when_they_are_there(survey):
    with pytest.raises(RuntimeError):
        survey.output._ebfs
    written = write_backend_output(survey)
    assert set(survey.output._ebfs) == {path for path, _ in written}


def test_the_files_of_every_generation_are_found(tmp_path, monkeypatch):
    survey = prepared_survey(tmp_path, monkeypatch, n_gens=(0, 1, 2))
    written = write_backend_output(survey)
    assert len(survey.output._ebfs) == 3 and set(survey.output._ebfs) == {path for path, _ in written}


def test_the_files_of_the_backend_of_another_output_are_not_found(survey, tmp_path):
    write_ebf(tmp_path / 'out' / 'another.0.ebf', {})
    write_backend_output(survey)
    assert len(survey.output._ebfs) == 1


def test_making_the_output_again_clears_the_files_of_the_backend_that_it_would_read(survey, tmp_path):
    written = write_backend_output(survey)
    other = write_ebf(tmp_path / 'out' / 'another.0.ebf', {})
    survey._prepare_survey_parameters_and_output(DEFAULT_CMD, (0,), output_dir=tmp_path / 'out', fsample=1)
    assert not written[0][0].exists() and other.exists()


def test_with_caching_the_files_of_the_backend_are_kept(tmp_path, monkeypatch):
    with pytest.warns(DeprecationWarning, match="caching mode"):
        survey = prepared_survey(tmp_path, monkeypatch, caching=True)
    written = write_backend_output(survey)
    survey._prepare_survey_parameters_and_output(DEFAULT_CMD, (0,), output_dir=tmp_path / 'out', fsample=1)
    assert written[0][0].exists()


@pytest.mark.parametrize("answer, removed", [('y', True), ('n', False), ('', False), ('yes', False)])
def test_without_being_forced_the_files_are_only_removed_for_a_y(answer, removed, survey, monkeypatch):
    path = write_backend_output(survey)[0][0]
    prompts = []
    monkeypatch.setattr('builtins.input', lambda prompt: prompts.append(prompt) or answer)
    survey.output._Output__clear_ebfs(force=False)
    assert prompts == [f"You are about to remove {path}, do I proceed? [y/N] "] and path.exists() is not removed


def test_a_file_of_the_backend_that_went_away_is_refused(survey):
    written = write_backend_output(survey, partitions=[0, 0, 1, 1, 1])
    survey.output._ebfs
    written[0][0].unlink()
    with pytest.raises(RuntimeError, match="Don't attempt creating an Output object on your own"):
        survey.output._Output__ebfs_partitions


# -- before there is anything to read ------------------------------------------------

@pytest.mark.parametrize("what", ["_vaex", "_vaex_per_partition"])
def test_an_output_that_has_read_nothing_refuses_to_give_its_catalogue(what, survey):
    with pytest.raises(RuntimeError, match="Don't attempt creating an Output object on your own"):
        getattr(survey.output, what)


def test_an_output_that_has_read_nothing_refuses_to_be_shown_or_indexed(survey):
    with pytest.raises(RuntimeError):
        repr(survey.output)
    with pytest.raises(RuntimeError):
        survey.output['px']
    with pytest.raises(RuntimeError):
        survey.output['px'] = 1


def test_an_attribute_that_does_not_exist_is_an_attribute_error_even_before_there_is_a_catalogue(survey):
    with pytest.raises(AttributeError):
        survey.output.nothing_of_the_sort
    assert 'export_keys' in dir(survey.output)


def test_what_is_not_implemented(survey):
    with pytest.raises(NotImplementedError):
        survey.output.save('somewhere')
    with pytest.raises(NotImplementedError):
        survey.output._path


def test_the_defaults_of_the_reading(survey):
    assert survey.output._max_pp_workers == 1 and survey.output._pp_auto_flush is True


if __name__ == '__main__':
    pass
