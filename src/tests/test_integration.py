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
Integration tests: what the compiled backend gives, through the whole of the package, for a small survey. They check what has to
hold whatever the numbers are (the layout of the files, the keys and the types of the catalogue, the limits that were asked for,
relations between the columns, the consistency of what is said of the parent particles), and that the files of the backend are laid
out as the synthetic ones that the other tests of Output are made with.
"""
import json

import ebf
import h5py
import numpy as np
import pytest

from ..galaxia_ananke._constants import PARENTPARTICLE_GROUP, PARTITIONING_GROUP, STARCATALOG_GROUP
from ..galaxia_ananke.Output import Output
from .real_backend import LIMITS, real_survey                 # noqa: F401 (the fixture)
from .synthetic import INTEGER_KEYS, dtype_of
from .utils import needs_backend

pytestmark = [pytest.mark.integration, needs_backend]

FLOAT32 = ['gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp', 'grav', 'lum', 'mact', 'mtip', 'teff']
PARENT_DATASETS = {'parentid', 'pos3', 'vel3', 'massinit', 'age', 'feh', 'h_cubic', 'star_count', 'minit_sum',
                   'parent_partition_count', 'parent_star_count', 'parent_minit_sum'}


def column(real_survey, key):
    return real_survey.output._vaex[key].to_numpy()


def read_h5(path):
    with h5py.File(path, 'r') as f5:
        return {'groups': {group: {key: f5[group][key][...] for key in f5[group]} for group in f5}, 'attrs': dict(f5.attrs)}


# -- what was made ------------------------------------------------------------------

def test_the_survey_has_stars_in_one_partition_and_one_file(real_survey):
    output = real_survey.output
    assert len(output._vaex) > 100
    assert list(output._hdf5s) == [0] and output._hdf5s[0].is_file()
    assert len(output._ebfs) == 1


def test_the_catalogue_has_the_keys_of_the_package_and_no_other(real_survey):
    output = real_survey.output
    assert len(output.export_keys) == 34 and len(output.catalogue_keys) == 44          # with the optional keys of the dummy input
    assert set(output._vaex.get_column_names()) == set(output.catalogue_keys)


def test_the_catalogue_has_the_types_of_its_columns(real_survey):
    for key in FLOAT32:
        assert column(real_survey, key).dtype == np.float32, key
    for key in INTEGER_KEYS:
        assert column(real_survey, key).dtype == np.int64, key


def test_the_file_has_the_groups_the_attributes_and_the_datasets_that_it_is_documented_to_have(real_survey):
    h5 = read_h5(real_survey.output._hdf5s[0])
    assert set(h5['groups']) == {STARCATALOG_GROUP, PARENTPARTICLE_GROUP, PARTITIONING_GROUP}
    assert set(h5['groups'][STARCATALOG_GROUP]) == set(real_survey.output.catalogue_keys)
    assert set(h5['groups'][PARENTPARTICLE_GROUP]) == PARENT_DATASETS
    assert h5['groups'][PARTITIONING_GROUP]['partition_id'].tolist() == [0] and h5['attrs']['Cat_partition_id'] == 0
    for prefix in ('Input_', 'Survey_', 'SurveyParam_'):
        assert any(key.startswith(prefix) for key in h5['attrs']), prefix
    assert {'Header', 'All_partition_ids', 'Galaxia_ebf_output_logs'} <= set(h5['attrs'])


def test_the_attributes_of_the_file_are_the_metadata_of_the_output(real_survey):
    attrs = read_h5(real_survey.output._hdf5s[0])['attrs']
    assert {key: attrs[key] for key in real_survey.output.all_metadata} == real_survey.output.all_metadata


def test_the_state_of_the_output_says_that_the_parent_data_was_exported(real_survey):
    output = real_survey.output
    assert json.loads((output.output_dir / f"{output.output_name}.state").read_text())['description'] == 'export_parent_data'


# -- what has to hold whatever the numbers are --------------------------------------

def test_the_stars_respect_the_limits_that_were_asked_for(real_survey):
    g, dmod = column(real_survey, 'gaia__dr2_g'), column(real_survey, 'dmod')
    assert np.isfinite(g).all() and g.max() <= LIMITS['app_mag_lim_hi'] + 1e-3          # apparent magnitude
    assert (g - dmod).max() <= LIMITS['abs_mag_lim_hi'] + 1e-3                          # absolute magnitude


def test_the_distance_modulus_is_that_of_the_distance(real_survey):
    rad, dmod = column(real_survey, 'rad'), column(real_survey, 'dmod')
    assert (rad > 0).all() and np.isfinite(rad).all()
    assert dmod == pytest.approx(5 * np.log10(rad) + 10, abs=1e-2)                      # the distance is in kpc


def test_the_parallax_is_the_inverse_of_the_distance(real_survey):
    assert column(real_survey, 'pi') == pytest.approx(1.0 / column(real_survey, 'rad'), rel=1e-5)


def test_the_angles_are_within_their_ranges(real_survey):
    glon, glat, ra, dec = (column(real_survey, key) for key in ('glon', 'glat', 'ra', 'dec'))
    assert (glon > -180).all() and (glon <= 180).all() and (np.abs(glat) <= 90).all()
    assert (ra >= 0).all() and (ra < 360).all() and (np.abs(dec) <= 90).all()


def test_the_temperatures_and_the_luminosities_are_no_longer_logarithms(real_survey):
    assert (column(real_survey, 'teff') > 1000).all() and (column(real_survey, 'lum') > 0).all()


def test_the_stars_are_not_more_massive_now_than_they_were_born_but_for_the_precision_of_the_backend(real_survey):
    # The current mass is the initial mass less what was lost, but not to the precision of the backend: it is found by interpolating
    # the isochrones, which gives up to 0.06% more than the initial mass for about one star in ten (1888 stars, 227 of them above).
    mact, minit = column(real_survey, 'mact'), column(real_survey, 'minit').astype(float)
    ratio = mact / minit
    assert (minit > 0).all() and (ratio <= 1.01).all() and (ratio > 0.5).all()
    assert 0.99 < np.median(ratio) <= 1.001


def test_the_stars_are_the_children_of_the_particles_that_were_given(real_survey):
    parentid = column(real_survey, 'parentid')
    assert parentid.min() >= 0 and parentid.max() < len(real_survey.particles['massinit'])
    assert (column(real_survey, 'partitionid') == 0).all()


# -- the parent particles -----------------------------------------------------------

def test_what_is_said_of_the_parent_particles_is_consistent_with_the_catalogue(real_survey):
    parents = read_h5(real_survey.output._hdf5s[0])['groups'][PARENTPARTICLE_GROUP]
    n_stars, minit = len(real_survey.output._vaex), column(real_survey, 'minit')
    assert parents['star_count'].sum() == n_stars and (parents['star_count'] >= 1).all()
    assert parents['minit_sum'].sum() == pytest.approx(minit.astype(float).sum(), rel=1e-5)
    assert len(set(parents['parentid'].tolist())) == len(parents['parentid'])           # each parent once
    assert set(parents['parentid'].tolist()) == set(column(real_survey, 'parentid').tolist())
    assert (parents['parent_partition_count'] == 1).all()                               # one partition
    assert (parents['parent_star_count'] == parents['star_count']).all()


def test_what_is_said_of_the_parent_particles_is_what_the_input_says_of_them(real_survey):
    parents = read_h5(real_survey.output._hdf5s[0])['groups'][PARENTPARTICLE_GROUP]
    ids, given = parents['parentid'], real_survey.particles
    for key in ('pos3', 'vel3', 'massinit', 'age', 'feh'):
        np.testing.assert_array_equal(parents[key], given[key][ids], err_msg=key)
    np.testing.assert_array_equal(parents['h_cubic'], real_survey.kernels[ids])


# -- the files of the backend are laid out as the synthetic ones --------------------

def test_the_file_of_the_backend_has_a_dataset_for_each_key_that_the_synthetic_files_have_with_the_same_type(real_survey):
    path = str(real_survey.output._ebfs[0])
    datasets = {key.decode().lstrip('/') for key in ebf._EbfMap.keys(path)}
    assert set(real_survey.output.export_keys) <= datasets and 'log' in datasets
    magnitudes = {key for photosystem in real_survey.output.photosystems for key in photosystem.to_export_keys}
    for key in real_survey.output.export_keys:
        assert ebf.getHeader(path, f"/{key}").get_dtype() == dtype_of(key, magnitudes), key


def test_the_log_of_the_backend_is_made_of_a_header_the_parameters_and_a_footer(real_survey):
    parsed = Output._Output__parse_galaxia_ebf_log(str(real_survey.output._ebfs[0].resolve()))
    assert isinstance(parsed['log_header'], str) and isinstance(parsed['log_footer'], str) and len(parsed) > 10
    assert float(parsed['fSample']) == pytest.approx(LIMITS['fsample'])


def test_the_first_parameter_of_the_log_of_the_backend_is_read_by_its_name(real_survey):
    # the backend writes it without the comment sign that the parser takes for its name when it is there (see test_output_export.py)
    path = real_survey.output._ebfs[0]
    parsed = Output._Output__parse_galaxia_ebf_log(str(path.resolve()))
    assert '#' not in parsed and parsed['outputFile'] == path.stem


if __name__ == '__main__':
    pass
