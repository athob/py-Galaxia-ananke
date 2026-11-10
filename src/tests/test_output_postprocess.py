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
Tests of the post-processing of the catalogue that Output does once it has read the files of the backend: the columns that the
pipelines make (checked against the pipelines applied one after the other to the table of the stars, which are tested by
themselves), what is flushed into the files, the group of the parent particles with the statistics of their stars, and what is said
of it. The stars are those of synthetic.py: six stars of three parents in two partitions, whose statistics are worked out by hand.
"""
import json
import warnings

import h5py
import numpy as np
import pandas as pd
import pytest

from ..galaxia_ananke import Output
from ..galaxia_ananke._constants import PARENTPARTICLE_GROUP, STARCATALOG_GROUP
from .synthetic import kernels, particles, prepared_survey, read_output, write_backend_output

PARENTS = [0, 0, 1, 1, 2, 2]
PARTITIONS = [0, 0, 0, 1, 1, 1]
NEW_OR_CHANGED = ['glon', 'glat', 'rad', 'mul', 'mub', 'vr', 'ra', 'dec', 'mura', 'mudec', 'pi', 'teff', 'lum']
POSTPROCESS = {'ra', 'dec', 'glon', 'glat', 'pi', 'mura', 'mudec', 'mul', 'mub', 'vr'}


def post_processed(tmp_path, monkeypatch, pp_auto_flush=True, max_pp_workers=1):
    output, written = read_output(tmp_path, monkeypatch, max_pp_workers=max_pp_workers, pp_auto_flush=pp_auto_flush,
                                  n_stars=6, partitions=PARTITIONS, parents=PARENTS)
    output.post_process_output()
    return output, written[0][1]


def group(path, name):
    with h5py.File(path, 'r') as f5:
        return {key: f5[name][key][...] for key in f5[name]}


def pipelines_applied_to(columns):
    """The stars as a table, after the pipelines that post-processing runs, applied to all of them at once."""
    df = pd.DataFrame(columns)
    Output._Output__pp_convert_cartesian_to_galactic(df)
    Output._Output__pp_convert_galactic_to_icrs(df)
    Output._Output__pp_last_conversions(df)
    return df


# -- the columns of the catalogue ---------------------------------------------------

def test_the_post_processed_files_have_every_key_of_a_catalogue(tmp_path, monkeypatch):
    output, _ = post_processed(tmp_path, monkeypatch)
    for path in output._hdf5s.values():
        assert set(group(path, STARCATALOG_GROUP)) == set(output.catalogue_keys)
    assert POSTPROCESS < set(output.catalogue_keys)


@pytest.mark.parametrize("pp_auto_flush", [True, False])
@pytest.mark.parametrize("max_pp_workers", [1, 2])
def test_the_columns_that_post_processing_makes_or_changes_are_those_of_the_pipelines(max_pp_workers, pp_auto_flush, tmp_path, monkeypatch):
    output, columns = post_processed(tmp_path, monkeypatch, pp_auto_flush=pp_auto_flush, max_pp_workers=max_pp_workers)
    expected = pipelines_applied_to(columns)
    for partition, path in output._hdf5s.items():
        stars = group(path, STARCATALOG_GROUP)
        rows = (np.asarray(columns['partitionid']) == partition)
        for key in NEW_OR_CHANGED:
            assert stars[key] == pytest.approx(expected[key].to_numpy()[rows], rel=1e-5, abs=1e-8), key


def test_the_columns_that_post_processing_does_not_change_are_left_as_they_were(tmp_path, monkeypatch):
    output, columns = post_processed(tmp_path, monkeypatch)
    for partition, path in output._hdf5s.items():
        stars = group(path, STARCATALOG_GROUP)
        rows = (np.asarray(columns['partitionid']) == partition)
        for key in ('px', 'py', 'pz', 'vx', 'vy', 'vz', 'minit', 'age', 'feh', 'parentid', 'partitionid'):
            assert stars[key].tolist() == columns[key][rows].tolist(), key


def test_the_distance_that_the_backend_gave_is_replaced_by_the_one_that_is_worked_out(tmp_path, monkeypatch):
    output, columns = post_processed(tmp_path, monkeypatch)
    assert set(columns['rad'].tolist()) == {99.0}
    x, y, z = columns['px'], columns['py'], columns['pz']
    assert np.concatenate([group(path, STARCATALOG_GROUP)['rad'] for path in output._hdf5s.values()]) == pytest.approx(np.sqrt(x**2 + y**2 + z**2), rel=1e-5)


def test_the_temperature_and_the_luminosity_of_the_files_are_no_longer_logarithms(tmp_path, monkeypatch):
    output, columns = post_processed(tmp_path, monkeypatch)
    teff = np.concatenate([group(path, STARCATALOG_GROUP)['teff'] for path in output._hdf5s.values()])
    lum = np.concatenate([group(path, STARCATALOG_GROUP)['lum'] for path in output._hdf5s.values()])
    assert teff == pytest.approx(10 ** columns['teff'].astype(float), rel=1e-5) and lum == pytest.approx(10 ** columns['lum'].astype(float), rel=1e-5)


def test_the_catalogue_that_is_loaded_afterwards_has_the_new_columns(tmp_path, monkeypatch):
    output, columns = post_processed(tmp_path, monkeypatch)
    assert POSTPROCESS <= set(output._vaex.get_column_names())
    assert len(output._vaex) == 6 and sorted(output._vaex_per_partition) == [0, 1]
    assert output._vaex['pi'].to_numpy() == pytest.approx(1.0 / output._vaex['rad'].to_numpy())


def test_a_catalogue_with_stars_is_post_processed_without_a_word_of_it_being_empty(tmp_path, monkeypatch):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        post_processed(tmp_path, monkeypatch)
    assert not [warning for warning in caught if 'no star' in str(warning.message)]


def test_what_was_post_processed_is_in_the_state_of_the_output(tmp_path, monkeypatch):
    output, _ = post_processed(tmp_path, monkeypatch)
    assert json.loads((output.output_dir / f"{output.output_name}.state").read_text())['description'] == 'export_parent_data'


def test_each_pipeline_says_that_it_runs_when_the_survey_is_verbose(tmp_path, monkeypatch, capsys):
    output, _ = read_output(tmp_path, monkeypatch, n_stars=6, partitions=PARTITIONS, parents=PARENTS)
    output.survey.verbose = True
    capsys.readouterr()
    output.post_process_output()
    said = [line for line in capsys.readouterr().out.splitlines() if line.startswith("Running ")]
    assert said == ["Running convert_cartesian_to_galactic post-processing pipeline",
                    "Running convert_galactic_to_icrs post-processing pipeline",
                    "Running last_conversions post-processing pipeline"]


# -- the parent particles -----------------------------------------------------------

def expected_parents(rows):
    """What the group of the parent particles holds for the partition whose parents are the particles of these rows."""
    p, k = particles(), kernels()
    return {'pos3': p['pos3'][rows], 'vel3': p['vel3'][rows], 'massinit': p['massinit'][rows], 'age': p['age'][rows],
            'feh': p['feh'][rows], 'h_cubic': k[rows]}


STATISTICS = {   # the stars of the partition, per parent and in all, as worked out by hand
    0: {'parentid': [0, 1], 'star_count': [2, 1], 'minit_sum': [3.0, 3.0], 'parent_partition_count': [1, 2],
        'parent_star_count': [2, 2], 'parent_minit_sum': [3.0, 7.0]},
    1: {'parentid': [1, 2], 'star_count': [1, 2], 'minit_sum': [4.0, 11.0], 'parent_partition_count': [2, 1],
        'parent_star_count': [2, 2], 'parent_minit_sum': [7.0, 11.0]},
}


@pytest.mark.parametrize("pp_auto_flush", [True, False])
def test_the_group_of_the_parent_particles_has_what_the_input_says_of_them_and_what_their_stars_are(pp_auto_flush, tmp_path, monkeypatch):
    output, _ = post_processed(tmp_path, monkeypatch, pp_auto_flush=pp_auto_flush)
    for partition, rows in {0: [0, 1], 1: [1, 2]}.items():
        parents = group(output._hdf5s[partition], PARENTPARTICLE_GROUP)
        assert set(parents) == {'parentid', 'pos3', 'vel3', 'massinit', 'age', 'feh', 'h_cubic', 'star_count', 'minit_sum',
                                'parent_partition_count', 'parent_star_count', 'parent_minit_sum'}
        for key, values in {**expected_parents(rows), **{k: v for k, v in STATISTICS[partition].items()}}.items():
            np.testing.assert_array_equal(parents[key], np.asarray(values), err_msg=f"partition {partition}, {key}")


def test_the_parents_of_a_partition_are_those_that_have_a_star_in_it_in_the_order_of_the_input(tmp_path, monkeypatch):
    output, _ = post_processed(tmp_path, monkeypatch)
    assert group(output._hdf5s[0], PARENTPARTICLE_GROUP)['parentid'].tolist() == [0, 1]
    assert group(output._hdf5s[1], PARENTPARTICLE_GROUP)['parentid'].tolist() == [1, 2]


def test_the_kernels_of_the_parents_are_scaled_by_the_factor_of_the_input(tmp_path, monkeypatch):
    survey = prepared_survey(tmp_path, monkeypatch, k_factor=2.0)
    write_backend_output(survey, n_stars=6, partitions=PARTITIONS, parents=PARENTS)
    survey.output.read_galaxia_output(None, 1, True)
    survey.output.post_process_output()
    assert group(survey.output._hdf5s[0], PARENTPARTICLE_GROUP)['h_cubic'].tolist() == (2.0 * kernels()[[0, 1]]).tolist()


def test_the_statistics_of_the_stars_per_parent_are_those_of_the_catalogue(tmp_path, monkeypatch):
    output, columns = post_processed(tmp_path, monkeypatch)
    stats = output._parent_stats_per_partitionid
    assert sorted(stats) == [0, 1] and stats[0].index.name == 'parentid'
    assert stats[0].loc[1, 'star_count'] == 1 and stats[1].loc[1, 'parent_star_count'] == 2 and stats[1].loc[2, 'minit_sum'] == 11.0


if __name__ == '__main__':
    pass
