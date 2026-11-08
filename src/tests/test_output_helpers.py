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
Tests of what Output is made of that does not need a survey: the longitude that is brought within bounds, the partitions that are
shared between processes, the keys of the properties that a catalogue has, and the documentation of the constructor. The
conversions of the coordinates are in test_output_pipelines.py, the HDF5 helpers in test_output_hdf5.py.
"""
import types

import numpy as np
import pytest

from ..galaxia_ananke import Output
from ..galaxia_ananke._constants import STARCATALOG_GROUP, PARENTPARTICLE_GROUP
from ..galaxia_ananke.photometry.PhotoSystem import PhotoSystem
from ..galaxia_ananke.Output import shift_g_lon     # the function, and not the module that the class takes its name from


def photo_system(name, mag_names):
    return PhotoSystem(types.SimpleNamespace(category='padova', name=name, mag_names=list(mag_names)))


GAIA = photo_system('GAIA__DR2', ['G', 'Gbp', 'Grp'])
LSST = photo_system('LSST', ['u', 'g'])


# -- the longitude ------------------------------------------------------------------

@pytest.mark.parametrize("longitude, shifted", [
    (0, 0), (45, 45), (-45, -45), (90, 90), (-90, -90), (180, 180), (-180, 180),         # -180 is 180: the interval is (-180, 180]
    (190, -170), (-190, 170), (270, -90), (359, -1), (360, 0), (720, 0),
])
def test_a_longitude_is_brought_within_minus_180_and_180(longitude, shifted):
    assert shift_g_lon(longitude) == pytest.approx(shifted)


def test_longitudes_are_shifted_together_and_always_end_up_within_the_bounds():
    longitudes = np.random.default_rng(3).uniform(-720, 720, 500)
    shifted = shift_g_lon(longitudes)
    assert np.all(shifted > -180) and np.all(shifted <= 180)
    assert np.allclose(np.cos(np.radians(shifted)), np.cos(np.radians(longitudes))) and np.allclose(np.sin(np.radians(shifted)), np.sin(np.radians(longitudes)))


# -- the partitions that a process takes --------------------------------------------

def consolidate(partitions, lengths, max_workers):
    result = Output._Output__consolidate_partitions_per_process(np.array(partitions), np.array(lengths), max_workers)
    return {int(process): [int(partition) for partition in group] for process, group in result.items()}


@pytest.mark.parametrize("partitions, lengths, max_workers, expected", [
    ([0, 1, 2, 3], [10, 40, 20, 30], 2, {0: [0, 2], 1: [1, 3]}),
    ([0, 1, 2, 3], [10, 40, 20, 30], 1, {0: [0, 2, 3, 1]}),                   # one process takes them all, from the smallest
    ([0, 1], [10, 30], 5, {0: [0], 1: [1]}),                                   # more workers than partitions: one each
    ([5], [0], 3, {0: [5]}),                                                   # nothing in it: still taken
    ([0, 1, 2, 3, 4, 5], [60, 10, 50, 20, 40, 30], 3, {0: [1, 3, 5], 1: [4], 2: [2, 0]}),   # the order alternates from a process to the next
    ([7, 8, 9], [100, 1, 10], 2, {0: [8, 9], 1: [7]}),                         # the labels of the partitions are not their position
    ([0, 1, 2], [1, 2, 3], 3, {0: [0], 1: [1], 2: [2]}),
])
def test_the_partitions_are_shared_between_the_processes_by_the_number_of_stars_in_them(partitions, lengths, max_workers, expected):
    assert consolidate(partitions, lengths, max_workers) == expected


def test_all_the_partitions_are_taken_once_whatever_the_number_of_workers():
    partitions, lengths = list(range(9)), [3, 70, 12, 45, 8, 60, 21, 90, 33]
    for workers in range(1, 12):
        taken = [partition for group in consolidate(partitions, lengths, workers).values() for partition in group]
        assert sorted(taken) == partitions


# -- the keys -----------------------------------------------------------------------

EXPORT = {'px', 'py', 'pz', 'vx', 'vy', 'vz', 'rad', 'dmod', 'mtip', 'mact', 'minit', 'age', 'grav', 'feh', 'teff', 'lum',
          'parentid', 'partid', 'partitionid'}
POSTPROCESS = {'ra', 'dec', 'glon', 'glat', 'pi', 'mura', 'mudec', 'mul', 'mub', 'vr'}


def test_the_keys_that_the_backend_gives():
    assert set(Output._export_keys) == EXPORT and len(Output._export_keys) == len(EXPORT)


def test_the_keys_that_are_made_afterwards():
    assert set(Output._postprocess_keys) == POSTPROCESS and len(Output._postprocess_keys) == len(POSTPROCESS)


def test_the_optional_keys_of_the_input_that_a_catalogue_can_have():
    keys = {key for key, description in Output._all_optional_properties}
    assert keys == {'dform', 'helium', 'carbon', 'nitrogen', 'oxygen', 'neon', 'magnesium', 'silicon', 'sulphur', 'calcium', 'alpha',
                    'satid'}                                                                     # satid is the id of the input
    assert Output._satindex == 'satid'


def test_the_names_of_the_columns_of_a_catalogue():
    assert (Output._pos, Output._vel, Output._cel, Output._gal) == (('px', 'py', 'pz'), ('vx', 'vy', 'vz'), ('ra', 'dec'), ('glon', 'glat'))
    assert (Output._rad, Output._dmod, Output._mtip, Output._mact, Output._minit) == ('rad', 'dmod', 'mtip', 'mact', 'minit')
    assert (Output._age, Output._grav, Output._feh, Output._teff, Output._lum) == ('age', 'grav', 'feh', 'teff', 'lum')
    assert (Output._parentid, Output._partitionid, Output._partid) == ('parentid', 'partitionid', 'partid')
    assert (Output._pi, Output._mu, Output._mugal, Output._vr) == ('pi', ('mura', 'mudec'), ('mul', 'mub'), 'vr')


def test_the_magnitudes_of_the_photometric_systems_are_exported_in_order():
    assert Output._compile_export_mag_names([GAIA, LSST]) == ('gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp', 'lsst_u', 'lsst_g')
    assert Output._compile_export_mag_names([]) == ()


def test_the_export_keys_are_those_of_the_backend_the_extra_ones_and_the_magnitudes():
    keys = Output._make_export_keys([GAIA], extra_keys=('dform', 'satid'))
    assert set(keys) == EXPORT | {'dform', 'satid', 'gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp'} and len(keys) == len(set(keys))


def test_the_catalogue_keys_are_the_export_keys_and_those_that_are_made_afterwards():
    keys = Output._make_catalogue_keys([GAIA, LSST], extra_keys=('dform',))
    assert set(keys) == EXPORT | POSTPROCESS | {'dform', 'gaia__dr2_g', 'gaia__dr2_gbp', 'gaia__dr2_grp', 'lsst_u', 'lsst_g'}
    assert len(keys) == len(set(keys))


def test_a_key_that_is_given_twice_is_only_there_once():
    assert sorted(Output._make_export_keys([GAIA], extra_keys=('px', 'dform', 'dform'))).count('px') == 1


# -- what an Output is made of ------------------------------------------------------

def test_what_an_output_gives_to_vaex():
    assert Output._vaex_under_list == ['_repr_html_']


def test_the_documentation_of_the_constructor_has_what_it_is_made_of_filled_in():
    documentation = Output.__init__.__doc__
    assert "{_output_properties}" not in documentation and "{_optional_properties}" not in documentation
    assert "{STARCATALOG_GROUP}" not in documentation and "{PARENTPARTICLE_GROUP}" not in documentation
    assert "via keys ``(px, py, pz)``" in documentation and "via key ``rad``" in documentation        # the ones that are made of three
    assert "via keys ``(glon, glat)``" in documentation and "via key ``satid``" in documentation
    assert f"``{STARCATALOG_GROUP}``" in documentation and f"``{PARENTPARTICLE_GROUP}``" in documentation


if __name__ == '__main__':
    pass
