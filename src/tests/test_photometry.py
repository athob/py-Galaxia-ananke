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
Tests of Photometry, the registry of the photometric systems that a survey can be given (a dictionary of the categories, each a
dictionary of systems, that is read from the installation once), and of the dictionary that it is made of, which takes names
with a slash. The registry of the installation is a singleton and is not touched: the ones that are tested are read again from
temporary directories, and the paths of the isochrones go to them too.
"""
import sys

import pytest

from ..galaxia_ananke import photometry
from ..galaxia_ananke._constants import CUSTOM_PHOTOCAT, LEGACY_PHOTOCAT
from ..galaxia_ananke.photometry.Isochrone import Isochrone
from ..galaxia_ananke.photometry.PhotoSystem import PhotoSystem
from ..galaxia_ananke.photometry.Photometry import Photometry, nested_dict
from .utils import custom_isochrone_data

PHOTOMETRY_MODULE = sys.modules[Photometry.__module__]
ISOCHRONE_MODULE = sys.modules[Isochrone.__module__]


@pytest.fixture
def isochrones_path(tmp_path, monkeypatch):
    """The directory where the registry looks for isochrones, and where the isochrones that are made go: a temporary one."""
    path = tmp_path / 'Isochrones'
    monkeypatch.setattr(PHOTOMETRY_MODULE, 'ISOCHRONES_PATH', path)
    monkeypatch.setattr(ISOCHRONE_MODULE, 'ISOCHRONES_PATH', path)
    monkeypatch.setattr(Photometry, 'instance', Photometry.instance)      # the registry that is read again takes this place: put back
    return path


def make_system(path, category, name, descriptor=True):
    directory = path / category / name
    directory.mkdir(parents=True)
    if descriptor:
        (directory / Isochrone._file_descriptor).write_text(f"{name} 8 6 2 G BP\n\n")
    return directory


def scan():
    """The registry read again from the directory, as it was when the package was imported (a singleton would give the old one)."""
    return Photometry.__new__(Photometry)


# -- the registry of the installation -----------------------------------------------

def test_the_registry_is_a_singleton_and_it_is_the_one_of_the_package():
    assert Photometry() is Photometry() and photometry.available_photo_systems is Photometry()


def test_the_registry_of_the_installation_has_the_system_that_surveys_use_by_default():
    system = photometry.available_photo_systems['padova/GAIA__DR2']
    assert isinstance(system, PhotoSystem) and system.key == 'padova/GAIA__DR2'
    assert {'G', 'Gbp', 'Grp'} <= set(system.mag_names)
    assert CUSTOM_PHOTOCAT in photometry.available_photo_systems


# -- what the registry is made of ---------------------------------------------------

def test_there_is_a_category_for_each_directory_and_a_system_for_each_isochrone_that_has_a_descriptor(isochrones_path):
    for category, name in (('padova', 'GAIA__DR2'), ('padova', 'LSST'), ('other', 'Thing'), (CUSTOM_PHOTOCAT, 'Mine')):
        make_system(isochrones_path, category, name)
    registry = scan()
    assert sorted(registry) == ['other', 'padova', CUSTOM_PHOTOCAT]
    assert sorted(registry['padova']) == ['GAIA__DR2', 'LSST'] and sorted(registry['other']) == ['Thing']
    assert all(isinstance(system, PhotoSystem) for category in registry.values() for system in category.values())
    assert registry['padova/LSST'] is registry['padova']['LSST'] and registry['padova/LSST'].key == 'padova/LSST'
    assert list(registry[CUSTOM_PHOTOCAT]) == ['Mine']


def test_an_isochrone_without_a_descriptor_is_not_a_photometric_system(isochrones_path):
    make_system(isochrones_path, 'padova', 'LSST')
    make_system(isochrones_path, 'padova', 'Empty', descriptor=False)
    assert list(scan()['padova']) == ['LSST']


def test_a_category_without_any_system_is_there_and_empty(isochrones_path):
    make_system(isochrones_path, 'other', 'Nothing', descriptor=False)
    assert scan()['other'] == {}


def test_the_tables_of_the_bolometric_corrections_are_not_systems(isochrones_path):
    make_system(isochrones_path, 'BolCorr', 'Table')
    make_system(isochrones_path, 'padova', 'LSST')
    assert 'BolCorr' not in scan()


def test_what_is_not_a_directory_of_a_directory_is_not_looked_at(isochrones_path):
    make_system(isochrones_path, 'padova', 'LSST')
    (isochrones_path / 'padova' / 'notes.txt').write_text('')
    (isochrones_path / 'readme.txt').write_text('')
    registry = scan()
    assert sorted(registry) == ['padova', CUSTOM_PHOTOCAT] and list(registry['padova']) == ['LSST']


def test_the_custom_category_is_there_even_when_nothing_is_in_it(isochrones_path):
    make_system(isochrones_path, 'padova', 'LSST')
    assert scan()[CUSTOM_PHOTOCAT] == {}


def test_without_any_isochrone_the_registry_is_the_empty_custom_category(isochrones_path):
    assert not isochrones_path.exists()
    registry = scan()
    assert list(registry) == [CUSTOM_PHOTOCAT] and registry[CUSTOM_PHOTOCAT] == {}


# -- changing the registry ----------------------------------------------------------

def test_the_padova_category_is_protected_once_it_is_there(isochrones_path):
    make_system(isochrones_path, LEGACY_PHOTOCAT, 'LSST')
    registry = scan()
    with pytest.raises(ValueError, match=f"The '{LEGACY_PHOTOCAT}' entry is protected"):
        registry[LEGACY_PHOTOCAT] = {}


def test_the_other_categories_can_be_replaced(isochrones_path):
    make_system(isochrones_path, LEGACY_PHOTOCAT, 'LSST')
    registry = scan()
    registry['other'] = {'one': 1}
    registry[CUSTOM_PHOTOCAT] = {'two': 2}
    assert registry['other'] == {'one': 1} and registry[CUSTOM_PHOTOCAT] == {'two': 2}


def test_the_padova_category_can_be_set_when_it_is_not_there(isochrones_path):
    registry = scan()
    registry[LEGACY_PHOTOCAT] = {'one': 1}
    assert registry[LEGACY_PHOTOCAT] == {'one': 1}


@pytest.mark.parametrize("key", [''.join(['pad', 'ova']), f"{LEGACY_PHOTOCAT}/Mine"], ids=["a name made another way", "a path"])
def test_the_padova_category_is_protected_whatever_the_way_it_is_named(key, isochrones_path):
    make_system(isochrones_path, LEGACY_PHOTOCAT, 'LSST')
    registry = scan()
    with pytest.raises(ValueError, match="is protected"):
        registry[key] = {}


@pytest.mark.parametrize("key", [f"{LEGACY_PHOTOCAT}_two", f"two_{LEGACY_PHOTOCAT}", f"{LEGACY_PHOTOCAT}_two/Mine"])
def test_a_category_that_is_only_named_like_the_padova_one_is_not_protected(key, isochrones_path):
    make_system(isochrones_path, LEGACY_PHOTOCAT, 'LSST')
    registry = scan()
    registry[key] = {'one': 1}                                # is not refused,
    assert key.split('/')[0] in registry                      # and its category is there


def test_nothing_can_be_deleted_from_the_registry(isochrones_path):
    make_system(isochrones_path, 'other', 'Thing')
    with pytest.raises(NotImplementedError):
        del scan()['other']


# -- the isochrones that are added --------------------------------------------------

def test_an_isochrone_is_added_in_the_custom_category_and_written_there(isochrones_path):
    registry = scan()
    registry.add_isochrone('mine', custom_isochrone_data())
    assert list(registry[CUSTOM_PHOTOCAT]) == ['mine']
    assert (isochrones_path / CUSTOM_PHOTOCAT / 'mine' / 'IsoFileDescriptor.txt').is_file()


def test_a_name_with_a_slash_cannot_be_added(isochrones_path):
    with pytest.raises(ValueError):
        scan().add_isochrone('a/b', custom_isochrone_data())
    assert not (isochrones_path / CUSTOM_PHOTOCAT).exists()


def test_an_isochrone_that_was_added_can_only_be_added_again_over_the_old_one(isochrones_path):
    registry = scan()
    registry.add_isochrone('mine', custom_isochrone_data())
    with pytest.raises(FileExistsError):
        registry.add_isochrone('mine', custom_isochrone_data())
    registry.add_isochrone('mine', custom_isochrone_data(), overwrite=True)
    assert list(registry[CUSTOM_PHOTOCAT]) == ['mine']


def test_an_isochrone_that_was_added_is_a_photometric_system_like_the_other_entries(isochrones_path):
    registry = scan()
    registry.add_isochrone('mine', custom_isochrone_data())
    system = registry[CUSTOM_PHOTOCAT]['mine']
    assert isinstance(system, PhotoSystem) and (system.category, system.name) == (CUSTOM_PHOTOCAT, 'mine')
    assert registry[f"{CUSTOM_PHOTOCAT}/mine"] is system


# -- the dictionary with slashes ----------------------------------------------------

def test_a_name_with_slashes_goes_down_the_dictionaries():
    nested = nested_dict({'a': {'b': {'c': 1}}, 'd': 2})
    assert nested['a/b/c'] == 1 and nested['a/b'] == {'c': 1} and nested['d'] == 2


@pytest.mark.parametrize("key", ['x', 'a/x', 'a/b/x'])
def test_a_name_that_is_not_there_is_a_key_error(key):
    with pytest.raises(KeyError):
        nested_dict({'a': {'b': {'c': 1}}})[key]


def test_a_name_without_a_slash_is_set_as_in_any_dictionary():
    nested = nested_dict()
    nested['a'] = {'b': 1}
    assert nested == {'a': {'b': 1}}


@pytest.mark.parametrize("start", [{}, {'a': {'x': 1}}], ids=["a new category", "an existing category"])
def test_a_name_with_a_slash_can_be_set(start):
    nested = nested_dict(start)
    nested['a/b'] = 2
    assert nested['a/b'] == 2 and nested['a'] == {**start.get('a', {}), 'b': 2}


def test_a_name_with_several_slashes_makes_the_dictionaries_on_the_way():
    nested = nested_dict({'d': 2})
    nested['a/b/c'] = 3
    nested['a/b/e'] = 4
    assert nested == {'d': 2, 'a': {'b': {'c': 3, 'e': 4}}} and nested['a/b/c'] == 3


def test_a_name_with_a_slash_is_set_in_the_dictionary_that_is_there_and_not_in_a_copy_of_it():
    inner = {'x': 1}
    nested = nested_dict({'a': inner})
    nested['a/b'] = 2
    assert inner == {'x': 1, 'b': 2} and nested['a'] is inner


def test_a_system_can_be_set_in_a_category_of_the_registry_by_its_path(isochrones_path):
    registry = scan()
    registry[f"{CUSTOM_PHOTOCAT}/Mine"] = 'a system'
    assert registry[CUSTOM_PHOTOCAT] == {'Mine': 'a system'} and registry[f"{CUSTOM_PHOTOCAT}/Mine"] == 'a system'


# -- clearing the isochrones --------------------------------------------------------

def make_installation(isochrones_path):
    make_system(isochrones_path, 'padova', 'LSST')
    make_system(isochrones_path, CUSTOM_PHOTOCAT, 'Mine')


def test_nothing_is_removed_unless_it_is_asked_hard(isochrones_path):
    make_installation(isochrones_path)
    scan()._purge()
    scan()._purge(hard=False, yes=True)
    assert sorted(path.name for path in isochrones_path.iterdir()) == ['padova', CUSTOM_PHOTOCAT]


def test_everything_is_removed_when_it_is_asked_hard_and_yes(isochrones_path):
    make_installation(isochrones_path)
    scan()._purge(hard=True, yes=True)
    assert isochrones_path.is_dir() and list(isochrones_path.iterdir()) == []


@pytest.mark.parametrize("answer, removed", [('y', True), ('yes', True), ('YES', True), ('n', False), ('', False), ('maybe', False)])
def test_what_is_removed_when_it_is_asked_hard_depends_on_the_answer(answer, removed, isochrones_path, monkeypatch):
    make_installation(isochrones_path)
    monkeypatch.setattr('builtins.input', lambda prompt: answer)
    scan()._purge(hard=True)
    assert (list(isochrones_path.iterdir()) == []) is removed


if __name__ == '__main__':
    pass
