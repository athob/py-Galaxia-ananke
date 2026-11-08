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
Tests of the marking of the properties that go into the metadata of a class, mark_metadata_prop and
collect_metadata_marked_properties, which Input, Survey and Output use to write what they are in their files.
"""
import pytest

from ..galaxia_ananke._builtin_utils import METADATA_MARK, collect_metadata_marked_properties, mark_metadata_prop


@collect_metadata_marked_properties
class Thing:
    def __init__(self, base):
        self.base = base

    @property
    @mark_metadata_prop
    def single(self):               # marked before @property
        return self.base

    @mark_metadata_prop
    @property
    def double(self):               # marked after @property (and defined after single, which is after it in lexical order)
        return 2 * self.base

    @property
    def hidden(self):               # a property that is not marked
        return -1

    plain = 5                       # an attribute that is not a property


# -- the marking --------------------------------------------------------------------

def test_a_function_is_marked_and_returned_itself():
    def getter(self):
        return 1
    assert mark_metadata_prop(getter) is getter
    assert getattr(getter, METADATA_MARK) is True


def test_a_property_is_returned_itself_and_its_getter_is_marked():
    prop = property(lambda self: 1)
    assert mark_metadata_prop(prop) is prop
    assert getattr(prop.fget, METADATA_MARK) is True


def test_a_property_without_a_getter_is_accepted_and_marks_nothing():
    prop = property()
    assert mark_metadata_prop(prop) is prop


def test_both_orders_of_the_decorators_mark_the_property():
    assert getattr(Thing.double.fget, METADATA_MARK) is True
    assert getattr(Thing.single.fget, METADATA_MARK) is True
    assert not hasattr(Thing.hidden.fget, METADATA_MARK)


# -- the collection -----------------------------------------------------------------

def test_the_names_are_those_of_the_marked_properties_in_lexical_order():
    assert Thing(1)._metadata_names == ['double', 'single']


def test_a_getter_whose_mark_is_false_is_not_collected():
    @collect_metadata_marked_properties
    class Off:
        @mark_metadata_prop
        @property
        def value(self):
            return 1
    setattr(Off.value.fget, METADATA_MARK, False)
    assert Off()._metadata_names == []


def test_the_metadata_holds_the_current_value_of_each_marked_property():
    thing = Thing(3)
    assert thing._metadata == {'double': 6, 'single': 3}
    thing.base = 10
    assert thing._metadata == {'double': 20, 'single': 10}


def test_two_instances_have_their_own_values():
    assert Thing(1)._metadata != Thing(2)._metadata


def test_the_names_are_made_once_per_instance():
    thing = Thing(1)
    assert thing._metadata_names is thing._metadata_names


def test_the_decorator_returns_the_class_itself():
    class Plain:
        pass
    assert collect_metadata_marked_properties(Plain) is Plain


def test_a_class_without_marked_properties_has_empty_metadata():
    @collect_metadata_marked_properties
    class Nothing:
        @property
        def value(self):
            return 1
    assert (Nothing()._metadata_names, Nothing()._metadata) == ([], {})


def test_a_property_without_a_getter_is_not_in_the_metadata():
    @collect_metadata_marked_properties
    class Odd:
        broken = property()
        @mark_metadata_prop
        @property
        def fine(self):
            return 1
    assert Odd()._metadata == {'fine': 1}


def test_a_subclass_has_the_marked_properties_of_its_own_class_only():
    # How it is today: the names come from the class of the instance, not from its parents. A subclass that is
    # to have its own metadata is decorated as well (and then still does not have the ones of its parent).
    class Child(Thing):
        @property
        @mark_metadata_prop
        def extra(self):
            return 7
    assert Child(1)._metadata_names == ['extra']
    assert Child(1)._metadata == {'extra': 7}


if __name__ == '__main__':
    pass
