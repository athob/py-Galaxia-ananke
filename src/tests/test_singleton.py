#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
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
from ..galaxia_ananke._builtin_utils import Singleton


def test_singleton():
    assert isinstance(Singleton, type)
    class IsSingletonClass(metaclass=Singleton):
        pass
    temp1 = IsSingletonClass()
    temp2 = IsSingletonClass()
    assert temp1 == temp2


def test_a_singleton_class_has_one_instance():
    class Alone(metaclass=Singleton):
        pass
    assert Alone() is Alone()


def test_the_first_call_decides_and_the_later_arguments_are_ignored():
    class WithValue(metaclass=Singleton):
        def __init__(self, value):
            self.value = value
    first = WithValue(1)
    assert WithValue(2) is first
    assert first.value == 1


def test_each_singleton_class_has_its_own_instance():
    class One(metaclass=Singleton):
        pass
    class Other(metaclass=Singleton):
        pass
    assert One() is not Other()
    assert isinstance(One(), One) and isinstance(Other(), Other)


def test_a_subclass_is_another_singleton_than_its_parent():
    class Parent(metaclass=Singleton):
        pass
    class Child(Parent):
        pass
    assert Child() is not Parent()
    assert Child() is Child()


if __name__ == '__main__':
    pass
