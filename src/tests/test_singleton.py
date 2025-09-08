#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
#
# This file is part of py-Galaxia-ananke:
# <https://github.com/athob/py-Galaxia-ananke>.
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


if __name__ == '__main__':
    pass
