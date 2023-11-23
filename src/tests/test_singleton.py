#!/usr/bin/env python
from ..Galaxia_ananke.utils import Singleton


def test_singleton():
    assert isinstance(Singleton, type)
    class IsSingletonClass(metaclass=Singleton):
        pass
    temp1 = IsSingletonClass()
    temp2 = IsSingletonClass()
    assert temp1 == temp2


if __name__ == '__main__':
    pass
