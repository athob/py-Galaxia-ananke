#!/usr/bin/env python
import pytest

from ..Galaxia_ananke.utils import compare_given_and_required


def test():
    given = ['A','B','C','D']
    required = ['A','C']
    optional = ['B','D']
    extra = ['E']
    compare_given_and_required(given, required, optional)
    compare_given_and_required(given, required, optional+extra)
    with pytest.raises(ValueError):
        compare_given_and_required(given)
    with pytest.raises(ValueError):
        compare_given_and_required(given, required)
    with pytest.raises(ValueError):
        compare_given_and_required(given, required+extra, optional)
    with pytest.raises(ValueError):
        compare_given_and_required(given+extra, required, optional)


if __name__ == '__main__':
    pass
