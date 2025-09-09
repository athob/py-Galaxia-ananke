#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 (AGPL-3.0).
# 
# The full copyright notice, including terms governing use, modification,
# and redistribution, is contained in the files LICENSE and COPYRIGHT,
# which can be found at the root of the source code distribution tree:
# - LICENSE <https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE>
# - COPYRIGHT <https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT>
#
import pytest

from ..galaxia_ananke._builtin_utils import compare_given_and_required


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
