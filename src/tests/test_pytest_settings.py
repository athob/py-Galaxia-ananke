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
Tests of the settings of pytest in pyproject.toml (the ``filterwarnings`` entries). They only need the settings, which are in
force while a test runs, and pytest shows a DeprecationWarning that they do not ignore.
"""
import warnings


def test_the_settings_ignore_the_fork_warning_of_multiprocess_and_that_of_no_other_module():
    message = "This process (pid=1) is multi-threaded, use of fork() may lead to deadlocks in the child."
    with warnings.catch_warnings(record=True) as caught:
        for filename in ('multiprocess/popen_fork.py', 'another/module.py'):
            module = filename[:-3].replace('/', '.')
            warnings.warn_explicit(message, DeprecationWarning, filename, 1, module=module, registry={})
    assert [w.filename for w in caught] == ['another/module.py']
