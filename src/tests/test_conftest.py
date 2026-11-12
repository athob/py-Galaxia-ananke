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
Tests of conftest.py, the seeding that every test starts from. They run pytest, with that very file, on a project
of two tests that print what they draw, so that what is checked is what a run does and not a copy of the logic.
"""
import importlib.util
import os
import pathlib
import re
import subprocess
import sys
import textwrap

CONFTEST = pathlib.Path(__file__).with_name("conftest.py")

DRAWS = textwrap.dedent('''
    import random

    import numpy as np


    def test_first():
        print("DRAW first", random.random(), np.random.random())


    def test_second():
        print("DRAW second", random.random(), np.random.random())


    def test_fails():
        assert False
''')


def run_pytest(tmp_path, *selection, shift=None):
    """Run pytest, with the conftest.py of the suite, on a project made of DRAWS: (the output, the return code)."""
    (tmp_path / "conftest.py").write_text(CONFTEST.read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "test_draws.py").write_text(DRAWS, encoding="utf-8")
    env = {key: value for key, value in os.environ.items() if key != "GALAXIA_TEST_SEED"}
    if shift is not None:
        env["GALAXIA_TEST_SEED"] = str(shift)
    done = subprocess.run([sys.executable, "-m", "pytest", "-s", "-q", "-p", "no:cacheprovider", f"--rootdir={tmp_path}",
                           *selection], capture_output=True, text=True, env=env, cwd=tmp_path)
    return done.stdout + done.stderr, done.returncode


def draws(tmp_path, shift=None):
    """What the two tests drew: {'first': (python, numpy), 'second': (python, numpy)}."""
    output, code = run_pytest(tmp_path, "test_draws.py::test_first", "test_draws.py::test_second", shift=shift)
    assert code == 0, output
    return {name: (python, numpy) for name, python, numpy in re.findall(r"DRAW (\w+) (\S+) (\S+)", output)}


def test_a_test_draws_the_same_values_in_every_run(tmp_path):
    first_run = draws(tmp_path)
    assert len(first_run) == 2
    assert draws(tmp_path) == first_run


def test_two_tests_do_not_draw_the_same_values(tmp_path):
    drawn = draws(tmp_path)
    assert drawn["first"][0] != drawn["second"][0], "Python's random state is not seeded apart for each test"
    assert drawn["first"][1] != drawn["second"][1], "the random state of NumPy is not seeded apart for each test"


def test_the_seed_can_be_shifted(tmp_path):
    default, shifted = draws(tmp_path), draws(tmp_path, shift=1)
    assert shifted != default
    assert shifted == draws(tmp_path, shift=1), "a shifted run is not reproducible either"


def test_a_failing_test_reports_its_seed(tmp_path):
    output, code = run_pytest(tmp_path, "test_draws.py::test_fails")
    assert code == 1
    assert "random seed" in output


# -- the files that a test leaves where pytest runs from -------------------------------

LITTER = textwrap.dedent('''
    import pathlib


    def test_leaves_a_file():
        pathlib.Path("left_behind.txt").write_text("")


    def test_leaves_the_same_file_too():
        pathlib.Path("left_behind.txt").write_text("")


    def test_works_in_its_own_directory(tmp_path):
        (tmp_path / "mine.txt").write_text("")


    def test_leaves_what_is_hidden():
        pathlib.Path(".hidden").write_text("")
''')


def run_litter(tmp_path, *selection):
    """Run pytest, with the conftest.py of the suite, on a project made of LITTER, from its directory: (the output, the return code)."""
    (tmp_path / "conftest.py").write_text(CONFTEST.read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "test_litter.py").write_text(LITTER, encoding="utf-8")
    done = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", f"--rootdir={tmp_path}", *selection],
                          capture_output=True, text=True, cwd=tmp_path)
    return done.stdout + done.stderr, done.returncode


def load_conftest():
    spec = importlib.util.spec_from_file_location("conftest_under_test", CONFTEST)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_a_test_that_leaves_a_file_where_pytest_runs_from_fails_with_its_name_and_the_file(tmp_path):
    output, code = run_litter(tmp_path)
    assert code != 0 and "2 errors" in output
    for name in ("test_leaves_a_file", "test_leaves_the_same_file_too"):          # the second one is reported although the file is the same
        assert f"ERROR at teardown of {name}" in output
    assert "left_behind.txt" in output and "where pytest runs from" in output
    for name in ("test_works_in_its_own_directory", "test_leaves_what_is_hidden"):
        assert f"ERROR at teardown of {name}" not in output


def test_the_file_that_a_test_left_is_removed_but_not_what_is_hidden(tmp_path):
    run_litter(tmp_path)
    assert not (tmp_path / "left_behind.txt").exists() and (tmp_path / ".hidden").exists()


def test_the_tests_that_work_in_their_own_directory_or_leave_what_is_hidden_pass(tmp_path):
    output, code = run_litter(tmp_path, "-k", "own_directory or hidden")
    assert code == 0 and "2 passed" in output


def test_what_is_left_in_a_directory_leaves_out_what_was_there_what_is_hidden_and_the_caches_of_python(tmp_path):
    (tmp_path / "there.txt").write_text("")
    before = {path.name for path in tmp_path.iterdir()}
    for name in ("new.txt", ".hidden", ".coverage.host.1.abc"):
        (tmp_path / name).write_text("")
    (tmp_path / "__pycache__").mkdir()
    (tmp_path / "new_directory").mkdir()
    assert [path.name for path in load_conftest().left_in(tmp_path, before)] == ["new.txt", "new_directory"]
