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
Tests of State, which tracks the steps of a survey and, when caching, skips the ones that a previous run did.
"""
import json
import types

import pytest

from ..galaxia_ananke._builtin_utils import State, _BaseState


def make_state(tmp_path, caching=True, name="survey"):
    return State(tmp_path / name, types.SimpleNamespace(caching=caching, verbose=False))


# -- the state itself ---------------------------------------------------------------

def test_a_new_state_is_the_initial_one():
    state = _BaseState()
    assert (state.register, state.description) == ([0], "initial")


def test_new_states_do_not_share_their_register():
    first, second = _BaseState(), _BaseState()
    first.update()
    assert (first.register, second.register) == ([1], [0])


@pytest.mark.parametrize("levels, register", [
    ([0], [1]),
    ([0, 0], [2]),
    ([0, 1], [1, 1]),
    ([0, 1, 1], [1, 2]),
    ([0, 1, 1, 0], [2]),            # a new step at a level forgets what was counted below it
    ([2], [0, 0, 1]),               # a level that does not exist yet is made, with zeros above it
    ([0, 2, 1], [1, 1]),
])
def test_update_counts_a_step_at_its_level(levels, register):
    state = _BaseState()
    for level in levels:
        state.update(level=level)
    assert state.register == register


def test_update_sets_the_description():
    state = _BaseState()
    state.update("exporting")
    assert state.description == "exporting"
    state.update()
    assert state.description == ""


def test_a_state_goes_through_a_dictionary_and_back():
    state = _BaseState(register=[3, 1], description="exporting")
    assert state.asdict() == {"register": [3, 1], "description": "exporting"}
    assert _BaseState.fromdict(state.asdict()) == state


@pytest.mark.parametrize("smaller, larger", [
    (([0], "a"), ([1], "a")),
    (([1], "a"), ([1, 0], "a")),    # a longer register is further
    (([1, 2], "a"), ([2], "a")),
    (([1], "a"), ([1], "b")),       # at the same register, the description decides
])
def test_states_are_ordered_by_register_then_by_description(smaller, larger):
    small, large = _BaseState(*smaller), _BaseState(*larger)
    assert small < large and small <= large and large > small and large >= small
    assert small != large and not large < small


def test_equal_states_are_equal_and_neither_smaller_nor_larger():
    one, other = _BaseState([1], "a"), _BaseState([1], "a")
    assert one == other and one <= other and one >= other
    assert not one < other and not one > other


# -- the state and its file ---------------------------------------------------------

def test_a_state_makes_its_file_with_the_initial_state(tmp_path):
    state = make_state(tmp_path)
    assert state.file_path == (tmp_path / "survey.state").resolve()
    assert json.loads(state.file_path.read_bytes()) == {"register": [0], "description": "initial"}


def test_the_state_extension_replaces_the_last_one(tmp_path):
    # Output ends the base of its files with the number of the partition: that is what .state replaces
    assert State(tmp_path / "a.b.0", types.SimpleNamespace(caching=True, verbose=False)).file_path.name == "a.b.state"


def test_a_state_does_not_overwrite_the_file_it_finds(tmp_path):
    (tmp_path / "survey.state").write_text(json.dumps({"register": [5], "description": "later"}))
    state = make_state(tmp_path)
    assert json.loads(state.file_path.read_text()) == {"register": [5], "description": "later"}


def test_caching_and_verbose_are_those_of_the_parent(tmp_path):
    state = State(tmp_path / "survey", types.SimpleNamespace(caching="the caching", verbose="the verbose"))
    assert (state.caching, state.verbose) == ("the caching", "the verbose")


def test_without_caching_a_step_always_runs_and_is_recorded(tmp_path):
    calls = []
    for _ in range(2):
        state = make_state(tmp_path, caching=False)
        state.check_state_file_before_running("step")(lambda: calls.append(1))()
    assert len(calls) == 2
    assert json.loads(state.file_path.read_text()) == {"register": [1], "description": "step"}


def test_with_caching_a_step_that_was_done_is_skipped(tmp_path):
    calls = []
    for _ in range(2):
        state = make_state(tmp_path, caching=True)
        state.check_state_file_before_running("step")(lambda: calls.append(1))()
    assert calls == [1]


def test_with_caching_only_the_steps_after_the_last_one_done_run(tmp_path):
    done = []

    def pipeline(steps):
        state = make_state(tmp_path, caching=True)
        for step in steps:
            state.check_state_file_before_running(step)(lambda step=step: done.append(step))()

    pipeline(["a", "b"])
    pipeline(["a", "b", "c"])
    assert done == ["a", "b", "c"]


def test_with_caching_a_state_file_that_vanished_means_that_nothing_was_done(tmp_path):
    calls = []
    state = make_state(tmp_path, caching=True)
    state.file_path.unlink()
    state.check_state_file_before_running("step")(lambda: calls.append(1))()
    assert calls == [1]
    assert state.file_path.exists()                                  # and the step is recorded again


def test_a_step_at_another_level_is_recorded_at_that_level(tmp_path):
    state = make_state(tmp_path, caching=False)
    state.check_state_file_before_running("step", level=1)(lambda: None)()
    assert json.loads(state.file_path.read_text()) == {"register": [0, 1], "description": "step"}


def test_the_wrapper_gives_its_arguments_to_the_step_and_returns_nothing(tmp_path):
    seen = []
    state = make_state(tmp_path, caching=False)
    wrapped = state.check_state_file_before_running("step")(lambda *args, **kwargs: seen.append((args, kwargs)) or "ignored")
    assert wrapped(1, 2, key=3) is None
    assert seen == [((1, 2), {"key": 3})]


if __name__ == '__main__':
    pass
