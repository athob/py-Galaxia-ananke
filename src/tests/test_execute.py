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
Tests of execute, which runs the commands of the backend and prints what they say. Every command here is a small
Python script of its own, run by the interpreter that runs the tests.
"""
import subprocess
import sys
import textwrap

import pytest

from ..galaxia_ananke._builtin_utils import execute


def write_script(tmp_path, name, body):
    script = tmp_path / name
    script.write_text(textwrap.dedent(body))
    return script


@pytest.fixture
def say(tmp_path):
    """The command that prints the word it is given."""
    script = write_script(tmp_path, "say.py", "import sys\nprint(sys.argv[1])\n")
    return lambda word: f"{sys.executable} {script} {word}"


@pytest.fixture
def fail(tmp_path):
    """The command that prints a line and exits with the status it is given."""
    script = write_script(tmp_path, "fail.py", "import sys\nprint('before')\nsys.exit(int(sys.argv[1]))\n")
    return lambda status: f"{sys.executable} {script} {status}"


def test_a_command_is_announced_and_its_output_tagged(say, capsys):
    command = say("hello")
    execute(command)
    assert capsys.readouterr().out.splitlines() == [f"Executing JOB 1/1 = {command}", "JOB 1/1 | hello"]


def test_a_list_of_one_command_is_the_same_as_the_command(say, capsys):
    command = say("hello")
    execute(command)
    alone = capsys.readouterr().out
    execute([command])
    assert capsys.readouterr().out == alone


def test_what_a_command_writes_to_its_error_stream_is_printed_too(tmp_path, capsys):
    # the backend reports its own failures there (an uncaught C++ exception, for one)
    script = write_script(tmp_path, "oops.py", "import sys\nprint('oops', file=sys.stderr)\n")
    execute(f"{sys.executable} {script}")
    assert "JOB 1/1 | oops" in capsys.readouterr().out


def test_nothing_is_printed_when_not_verbose(say, capsys):
    execute(say("hello"), verbose=False)
    assert capsys.readouterr().out == ""


def test_several_commands_are_announced_first_and_tagged_each(say, capsys):
    commands = [say("one"), say("two")]
    execute(commands)
    lines = capsys.readouterr().out.splitlines()
    assert lines[:2] == [f"Executing JOB 1/2 = {commands[0]}", f"Executing JOB 2/2 = {commands[1]}"]
    assert sorted(lines[2:]) == ["JOB 1/2 | one", "JOB 2/2 | two"]


def test_the_tags_are_as_wide_as_the_number_of_commands_needs(say, capsys):
    commands = [say(f"word{i}") for i in range(10)]
    execute(commands)
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == f"Executing JOB  1/10 = {commands[0]}"
    assert lines[9] == f"Executing JOB 10/10 = {commands[9]}"
    assert "JOB  1/10 | word0" in lines and "JOB 10/10 | word9" in lines


def test_a_failing_command_raises_with_its_status_and_after_its_output(fail, capsys):
    command = fail(3)
    with pytest.raises(subprocess.CalledProcessError) as error:
        execute(command)
    assert (error.value.returncode, error.value.cmd) == (3, command)
    assert "JOB 1/1 | before" in capsys.readouterr().out


def test_the_failing_command_among_several_is_the_one_reported(say, fail, capsys):
    ok, failing = say("fine"), fail(5)
    with pytest.raises(subprocess.CalledProcessError) as error:
        execute([ok, failing])
    assert (error.value.returncode, error.value.cmd) == (5, failing)
    assert "JOB 1/2 | fine" in capsys.readouterr().out


def test_nothing_is_raised_when_the_commands_succeed_even_if_verbose_is_off(say):
    assert execute([say("a"), say("b")], verbose=False) is None


def test_the_other_arguments_go_to_popen(tmp_path, capsys):
    there = tmp_path / "there"
    there.mkdir()
    script = write_script(tmp_path, "where.py", "import os\nprint(os.getcwd())\n")
    execute(f"{sys.executable} {script}", cwd=there)
    assert capsys.readouterr().out.splitlines()[-1] == f"JOB 1/1 | {there.resolve()}"


def test_commands_run_at_the_same_time_whatever_max_workers_is(tmp_path):
    # Each command waits for the flag of the other: they only end well if both are running together, and
    # max_workers, which is not used, does not change that.
    meet = write_script(tmp_path, "meet.py", '''
        import pathlib, sys, time
        mine, other = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
        mine.touch()
        for _ in range(100):
            if other.exists():
                sys.exit(0)
            time.sleep(0.1)
        sys.exit(1)
    ''')
    a, b = tmp_path / "a.flag", tmp_path / "b.flag"
    execute([f"{sys.executable} {meet} {a} {b}", f"{sys.executable} {meet} {b} {a}"], max_workers=1, verbose=False)


@pytest.mark.xfail(raises=subprocess.CalledProcessError,
                   reason="execute splits a command on its spaces, so a path that holds one is cut in two")
def test_a_script_whose_path_holds_a_space_runs(tmp_path, capsys):
    folder = tmp_path / "a folder"
    folder.mkdir()
    script = write_script(folder, "say.py", "print('hello')\n")
    execute(f"{sys.executable} {script}")
    assert "JOB 1/1 | hello" in capsys.readouterr().out


if __name__ == '__main__':
    pass
