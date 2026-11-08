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
import os
import pathlib

import pytest

from ..galaxia_ananke._builtin_utils import make_symlink
from .utils import in_tmp_wd


@in_tmp_wd
def test_a_symlink_to_the_file_is_made_in_the_directory():
    file_path = 'file_path'
    file_content = "test"
    pathlib.Path(file_path).write_text(file_content)
    dest_dir = 'dest_dir'
    pathlib.Path(dest_dir).mkdir()
    make_symlink(file_path, dest_dir)
    symlink = pathlib.Path(dest_dir)/file_path
    assert symlink.is_symlink()
    assert symlink.resolve() == pathlib.Path(file_path).absolute()


def test_the_symlink_has_the_name_of_the_file_and_points_to_its_absolute_path(tmp_path):
    (tmp_path / "source").mkdir()
    (tmp_path / "dest").mkdir()
    file_path = tmp_path / "source" / "particles.ebf"
    file_path.write_text("content")
    make_symlink(file_path, tmp_path / "dest")
    symlink = tmp_path / "dest" / "particles.ebf"
    assert os.readlink(symlink) == str(file_path.resolve())
    assert symlink.read_text() == "content"


def test_relative_paths_are_taken_from_where_the_call_is_made(tmp_path, monkeypatch):
    (tmp_path / "dest").mkdir()
    (tmp_path / "file").write_text("content")
    monkeypatch.chdir(tmp_path)
    make_symlink("file", "dest")
    assert os.readlink(tmp_path / "dest" / "file") == str((tmp_path / "file").resolve())


def test_a_symlink_that_is_already_there_is_replaced(tmp_path):
    (tmp_path / "dest").mkdir()
    old, new = tmp_path / "old" / "file", tmp_path / "new" / "file"
    for path, content in ((old, "old"), (new, "new")):
        path.parent.mkdir()
        path.write_text(content)
    make_symlink(old, tmp_path / "dest")
    make_symlink(new, tmp_path / "dest")
    assert (tmp_path / "dest" / "file").read_text() == "new"


def test_a_broken_symlink_that_is_already_there_is_replaced(tmp_path):
    (tmp_path / "dest").mkdir()
    (tmp_path / "dest" / "file").symlink_to(tmp_path / "nothing")
    (tmp_path / "file").write_text("content")
    make_symlink(tmp_path / "file", tmp_path / "dest")
    assert (tmp_path / "dest" / "file").read_text() == "content"


def test_a_regular_file_that_is_already_there_is_replaced_by_the_symlink(tmp_path):
    (tmp_path / "dest").mkdir()
    (tmp_path / "dest" / "file").write_text("a copy")
    (tmp_path / "file").write_text("the file")
    make_symlink(tmp_path / "file", tmp_path / "dest")
    assert (tmp_path / "dest" / "file").is_symlink()
    assert (tmp_path / "dest" / "file").read_text() == "the file"


def test_a_symlink_of_a_symlink_points_to_the_file_itself(tmp_path):
    (tmp_path / "dest").mkdir()
    (tmp_path / "file").write_text("content")
    (tmp_path / "link").symlink_to(tmp_path / "file")
    make_symlink(tmp_path / "link", tmp_path / "dest")
    assert os.readlink(tmp_path / "dest" / "file") == str((tmp_path / "file").resolve())


def test_a_directory_in_the_way_is_an_error_and_is_left_alone(tmp_path):
    (tmp_path / "dest" / "file").mkdir(parents=True)
    (tmp_path / "file").write_text("content")
    with pytest.raises((IsADirectoryError, PermissionError)):              # which one depends on the system
        make_symlink(tmp_path / "file", tmp_path / "dest")
    assert (tmp_path / "dest" / "file").is_dir() and not (tmp_path / "dest" / "file").is_symlink()


def test_a_directory_that_does_not_exist_is_a_file_not_found_error(tmp_path):
    (tmp_path / "file").write_text("content")
    with pytest.raises(FileNotFoundError):
        make_symlink(tmp_path / "file", tmp_path / "nowhere")


if __name__ == '__main__':
    pass
