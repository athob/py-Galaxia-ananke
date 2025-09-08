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
import pathlib

from ..galaxia_ananke._builtin_utils import make_symlink
from .utils import in_tmp_wd


@in_tmp_wd
def test():
    file_path = 'file_path'
    file_content = "test"
    pathlib.Path(file_path).write_text(file_content)
    dest_dir = 'dest_dir'
    pathlib.Path(dest_dir).mkdir()
    make_symlink(file_path, dest_dir)
    symlink = pathlib.Path(dest_dir)/file_path
    assert symlink.is_symlink()
    assert symlink.resolve() == pathlib.Path(file_path).absolute()


if __name__ == '__main__':
    pass
