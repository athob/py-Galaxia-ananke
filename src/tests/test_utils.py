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
Tests of the utilities that need pandas: RecordingDataFrame, which Output uses to find out which columns of a survey a function
of the user reads, and hash_iterable, whose value names the files of the cache of an input.
"""
import warnings

import pandas as pd
import pytest

from ..galaxia_ananke.utils import RecordingDataFrame, hash_iterable


# -- hash_iterable ------------------------------------------------------------------
# The values are pinned: they name the files of the cache, so a change of them is a change of what a cache is, which
# leaves the caches of every earlier run unused.

@pytest.mark.parametrize("elements, expected", [
    ([], b"e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),               # the SHA-256 of nothing
    ([b"a"], b"da3811154d59c4267077ddd8bb768fa9b06399c486e1fc00485116b57c9872f5"),
    ([b"a", b"b"], b"29efd862964bb39dacc7a263c5db9f8187f6f7264ef9190b42aa8f3099bf9ad6"),
    ([b"b", b"a"], b"0d6b8562ada6532fbb83d290ef589552c0cdd3595596402f5d2ff2aa3368ffea"),
    ([b"ab"], b"6a866e4029d510139186c362675bb993a91e9132c867722faa519c5711ac1807"),
])
def test_the_hash_of_some_elements_is_what_it_always_was(elements, expected):
    assert hash_iterable(elements) == expected


def test_the_hash_is_the_hexadecimal_digest_as_bytes():
    digest = hash_iterable([b"anything"])
    assert isinstance(digest, bytes) and len(digest) == 64
    assert bytes.fromhex(digest.decode("ascii"))                 # it is hexadecimal


def test_the_order_of_the_elements_matters_and_so_do_their_boundaries():
    assert hash_iterable([b"a", b"b"]) != hash_iterable([b"b", b"a"])
    assert hash_iterable([b"a", b"b"]) != hash_iterable([b"ab"])


def test_any_iterable_of_buffers_can_be_hashed():
    assert (hash_iterable(element for element in [b"a", b"b"])
            == hash_iterable((b"a", b"b"))
            == hash_iterable([bytearray(b"a"), memoryview(b"b")])
            == hash_iterable([b"a", b"b"]))


def test_an_element_that_is_not_bytes_is_a_type_error():
    with pytest.raises(TypeError):
        hash_iterable(["text"])


# -- RecordingDataFrame -------------------------------------------------------------

@pytest.fixture
def frame():
    return RecordingDataFrame({"age": [1, 2], "feh": [3, 4], "mass": [5, 6]})


@pytest.mark.xfail(raises=UserWarning, reason="making a frame sets an attribute that pandas takes for a column that someone "
                   "forgot to make, and warns: whoever post-processes a survey with a function of their own sees it")
def test_a_recording_frame_is_made_without_a_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        RecordingDataFrame({"age": [1, 2]})


def test_a_new_frame_has_read_nothing(frame):
    assert frame.record_of_all_used_keys == set()


def test_a_frame_is_a_dataframe_that_gives_what_a_dataframe_gives(frame):
    assert isinstance(frame, pd.DataFrame)
    assert frame["age"].tolist() == [1, 2]
    assert frame[["age", "mass"]].columns.tolist() == ["age", "mass"]


def test_a_column_that_is_read_is_recorded(frame):
    frame["feh"]
    assert frame.record_of_all_used_keys == {"feh"}


def test_the_columns_that_are_read_together_are_recorded(frame):
    frame[["age", "mass"]]
    assert frame.record_of_all_used_keys == {"age", "mass"}


def test_reading_a_column_twice_records_it_once_and_the_record_is_a_set(frame):
    frame["age"]
    frame["age"]
    assert frame.record_of_all_used_keys == {"age"}
    assert isinstance(frame.record_of_all_used_keys, set)


def test_looking_at_the_frame_does_not_record_anything(frame):
    len(frame), frame.columns, frame.shape
    assert frame.record_of_all_used_keys == set()


def test_two_frames_record_apart():
    first, second = RecordingDataFrame({"age": [1]}), RecordingDataFrame({"age": [1]})
    first["age"]
    assert first.record_of_all_used_keys == {"age"} and second.record_of_all_used_keys == set()


def test_a_function_of_the_user_on_an_empty_frame_shows_which_columns_it_reads():
    # How Output uses it: the function is run on an empty frame that has all the columns, to learn which ones it needs
    frame = RecordingDataFrame([], columns=["px", "py", "pz"], dtype=float)
    def function(df):
        df["radius"] = df["px"] + df["py"]
    function(frame)
    assert frame.record_of_all_used_keys == {"px", "py"}


if __name__ == '__main__':
    pass
