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
Tests of how an Output with a catalogue is used like the vaex table that it holds: columns by index and by attribute, the methods
of the table, how it is shown, completed and listed. The catalogue is made of the HDF5 files of synthetic output of the backend.
"""
import pytest

from .synthetic import read_output

EXPORT_COUNT = 22            # the keys that the backend gives and the three magnitudes of Gaia


@pytest.fixture
def made(tmp_path, monkeypatch):
    output, written = read_output(tmp_path, monkeypatch, n_stars=6, partitions=[0, 0, 0, 1, 1, 1], pp_auto_flush=False)
    return output, written[0][1]


def test_a_column_is_reached_by_its_name_and_by_attribute(made):
    output, columns = made
    assert output['px'].to_numpy().tolist() == columns['px'].tolist()
    assert output.px.to_numpy().tolist() == columns['px'].tolist()
    assert output['minit'].to_numpy().tolist() == columns['minit'].tolist()


def test_the_methods_of_the_table_are_the_methods_of_the_output(made):
    output, columns = made
    assert output.shape == (6, EXPORT_COUNT) and output.count() == 6
    assert set(output.get_column_names()) == set(output.export_keys)
    assert len(output.to_pandas_df()) == 6 and output.to_pandas_df()['px'].tolist() == columns['px'].tolist()


def test_an_output_is_shown_as_its_table(made):
    output, _ = made
    assert repr(output) == repr(output._vaex)


def test_a_column_that_is_set_goes_to_the_table(made):
    output, columns = made
    output['double_px'] = output.px * 2
    assert output._vaex['double_px'].to_numpy().tolist() == (columns['px'] * 2).tolist()
    assert 'double_px' in output.get_column_names()


def test_what_is_listed_is_what_the_table_has_and_what_the_output_has(made):
    output, _ = made
    listed = output.__dir__()                  # dir() sorts what it is given: only the method itself can tell whether it sorts
    assert {'px', 'minit', 'to_pandas_df', 'get_column_names'} <= set(listed)           # the table
    assert {'export_keys', 'output_dir', 'read_galaxia_output', 'post_process_output'} <= set(listed)     # the output
    assert listed == sorted(listed) and set(dir(output)) == set(listed)


def test_the_table_is_shown_as_html_in_notebooks(made):
    output, _ = made
    assert isinstance(output._repr_html_(), str) and 'px' in output._repr_html_() and '_repr_html_' in dir(output)


def test_an_attribute_that_nothing_has_is_an_attribute_error(made):
    output, _ = made
    with pytest.raises(AttributeError):
        output.nothing_of_the_sort


def test_the_private_attributes_of_the_table_are_not_those_of_the_output(made):
    output, _ = made
    with pytest.raises(AttributeError):
        output._get_column_names_that_do_not_exist


def test_the_length_of_an_output_is_the_number_of_stars(made):
    output, _ = made
    assert len(output) == 6 == output.shape[0] == len(output._vaex)
    assert bool(output) is True                          # as for any container, what has something in it is true


def test_the_length_and_the_truth_of_an_output_need_its_catalogue_to_be_loaded(tmp_path, monkeypatch):
    output, _ = read_output(tmp_path, monkeypatch, n_stars=6, partitions=[0, 0, 0, 1, 1, 1], pp_auto_flush=True)       # nothing loaded yet
    for ask in (len, bool):
        with pytest.raises(RuntimeError, match="Don't attempt creating an Output object on your own"):
            ask(output)


if __name__ == '__main__':
    pass
