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
"""
Module miscellaneous utilities
"""
from typing import Any, Protocol, Iterable, Set, Union
import hashlib

import pandas as pd
from numpy.typing import ArrayLike

from ._builtin_utils import *
from ._constants import HASH_ENCODING
from ._patch_vaex import apply_vaex_patch


__all__ = ['classproperty', 'CallableDFtoNone', 'CallableDFtoInt', 'RecordingDataFrame', 'Singleton', 'State', 'execute', 'format_docstring', 'make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'common_entries', 'tuple_itemgetter', 'lexicalorder_dict', 'in_description_order', 'mark_metadata_prop', 'collect_metadata_marked_properties', 'metadata_dicts_to_consolidated_json', 'hash_iterable']


apply_vaex_patch()

class CallableDFtoNone(Protocol):
    """
    The type of a function that takes the table of a partition of the catalogue, and then other positional arguments, and
    gives nothing: the post processing pipelines that ``Output.apply_post_process_pipeline_and_flush`` applies, which change
    the table in place.
    """
    def __call__(self, df: pd.DataFrame, *args: Any) -> None:  # TODO change DataFrame typing annotation to a "DataFrameLike" type if such exists (similar to ArrayLike)
        """
        Process the table of a partition of the catalogue in place.

        Parameters
        ----------
        df : vaex.dataframe.DataFrame
            The table of a partition: pipelines are given a vaex DataFrame, the one opened on the HDF5 file of the partition
            or the one held in memory. The type hint says pandas.DataFrame for now (a note in the code says that it is to
            change): what the pipelines of the package do with the table works the same with a pandas one, which is what the
            tests give them.
        *args
            Other positional arguments, as the user gave them.
        """
        pass


class CallableDFtoInt(Protocol):
    """
    The type of a function that takes the table of the catalogue and gives the partition of each of its stars: the
    partitioning rules that ``Survey.make_survey`` takes.

    The rule is first applied to an empty table, to know the columns that it uses, and then to a table made of those
    columns.
    """
    def __call__(self, df: pd.DataFrame, *args: Any) -> ArrayLike:  # TODO change DataFrame typing annotation to a "DataFrameLike" type if such exists (similar to ArrayLike)
        """
        Give the partition of each star of the table of the catalogue.

        Parameters
        ----------
        df : pandas.DataFrame
            The table.
        *args
            Other positional arguments.

        Returns
        -------
        partition : array_like of int
            The number of the partition of each star, one for each row of the table. A star with a negative number is in no
            partition, and is left out of the catalogue.
        """
        pass


class RecordingDataFrame(pd.DataFrame):
    """
    Pandas DataFrame that records all the keys that are used to get items from it.

    The record is what a function of the table needs from it: it is read, for a function applied to an empty table, to know
    the columns that the function uses.
    """
    _record_of_all_used_keys: Set[Any]
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """
        Make the table and an empty record.

        Parameters
        ----------
        *args, **kwargs
            Given to ``pandas.DataFrame``.
        """
        super().__init__(*args, **kwargs)
        object.__setattr__(self, '_record_of_all_used_keys', set())      # not through pandas, which takes a new attribute for a column that someone forgot to make, and warns
    def _add_to_record_of_all_used_keys(self, keys: Union[str, Iterable[Any]]) -> None:
        """
        Add keys to the record.

        Parameters
        ----------
        keys : str or iterable
            A key, if it is a string, or the keys.
        """
        if isinstance(keys, str):
            keys = [keys]
        for key in keys:
            self._record_of_all_used_keys.add(key)
    def __getitem__(self, key: Any) -> Any:
        """
        Record the key, or the keys, and get the item as pandas does.

        Parameters
        ----------
        key : object
            What the item is got with, as ``pandas.DataFrame.__getitem__`` takes it.

        Returns
        -------
        item : object
            What pandas gives.
        """
        self._add_to_record_of_all_used_keys(key)
        return super().__getitem__(key)
    # def __setitem__(self, key, value):
    #     self._add_to_record_of_all_used_keys(key)
    #     super().__setitem__(key, value)
    # def __delitem__(self, key):
    #     self._add_to_record_of_all_used_keys(key)
    #     super().__delitem__(key)
    @property
    def record_of_all_used_keys(self) -> Set[Any]:
        """The set of the keys that were used to get items of the table."""
        return self._record_of_all_used_keys


def hash_iterable(iterable: Iterable) -> bytes:
    """
    Give a hash of an iterable of arrays.

    The hash is part of the names of the files of the cache of an input: a change in how it is computed leaves the caches of
    the earlier runs unused.

    Parameters
    ----------
    iterable : iterable
        The arrays, or anything else that ``hashlib`` hashes as bytes.

    Returns
    -------
    digest : bytes
        The hexadecimal digest, encoded, of the hexadecimal digests of the elements.
    """
    return bytes(hashlib.sha256(
        bytes('\n'.join([hashlib.sha256(element).hexdigest()
                         for element in iterable]),
                HASH_ENCODING)
        ).hexdigest(), HASH_ENCODING)


if __name__ == '__main__':
    raise NotImplementedError()
