#!/usr/bin/env python
"""
Module miscellaneous utilities
"""
from typing import Any, Protocol, Iterable
import hashlib

import pandas as pd

from ._builtin_utils import *
from ._constants import HASH_ENCODING


__all__ = ['classproperty', 'CallableDFtoNone', 'CallableDFtoInt', 'RecordingDataFrame', 'Singleton', 'State', 'execute', 'make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'common_entries', 'lexicalorder_dict', 'hash_iterable']


class CallableDFtoNone(Protocol):
    def __call__(self, df: pd.DataFrame, *args: Any) -> None:  # TODO change DataFrame typing annotation to a "DataFrameLike" type if such exists (similar to ArrayLike)
        pass


class CallableDFtoInt(Protocol):
    def __call__(self, df: pd.DataFrame, *args: Any) -> int:  # TODO change DataFrame typing annotation to a "DataFrameLike" type if such exists (similar to ArrayLike)
        pass


class RecordingDataFrame(pd.DataFrame):
    """
    Pandas DataFrame that records all its used keys from getitem
    """
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._record_of_all_used_keys = set()
    def _add_to_record_of_all_used_keys(self, keys):
        if isinstance(keys, str):
            keys = [keys]
        for key in keys:
            self._record_of_all_used_keys.add(key)
    def __getitem__(self, key):
        self._add_to_record_of_all_used_keys(key)
        return super().__getitem__(key)
    # def __setitem__(self, key, value):
    #     self._add_to_record_of_all_used_keys(key)
    #     super().__setitem__(key, value)
    # def __delitem__(self, key):
    #     self._add_to_record_of_all_used_keys(key)
    #     super().__delitem__(key)
    @property
    def record_of_all_used_keys(self):
        return self._record_of_all_used_keys


def hash_iterable(iterable: Iterable) -> bytes:
    return bytes(hashlib.sha256(
        bytes('\n'.join([hashlib.sha256(element).hexdigest()
                         for element in iterable]),
                HASH_ENCODING)
        ).hexdigest(), HASH_ENCODING)


if __name__ == '__main__':
    pass
