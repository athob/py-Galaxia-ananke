#!/usr/bin/env python
"""
Module miscellaneous utilities
"""
from typing import Any, Protocol

import pandas as pd

from ._builtin_utils import *


__all__ = ['CallableDFtoNone', 'Singleton', 'State', 'execute', 'make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'common_entries']


class CallableDFtoNone(Protocol):
    def __call__(self, df: pd.DataFrame, *args: Any) -> None:  # TODO change DataFrame typing annotation to a "DataFrameLike" type if such exists (similar to ArrayLike)
        pass


if __name__ == '__main__':
    pass
