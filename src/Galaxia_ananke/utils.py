#!/usr/bin/env python
"""
Module miscellaneous utilities
"""
from typing import Union, List, Dict, Protocol
from itertools import zip_longest
import subprocess
import pathlib

import pandas as pd


__all__ = ['CallableDFtoNone', 'make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'Singleton', 'execute', 'common_entries']


class CallableDFtoNone(Protocol):
    def __call__(self, df: pd.DataFrame) -> None:  # TODO change DataFrame typing annotation to a "DataFrameLike" type if such exists (similar to ArrayLike)
        pass


def make_symlink(file_path, dest_dir):
    file_path = pathlib.Path(file_path).resolve()
    dest_dir = pathlib.Path(dest_dir).resolve()
    symlink_name = dest_dir / file_path.name
    try:
        symlink_name.unlink()
    except FileNotFoundError:
        pass
    symlink_name.symlink_to(file_path)


def compare_given_and_required(given, required=set(), optional=set(), error_message="Given particle data covers wrong set of keys"):
    given = set(given)
    required = set(required)
    optional = set(optional)
    if given-optional != required:
        missing = required.difference(given)
        missing = f"misses {missing}" if missing else ""
        extra = given.difference(required.union(optional))
        extra = f"misincludes {extra}" if extra else ""
        raise ValueError(f"{error_message}: {missing}{' & ' if missing and extra else ''}{extra}")


def confirm_equal_length_arrays_in_dict(dictionary: dict, control: str = None, error_message_dict_name: str = ''):
    if control is None and dictionary:
        control = list(dictionary.keys())[0]
    wrong_keys = []
    for key in set(dictionary.keys()) - {control}:
        if len(dictionary[key]) != len(dictionary[control]): wrong_keys.append(key)
    if wrong_keys:
        raise ValueError(f"Array{'' if len(wrong_keys)==1 else 's'} representing propert{'y' if len(wrong_keys)==1 else 'ies'} {set(wrong_keys)} in the provided {error_message_dict_name + bool(error_message_dict_name)*' '}input dictionary do{'es' if len(wrong_keys)==1 else ''} not have the same length as property {control}.")


class Singleton(type):
    """
    Singleton metaclass. Directly taken from
    https://stackoverflow.com/questions/6760685/creating-a-singleton-in-python
    """
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


# class CachedInstance(type):
#     _instances = {}
#     def __call__(cls, *args, **kwargs):
#         index = cls, args
#         if index not in cls._instances:
#             cls._instances[index] = super(CachedInstance, cls).__call__(*args, **kwargs)
#         return cls._instances[index]


def _execute_generator(cmds: List[str], **kwargs):
    popens = [subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True, **kwargs)
              for cmd in cmds]
    master_stdout_readline_iter = zip_longest(*(iter(popen.stdout.readline, "") for popen in popens))
    for stdout_lines in master_stdout_readline_iter:
        for tag, stdout_line in enumerate(stdout_lines, start=1):
            yield tag, stdout_line
    return_codes = [popen.wait() for popen in popens if popen.stdout.close() is None]
    for return_code, cmd in zip(return_codes, cmds):
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)


def execute(cmds: Union[str, List[str]], verbose: bool = True, **kwargs):
    """
    Run the commands described by cmds, and use
    verbose kwarg to redirect output/error stream
    to python output stream.
    Adapted from https://stackoverflow.com/a/4417735
    """
    if isinstance(cmds, str):
        cmds: List[str] = [cmds]
    n_cmds = len(cmds)
    len_tags = len(str(n_cmds))+1
    if verbose:
        for proc_id, cmd in enumerate(cmds, start=1):
            print(f"Executing JOB{proc_id: {len_tags}d}/{n_cmds} = {cmd}")
    for proc_id, stdout_line in _execute_generator(cmds, **kwargs):
        print(f"JOB{proc_id: {len_tags}d}/{n_cmds} | {stdout_line}", end="") if verbose else None


def common_entries(*dcts: Dict):
    """
    common_entries function equivalent to zip in dictionaries. Directly taken from
    https://stackoverflow.com/a/16458780
    """
    if not dcts:
        return
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)


if __name__ == '__main__':
    pass
