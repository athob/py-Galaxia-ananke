#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2022  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 (AGPL-3.0).
# 
# The full copyright notice, including terms governing use, modification,
# and redistribution, is contained in the files LICENSE and COPYRIGHT,
# which can be found at the root of the source code distribution tree:
# - LICENSE <https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE>
# - COPYRIGHT <https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT>
#
"""
Module utilities using built-in implementation
"""
import sys
from types import ModuleType
from typing import Type, TypeVar, Any, Union, List, Dict, OrderedDict, Callable
from typing_extensions import Self, ParamSpec
from collections import OrderedDict as ODict
from collections import Counter, defaultdict
from functools import cached_property, total_ordering
from itertools import zip_longest
import dataclasses as dc
import subprocess
import importlib.util
import pathlib
import json
import re


__all__ = ['Singleton', 'classproperty', 'State', 'execute', 'make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'common_entries', 'get_version_of_command', 'lexicalorder_dict', 'import_source_file', 'mark_metadata_prop', 'collect_metadata_marked_properties', 'metadata_dicts_to_consolidated_json']


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


_BS = TypeVar('_BS', bound='_BaseState')


@total_ordering
@dc.dataclass()
class _BaseState:
    """
    Base dataclass to track an evolving state.
    """
    register: List[int] = dc.field(default_factory=lambda: [0])
    description: str = "initial"
    @classmethod
    def fromdict(cls: Type[_BS], dictionary: Dict[str, Any]) -> _BS:
        return cls(**dictionary)
    def asdict(self) -> Dict[str, Any]:
        return dc.asdict(self)
    def __increment_register(self, level: int) -> None:
        self.register += max(level+1 - len(self.register), 0)*[0]
        del self.register[level+1:]
        self.register[level] += 1
    def update(self, description: str = "", level: int = 0) -> None:
        self.__increment_register(level)
        self.description = description
    def __eq__(self, other: Self) -> bool:
        return (self.register, self.description) == (other.register, other.description)
    def __lt__(self, other: Self) -> bool:
        return (self.register, self.description) < (other.register, other.description)


Param = ParamSpec("Param")


class State(_BaseState):
    """
    Class to track an evolving state with a file.
    """
    def __init__(self, filebase: Union[None, str, pathlib.Path], parent, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs) 
        self.filebase = filebase
        self.parent = parent
        if not self.file_path.exists():
            self.__writefile()
    @property
    def caching(self):
        return self.parent.caching
    @property
    def verbose(self):
        return self.parent.verbose
    @property
    def file_path(self) -> pathlib.Path:
        return pathlib.Path(self.filebase).with_suffix(".state").resolve()
    def __writefile(self) -> None:
        self.file_path.write_bytes(bytes(json.dumps(self.asdict()),'ascii'))
    def __readfile(self) -> _BaseState:
        file_path = self.file_path
        if file_path.exists():
            return _BaseState.fromdict(json.loads(self.file_path.read_bytes()))
        else:
            return _BaseState()
    @property
    def __is_behind_state_of_file(self) -> bool:
        return self <= self.__readfile()
    def check_state_file_before_running(self, *args, **kwargs):
        def decorator(func: Callable[Param, None]) -> Callable[Param, None]:
            def wrapper(*w_args, **w_kwargs) -> None:
                self.update(*args, **kwargs)
                if (not self.__is_behind_state_of_file
                    if self.caching else True):
                    func(*w_args, **w_kwargs)
                    self.__writefile()
            return wrapper
        return decorator


class classproperty(object):
    """
    Credit https://stackoverflow.com/a/5192374
    """
    def __init__(self, f):
        self.f = f
    def __get__(self, _, owner):
        return self.f(owner)


def _execute_generator(cmds: List[str], max_workers: int = None, **kwargs):
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


def execute(cmds: Union[str, List[str]], max_workers: int = None, verbose: bool = True, **kwargs):
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
    for proc_id, stdout_line in _execute_generator(cmds, max_workers=max_workers, **kwargs):
        print(f"JOB{proc_id: {len_tags}d}/{n_cmds} | {stdout_line}", end="") if verbose else None


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


def common_entries(*dcts: Dict):
    """
    common_entries function equivalent to zip in dictionaries. Directly taken from
    https://stackoverflow.com/a/16458780
    """
    if not dcts:
        return
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)


def get_version_of_command(cmd):
    return re.findall("((?:[0-9]+\.)+[0-9]+)",
                      str(subprocess.check_output([cmd, '--version'])))[0]


def lexicalorder_dict(dictionary: Dict[str, Any]) -> OrderedDict[str, Any]:
    return ODict({k: dictionary[k] for k in sorted(dictionary)})


def import_source_file(module_name, file_path) -> ModuleType:
    # based on https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


METADATA_MARK = '_marked_for_metadata'


def mark_metadata_prop(prop_or_getter: Union[Callable, property]) -> Union[Callable, property]:
    """
    Marks a property getter so that the corresponding property can be
    included in the parent Class collected dictionary of metadata
    properties. Can be applied before OR after @property without causing
    an AttributeError.
    """
    if isinstance(prop_or_getter, property):
        # Applied after @property – mark the getter instead.
        if prop_or_getter.fget is not None:
            mark_metadata_prop(prop_or_getter.fget)
    else:
        # Applied before @property – mark the function.
        setattr(prop_or_getter, METADATA_MARK, True)
    return prop_or_getter


C = TypeVar('C', bound=type)   # Generic for any class


def collect_metadata_marked_properties(cls: C) -> C:
    """
    Class decorator that adds a `_metadata` property.
    This property returns a dictionary whose keys are the names of
    all properties decorated with `@mark_for_metadata` and whose values
    are the current values of those properties on the instance.
    """
    # Define property that gather names of all marked properties.
    def _metadata_names(self: Self) -> List[str]:
        return [
            name for name, attr in self.__class__.__dict__.items()
            if ((getattr(attr.fget, METADATA_MARK)
                 if hasattr(attr.fget, METADATA_MARK) else False)
                if isinstance(attr, property) else False)
        ]
    # Add it to the class.
    cls._metadata_names = cached_property(_metadata_names)
    cls._metadata_names.__set_name__(cls, '_metadata_names')
    # Define property that builds the dictionary.
    def _metadata(self: Self) -> Dict[str, Any]:
        return {name: getattr(self, name) for name in self._metadata_names}
    # Add it to the class.
    cls._metadata = property(_metadata)
    return cls


def metadata_dicts_to_consolidated_json(*args, **kwargs):
    # Create compact metadata
    metadata = create_compact_metadata(*args, **kwargs)
    json_string = json.dumps(metadata, indent=2)
    return json_string


def create_compact_metadata(metadata_dicts: List[Dict[str, Union[str,int,float]]],
                            metadata_name: str = "") -> Dict[str, Union[str,int,List[str],List[Dict[str, Union[str,int,float]]]]]:
    """Create a compact JSON representation that avoids redundant data"""
    if not metadata_dicts:
        return {}
    # Find common values across all metadata dictionaries
    all_keys = get_majority_voting_ordered_common_keys(metadata_dicts)
    common_values = {}
    varying_keys = []
    for key in all_keys:
        values = [metadata_dict.get(key) for metadata_dict in metadata_dicts]
        # Check if all values are the same (and key exists in all metadata dicts)
        if all(v == values[0] for v in values) and None not in values:
            common_values[key] = values[0]
        else:
            varying_keys.append(key)
    # Create per-metadata variations (only store what changes from common)
    variations = []
    for metadata_dict in metadata_dicts:
        variation = {}
        for key in varying_keys:
            if key in metadata_dict and metadata_dict[key] != common_values.get(key):
                variation[key] = metadata_dict[key]
        variations.append(variation)
    # Build compact metadata structure
    metadata = {
        'metadata_name':  metadata_name,
        'schema_version': '1.0',
        'common_values':  common_values,
        'varying_keys':   varying_keys,
        'metadata_count': len(metadata_dicts),
        'variations':     variations
    }
    return metadata


def get_majority_voting_ordered_common_keys(metadata_dicts: List[Dict[str, Union[str,int,float]]]):
    """Use the most common ordering pattern across logs"""
    if not metadata_dicts:
        return []
    # Get all unique keys
    all_keys = set().union(*(metadata_dict.keys() for metadata_dict in metadata_dicts))
    # Build position frequency for each key
    position_freq = defaultdict(Counter)
    for metadata_dict in metadata_dicts:
        for pos, key in enumerate(metadata_dict.keys()):
            position_freq[key][pos] += 1
    # Assign each key to its most common position
    key_positions = {}
    for key, freq_dict in position_freq.items():
        most_common_pos = freq_dict.most_common(1)[0][0]
        key_positions[key] = most_common_pos
    # Sort by position, with ties broken by first appearance
    ordered_keys = sorted(
        all_keys,
        key=lambda k: (key_positions[k], min(i for i, metadata_dict in enumerate(metadata_dicts) if k in metadata_dict))
    )
    return ordered_keys


def reconstruct_log(compact_metadata: Dict[str, Union[str,int,List[str],List[Dict[str, Union[str,int,float]]]]],
                    metadata_index: int) -> Dict[str, Union[str,int,float]]:
    """Reconstruct original metadata dictionary from compact metadata"""
    if metadata_index >= compact_metadata['metadata_count']:
        raise IndexError(f"Metadata index {metadata_index} out of range")
    # Start with common values
    metadata_dict = compact_metadata['common_values'].copy()
    # Apply variations for this specific metadata
    variation = compact_metadata['variations'][metadata_index]
    metadata_dict.update(variation)
    return metadata_dict


if __name__ == '__main__':
    raise NotImplementedError()
