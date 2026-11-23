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
Module utilities using built-in implementation
"""
import sys
from types import ModuleType
from typing import (Type, TypeVar, Any, Union, Optional, List, Dict, DefaultDict, OrderedDict, Callable, Iterator, Iterable,
                    Hashable, Mapping, Tuple, TypedDict, IO, cast)
if sys.version_info >= (3, 10):
    from typing import ParamSpec
else:
    from typing_extensions import ParamSpec
from collections import OrderedDict as ODict
from collections import Counter, defaultdict
from functools import cached_property, total_ordering
from itertools import zip_longest
from operator import itemgetter
import dataclasses as dc
import shlex
import subprocess
import importlib.util
from importlib.abc import Loader
from importlib.machinery import ModuleSpec
import pathlib
import json
import re


__all__ = ['Singleton', 'classproperty', 'State', 'execute', 'make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'common_entries', 'tuple_itemgetter', 'get_version_of_command', 'lexicalorder_dict', 'import_source_file', 'mark_metadata_prop', 'collect_metadata_marked_properties', 'metadata_dicts_to_consolidated_json']


class Singleton(type):
    """
    Singleton metaclass. Directly taken from
    https://stackoverflow.com/questions/6760685/creating-a-singleton-in-python
    """
    _instances: Dict[type, Any] = {}
    def __call__(cls, *args: Any, **kwargs: Any) -> Any:
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
    def __eq__(self, other: Any) -> bool:
        return (self.register, self.description) == (other.register, other.description)
    def __lt__(self, other: "_BaseState") -> bool:
        return (self.register, self.description) < (other.register, other.description)


Param = ParamSpec("Param")


class State(_BaseState):
    """
    Class to track an evolving state with a file.
    """
    def __init__(self, filebase: Union[str, pathlib.Path], parent: Any, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs) 
        self.filebase = filebase
        self.parent = parent
        if not self.file_path.exists():
            self.__writefile()
    @property
    def caching(self) -> bool:
        return self.parent.caching
    @property
    def verbose(self) -> bool:
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
    def check_state_file_before_running(self, *args: Any, **kwargs: Any
                                        ) -> Callable[[Callable[Param, None]], Callable[Param, None]]:
        def decorator(func: Callable[Param, None]) -> Callable[Param, None]:
            def wrapper(*w_args: Param.args, **w_kwargs: Param.kwargs) -> None:
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
    def __init__(self, f: Callable[[Any], Any]) -> None:
        self.f = f
    def __get__(self, _: Any, owner: Any) -> Any:
        return self.f(owner)


def _execute_generator(cmds: List[str], max_workers: Optional[int] = None, **kwargs: Any) -> Iterator[Tuple[int, str]]:
    arguments = [shlex.split(cmd) for cmd in cmds]          # all of them are read before any is run
    popens = [subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True, **kwargs)
              for argv in arguments]
    master_stdout_readline_iter = zip_longest(*(iter(cast(IO[str], popen.stdout).readline, "") for popen in popens))
    for stdout_lines in master_stdout_readline_iter:
        for tag, stdout_line in enumerate(stdout_lines, start=1):
            if stdout_line is not None:     # zip_longest gives None for a command that has no more lines
                yield tag, stdout_line
    return_codes = [popen.wait() for popen in popens if cast(IO[str], popen.stdout).close() is None]
    for return_code, cmd in zip(return_codes, cmds):
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)


def execute(cmds: Union[str, List[str]], max_workers: Optional[int] = None, verbose: bool = True, **kwargs: Any) -> None:
    """
    Run the commands described by cmds, and use
    verbose kwarg to redirect output/error stream
    to python output stream.
    Each command is a string that is read as a shell reads it (``shlex.split``): an argument that holds a space, a quote or
    anything that the shell would take apart has to be quoted, with ``shlex.quote`` or ``shlex.join``.
    Adapted from https://stackoverflow.com/a/4417735
    """
    commands: List[str] = [cmds] if isinstance(cmds, str) else cmds
    n_cmds = len(commands)
    len_tags = len(str(n_cmds))+1
    if verbose:
        for proc_id, cmd in enumerate(commands, start=1):
            print(f"Executing JOB{proc_id: {len_tags}d}/{n_cmds} = {cmd}")
    for proc_id, stdout_line in _execute_generator(commands, max_workers=max_workers, **kwargs):
        if verbose:
            # only the last line of a command can lack its end of line, and the next line printed would be glued to it
            print(f"JOB{proc_id: {len_tags}d}/{n_cmds} | {stdout_line}", end="" if stdout_line.endswith("\n") else "\n")


def make_symlink(file_path: Union[str, pathlib.Path], dest_dir: Union[str, pathlib.Path]) -> None:
    file_path = pathlib.Path(file_path).resolve()
    dest_dir = pathlib.Path(dest_dir).resolve()
    symlink_name = dest_dir / file_path.name
    try:
        symlink_name.unlink()
    except FileNotFoundError:
        pass
    symlink_name.symlink_to(file_path)


def compare_given_and_required(given: Iterable[Hashable], required: Iterable[Hashable] = set(),
                               optional: Iterable[Hashable] = set(),
                               error_message: str = "Given particle data covers wrong set of keys") -> None:
    given_keys = set(given)
    required_keys = set(required)
    optional_keys = set(optional)
    if given_keys-optional_keys != required_keys:
        missing_keys = required_keys.difference(given_keys)
        missing = f"misses {missing_keys}" if missing_keys else ""
        extra_keys = given_keys.difference(required_keys.union(optional_keys))
        extra = f"misincludes {extra_keys}" if extra_keys else ""
        raise ValueError(f"{error_message}: {missing}{' & ' if missing and extra else ''}{extra}")


def confirm_equal_length_arrays_in_dict(dictionary: Mapping[str, Any], control: Optional[str] = None,
                                        error_message_dict_name: str = '') -> None:
    if control is None and dictionary:
        control = list(dictionary.keys())[0]
    wrong_keys = []
    for key in set(dictionary.keys()) - {control}:
        if len(dictionary[key]) != len(dictionary[cast(str, control)]): wrong_keys.append(key)
    if wrong_keys:
        raise ValueError(f"Array{'' if len(wrong_keys)==1 else 's'} representing propert{'y' if len(wrong_keys)==1 else 'ies'} {set(wrong_keys)} in the provided {error_message_dict_name + bool(error_message_dict_name)*' '}input dictionary do{'es' if len(wrong_keys)==1 else ''} not have the same length as property {control}.")


def tuple_itemgetter(*items: Any) -> Callable[[Any], Tuple[Any, ...]]:
    """
    operator.itemgetter that always gives a tuple: itemgetter gives the item itself, and not a tuple of one, when it is given a single one
    """
    getter = itemgetter(*items)
    return getter if len(items) != 1 else (lambda obj: (getter(obj),))


def common_entries(*dcts: Dict[Any, Any]) -> Iterator[Tuple[Any, ...]]:
    """
    common_entries function equivalent to zip in dictionaries. Directly taken from
    https://stackoverflow.com/a/16458780
    """
    if not dcts:
        return
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)


def get_version_of_command(cmd: str) -> str:
    return re.findall(r"((?:[0-9]+\.)+[0-9]+)",
                      str(subprocess.check_output([cmd, '--version'])))[0]


def lexicalorder_dict(dictionary: Dict[str, Any]) -> OrderedDict[str, Any]:
    return ODict({k: dictionary[k] for k in sorted(dictionary)})


def import_source_file(module_name: str, file_path: Union[str, pathlib.Path]) -> ModuleType:
    # based on https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
    spec = cast(ModuleSpec, importlib.util.spec_from_file_location(module_name, file_path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    cast(Loader, spec.loader).exec_module(module)
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
    def _metadata_names(self: Any) -> List[str]:
        return sorted([
            name for name, attr in self.__class__.__dict__.items()
            if ((getattr(attr.fget, METADATA_MARK)
                 if hasattr(attr.fget, METADATA_MARK) else False)
                if isinstance(attr, property) else False)
        ])
    # Add it to the class (whose type does not declare the attributes that it gets).
    target: Any = cls
    target._metadata_names = cached_property(_metadata_names)
    target._metadata_names.__set_name__(cls, '_metadata_names')
    # Define property that builds the dictionary.
    def _metadata(self: Any) -> Dict[str, Any]:
        return {name: getattr(self, name) for name in self._metadata_names}
    # Add it to the class.
    target._metadata = property(_metadata)
    return cls


def metadata_dicts_to_consolidated_json(*args: Any, **kwargs: Any) -> str:
    # Create compact metadata
    metadata = create_compact_metadata(*args, **kwargs)
    json_string = json.dumps(metadata, indent=2)
    return json_string


class _CompactMetadata(TypedDict, total=False):
    # what create_compact_metadata gives: {} when it is given no metadata
    metadata_name: str
    schema_version: str
    common_values: Dict[str, Union[str,int,float]]
    varying_keys: List[str]
    metadata_count: int
    variations: List[Dict[str, Union[str,int,float]]]


def create_compact_metadata(metadata_dicts: List[Dict[str, Union[str,int,float]]],
                            metadata_name: str = "") -> _CompactMetadata:
    """Create a compact JSON representation that avoids redundant data"""
    if not metadata_dicts:
        return {}
    # Find common values across all metadata dictionaries
    all_keys = get_majority_voting_ordered_common_keys(metadata_dicts)
    common_values: Dict[str, Union[str,int,float]] = {}
    varying_keys: List[str] = []
    for key in all_keys:
        values = [metadata_dict.get(key) for metadata_dict in metadata_dicts]
        # Check if all values are the same (and key exists in all metadata dicts)
        if all(v == values[0] for v in values) and None not in values:
            common_values[key] = metadata_dicts[0][key]
        else:
            varying_keys.append(key)
    # Create per-metadata variations (only store what changes from common)
    variations: List[Dict[str, Union[str,int,float]]] = []
    for metadata_dict in metadata_dicts:
        variation: Dict[str, Union[str,int,float]] = {}
        for key in varying_keys:
            if key in metadata_dict and metadata_dict[key] != common_values.get(key):
                variation[key] = metadata_dict[key]
        variations.append(variation)
    # Build compact metadata structure
    metadata: _CompactMetadata = {
        'metadata_name':  metadata_name,
        'schema_version': '1.0',
        'common_values':  common_values,
        'varying_keys':   varying_keys,
        'metadata_count': len(metadata_dicts),
        'variations':     variations
    }
    return metadata


def get_majority_voting_ordered_common_keys(metadata_dicts: List[Dict[str, Union[str,int,float]]]) -> List[str]:
    """Use the most common ordering pattern across logs"""
    if not metadata_dicts:
        return []
    # Get all unique keys
    all_keys = set().union(*(metadata_dict.keys() for metadata_dict in metadata_dicts))
    # Build position frequency for each key
    position_freq: DefaultDict[str, Counter[int]] = defaultdict(Counter)
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


def reconstruct_log(compact_metadata: _CompactMetadata, metadata_index: int) -> Dict[str, Union[str,int,float]]:
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
