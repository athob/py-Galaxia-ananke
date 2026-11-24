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


__all__ = ['Singleton', 'classproperty', 'State', 'execute', 'format_docstring', 'make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'common_entries', 'tuple_itemgetter', 'get_version_of_command', 'lexicalorder_dict', 'in_description_order', 'import_source_file', 'mark_metadata_prop', 'collect_metadata_marked_properties', 'metadata_dicts_to_consolidated_json']


class Singleton(type):
    """
    Metaclass whose classes have one instance, which is made on the first call and given again on the next ones.

    Taken from https://stackoverflow.com/questions/6760685/creating-a-singleton-in-python
    """
    _instances: Dict[type, Any] = {}
    def __call__(cls, *args: Any, **kwargs: Any) -> Any:
        """
        Give the one instance of the class, making it on the first call.

        Parameters
        ----------
        *args, **kwargs
            Given to the class when the instance is made, on the first call; they are ignored afterwards.

        Returns
        -------
        instance : object
            The one instance of the class.
        """
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
    Dataclass that tracks a state which moves forward: a register of counters, one for each level of work, and a
    description.

    States are compared by their register first, and then by their description.

    Attributes
    ----------
    register : list of int
        The counters, one for each level: the work of a level counts the work of the levels above it.
    description : str
        What the last update said that it was done.
    """
    register: List[int] = dc.field(default_factory=lambda: [0])
    description: str = "initial"
    @classmethod
    def fromdict(cls: Type[_BS], dictionary: Dict[str, Any]) -> _BS:
        """
        Make a state from the values of its fields.

        Parameters
        ----------
        dictionary : dict
            The values of ``register`` and ``description``, as ``asdict`` gives them.

        Returns
        -------
        state : _BaseState
            The state.
        """
        return cls(**dictionary)
    def asdict(self) -> Dict[str, Any]:
        """
        Give the values of the fields of the state.

        Returns
        -------
        dictionary : dict
            The values of ``register`` and ``description``.
        """
        return dc.asdict(self)
    def __increment_register(self, level: int) -> None:
        """
        Count one more for a level, and forget the deeper ones.

        The register grows with zeros to hold the level if it is too short, the counters of the levels below it are deleted,
        and the counter of the level is increased by one.

        Parameters
        ----------
        level : int
            The level that moves forward.
        """
        self.register += max(level+1 - len(self.register), 0)*[0]
        del self.register[level+1:]
        self.register[level] += 1
    def update(self, description: str = "", level: int = 0) -> None:
        """
        Move the state forward at a level.

        Parameters
        ----------
        description : str, optional
            What is being done. Default to an empty string.
        level : int, optional
            The level of the work that is done: its counter is increased, and those of the deeper levels start again.
            Default to 0.
        """
        self.__increment_register(level)
        self.description = description
    def __eq__(self, other: Any) -> bool:
        """
        Tell whether another state has the same register and the same description.

        Parameters
        ----------
        other : object
            What to compare to, which has to have a ``register`` and a ``description``.

        Returns
        -------
        equal : bool
            Whether both are the same.
        """
        return (self.register, self.description) == (other.register, other.description)
    def __lt__(self, other: "_BaseState") -> bool:
        """
        Tell whether the state comes before another, by their registers and then by their descriptions.

        Parameters
        ----------
        other : _BaseState
            The state to compare to.

        Returns
        -------
        before : bool
            Whether the register of this state, and if both are equal its description, sorts before those of the other.
        """
        return (self.register, self.description) < (other.register, other.description)


Param = ParamSpec("Param")


class State(_BaseState):
    """
    State that tracks an evolving process with a file: the state is written in a file, and a process that is run again can
    tell what it did before.
    """
    def __init__(self, filebase: Union[str, pathlib.Path], parent: Any, *args: Any, **kwargs: Any) -> None:
        """
        Make the state, and write its file if there is none yet.

        Parameters
        ----------
        filebase : str or path
            The path that the file of the state is named after: the file has the suffix ``.state``.
        parent : object
            The object that the state belongs to, whose ``caching`` and ``verbose`` it reads (an Output).
        *args, **kwargs
            Given to the base state: its register and its description.
        """
        super().__init__(*args, **kwargs) 
        self.filebase = filebase
        self.parent = parent
        if not self.file_path.exists():
            self.__writefile()
    @property
    def caching(self) -> bool:
        """Whether the parent caches its results, which makes the state skip the work that the file says is done."""
        return self.parent.caching
    @property
    def verbose(self) -> bool:
        """Whether the parent prints what it does."""
        return self.parent.verbose
    @property
    def file_path(self) -> pathlib.Path:
        """The path of the file of the state: its base with the suffix ``.state``, made absolute."""
        return pathlib.Path(self.filebase).with_suffix(".state").resolve()
    def __writefile(self) -> None:
        """Write the state in its file, as JSON."""
        self.file_path.write_bytes(bytes(json.dumps(self.asdict()),'ascii'))
    def __readfile(self) -> _BaseState:
        """
        Read the state in the file.

        Returns
        -------
        state : _BaseState
            The state that the file holds, or a new one that has not moved if there is no file.
        """
        file_path = self.file_path
        if file_path.exists():
            return _BaseState.fromdict(json.loads(self.file_path.read_bytes()))
        else:
            return _BaseState()
    @property
    def __is_behind_state_of_file(self) -> bool:
        """Whether the state is not ahead of the one that the file holds, which means that the work is already done."""
        return self <= self.__readfile()
    def check_state_file_before_running(self, *args: Any, **kwargs: Any
                                        ) -> Callable[[Callable[Param, None]], Callable[Param, None]]:
        """
        Make a decorator that runs a function only if its work is not done yet.

        When the decorated function is called, the state is updated with the arguments given here. If the parent caches, the
        function runs only if the state is then ahead of the state in the file, and without caching it always runs. After it
        ran the state is written in its file.

        Parameters
        ----------
        *args, **kwargs
            Given to ``update``: the description and the level of the work.

        Returns
        -------
        decorator : callable
            The decorator, which gives the function that checks the state before running its own.
        """
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
    Descriptor that makes a property of a class, read from the class itself and not from an instance.

    Credit https://stackoverflow.com/a/5192374
    """
    def __init__(self, f: Callable[[Any], Any]) -> None:
        """
        Hold the function that gives the property.

        Parameters
        ----------
        f : callable
            The function, which is given the class.
        """
        self.f = f
    def __get__(self, _: Any, owner: Any) -> Any:
        """
        Give the value of the property for a class.

        Parameters
        ----------
        _ : object
            The instance that it is read from, which is not used.
        owner : type
            The class, which is given to the function.

        Returns
        -------
        value : object
            What the function gives for the class.
        """
        return self.f(owner)


def _execute_generator(cmds: List[str], max_workers: Optional[int] = None, **kwargs: Any) -> Iterator[Tuple[int, str]]:
    """
    Run commands at the same time, and give their output line by line as they print it, errors included.

    The commands are read as a shell reads them (``shlex.split``), all of them before any is run, and are not run through a
    shell.

    Parameters
    ----------
    cmds : list of str
        The commands.
    max_workers : int, optional
        Not used: all the commands are started at once.
    **kwargs
        Given to subprocess.Popen (for instance ``cwd``).

    Yields
    ------
    tag : int
        The number of the command, counted from 1.
    line : str
        A line that the command printed, its end of line included (the last line may not have one).

    Raises
    ------
    subprocess.CalledProcessError
        After the last line, if a command exited with an error: the one of the first command that did.
    """
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
    Run commands at the same time, and print their output with the number of the command, errors included.

    Each command is a string that is read as a shell reads it (``shlex.split``): an argument that holds a space, a quote or
    anything that a shell would take apart has to be quoted, with ``shlex.quote`` or ``shlex.join``. The commands are not
    run through a shell.

    Adapted from https://stackoverflow.com/a/4417735

    Parameters
    ----------
    cmds : str or list of str
        The command, or the commands.
    max_workers : int, optional
        Not used: all the commands are started at once.
    verbose : bool, optional
        Whether to print the commands and their output, errors included, as it comes, each line after ``JOB n/N |``. Default
        to True.
    **kwargs
        Given to subprocess.Popen (for instance ``cwd``).

    Raises
    ------
    subprocess.CalledProcessError
        Once the output is read, if a command exited with an error.
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


def format_docstring(function: Callable[..., Any], **values: Any) -> None:
    """
    Fill the placeholders of the docstring of a function with values.

    This is how the documentation of the package gets its defaults. A method is given as itself: its function holds the
    docstring. When Python runs with -OO there are no docstrings, and nothing is done.

    Parameters
    ----------
    function : callable
        The function, or the method, whose docstring is filled in place.
    **values
        The values of the placeholders, by name, as ``str.format`` takes them.
    """
    holder: Any = getattr(function, '__func__', function)
    if holder.__doc__ is not None:
        holder.__doc__ = holder.__doc__.format(**values)


def make_symlink(file_path: Union[str, pathlib.Path], dest_dir: Union[str, pathlib.Path]) -> None:
    """
    Make a symbolic link to a file, in a directory, with the name of the file.

    Both paths are resolved first, and a link of that name that is already there is replaced.

    Parameters
    ----------
    file_path : str or path
        The file that the link points to.
    dest_dir : str or path
        The directory where the link is made.

    Raises
    ------
    IsADirectoryError
        If a directory, and not a file or a link, already has the name that the link would have.
    """
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
    """
    Check that a collection of keys has the required ones, and the optional ones only besides.

    Parameters
    ----------
    given : iterable
        The keys that were given.
    required : iterable, optional
        The keys that have to be there. Default to none.
    optional : iterable, optional
        The keys that may be there too. Default to none.
    error_message : str, optional
        The beginning of the message of the error. Default to a message about the keys of the particle data.

    Raises
    ------
    ValueError
        If the given keys, the optional ones apart, are not the required ones: the message names the ones that are missing
        and the ones that are not expected.
    """
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
    """
    Check that all the arrays of a dictionary have the same length.

    Parameters
    ----------
    dictionary : dict
        The arrays, or anything that has a length, by name.
    control : str, optional
        The name of the array that the others are compared with. Default to the first one of the dictionary.
    error_message_dict_name : str, optional
        What the message of the error calls the dictionary, before ``input dictionary``. Default to nothing.

    Raises
    ------
    ValueError
        If an array has not the length of the control: the message names them.
    """
    if control is None and dictionary:
        control = list(dictionary.keys())[0]
    wrong_keys = []
    for key in set(dictionary.keys()) - {control}:
        if len(dictionary[key]) != len(dictionary[cast(str, control)]): wrong_keys.append(key)
    if wrong_keys:
        raise ValueError(f"Array{'' if len(wrong_keys)==1 else 's'} representing propert{'y' if len(wrong_keys)==1 else 'ies'} {set(wrong_keys)} in the provided {error_message_dict_name + bool(error_message_dict_name)*' '}input dictionary do{'es' if len(wrong_keys)==1 else ''} not have the same length as property {control}.")


def tuple_itemgetter(*items: Any) -> Callable[[Any], Tuple[Any, ...]]:
    """
    Make a function that gets items, always as a tuple.

    ``operator.itemgetter`` gives the item itself, and not a tuple of one, when it is given a single one: this one does not.

    Parameters
    ----------
    *items
        The items to get, as ``operator.itemgetter`` takes them.

    Returns
    -------
    getter : callable
        The function that gives, for an object, the tuple of its items.
    """
    getter = itemgetter(*items)
    return getter if len(items) != 1 else (lambda obj: (getter(obj),))


def common_entries(*dcts: Dict[Any, Any]) -> Iterator[Tuple[Any, ...]]:
    """
    Give the entries of dictionaries that have the same key, as ``zip`` does for lists.

    Taken from https://stackoverflow.com/a/16458780

    Parameters
    ----------
    *dcts : dict
        The dictionaries.

    Yields
    ------
    entry : tuple
        A key that all the dictionaries have, followed by its value in each of them, in the order of the dictionaries. The
        keys come in no particular order. Nothing comes if there is no dictionary.
    """
    if not dcts:
        return
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)


def get_version_of_command(cmd: str) -> str:
    """
    Give the version that a command reports with ``--version``.

    Parameters
    ----------
    cmd : str
        The command, which has to be on the PATH.

    Returns
    -------
    version : str
        The first of the dotted numbers in what the command prints, such as ``"4.3"`` for ``GNU Make 4.3``.

    Raises
    ------
    FileNotFoundError
        If the command does not exist.
    subprocess.CalledProcessError
        If the command exited with an error.
    IndexError
        If what the command prints has no dotted number.
    """
    return re.findall(r"((?:[0-9]+\.)+[0-9]+)",
                      str(subprocess.check_output([cmd, '--version'])))[0]


def lexicalorder_dict(dictionary: Dict[Any, Any]) -> OrderedDict[Any, Any]:
    """
    Give a dictionary with the same items, in the order of the sorted keys.

    Parameters
    ----------
    dictionary : dict
        The dictionary, whose keys can be sorted.

    Returns
    -------
    ordered : collections.OrderedDict
        The same items, sorted by key.
    """
    return ODict({k: dictionary[k] for k in sorted(dictionary)})


def in_description_order(properties: Iterable[Tuple[Any, str]]) -> List[Tuple[Any, str]]:
    """
    Put properties in the alphabetical order of their descriptions.

    A set has no order that stays the same from one run to the next (it depends on the hash seed), so a docstring that lists
    the properties of a set does it in this order.

    Parameters
    ----------
    properties : iterable of tuple
        The (key, description) pairs of the properties. The keys, which may be strings or tuples, only break a tie.

    Returns
    -------
    ordered : list of tuple
        The pairs, sorted by description.
    """
    return sorted(properties, key=lambda item: (item[1], str(item[0])))


def import_source_file(module_name: str, file_path: Union[str, pathlib.Path]) -> ModuleType:
    # based on https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
    """
    Import a Python source file as a module of the given name.

    The module is registered in ``sys.modules``, and the file is read on every call, whatever was imported under that name
    before. Based on https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly

    Parameters
    ----------
    module_name : str
        The name that the module is given.
    file_path : str or path
        The source file.

    Returns
    -------
    module : types.ModuleType
        The module.

    Raises
    ------
    AttributeError
        If the file is not one that Python imports (its suffix is not known), which leaves nothing to load it with.
    FileNotFoundError
        If the file does not exist.
    """
    spec = cast(ModuleSpec, importlib.util.spec_from_file_location(module_name, file_path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    cast(Loader, spec.loader).exec_module(module)
    return module


METADATA_MARK = '_marked_for_metadata'


_PG = TypeVar('_PG', bound=Union[Callable[..., Any], property])   # a getter or a property


def mark_metadata_prop(prop_or_getter: _PG) -> _PG:
    """
    Mark a property getter, so that its property goes in the metadata of a class decorated with
    collect_metadata_marked_properties.

    It can be applied before or after ``@property``: after it, the getter of the property is the one that is marked.

    Parameters
    ----------
    prop_or_getter : callable or property
        The getter to mark, or the property whose getter is.

    Returns
    -------
    prop_or_getter : callable or property
        What was given, which is marked and otherwise unchanged.
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
    Class decorator that gives a class the metadata of its marked properties.

    It adds a ``_metadata`` property, which gives a dictionary of the names of the properties marked with
    ``mark_metadata_prop`` and their current values, and a ``_metadata_names`` property that holds the sorted names.

    Parameters
    ----------
    cls : type
        The class to decorate.

    Returns
    -------
    cls : type
        The same class, with the two properties.
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
    """
    Give the compact form of the metadata of several sources as JSON.

    Parameters
    ----------
    *args, **kwargs
        Given to ``create_compact_metadata``.

    Returns
    -------
    json_string : str
        The compact metadata, as JSON indented by 2 spaces.
    """
    metadata = create_compact_metadata(*args, **kwargs)
    json_string = json.dumps(metadata, indent=2)
    return json_string


class _CompactMetadata(TypedDict, total=False):
    """
    The compact form of the metadata of several sources, which create_compact_metadata gives and reconstruct_log reads.

    It is an empty dictionary when there is no metadata to compact.

    Attributes
    ----------
    metadata_name : str
        The name that the metadata was given.
    schema_version : str
        The version of this form.
    common_values : dict
        The values that every source has, and that are the same.
    varying_keys : list of str
        The other keys, in the order that most of the sources have them.
    metadata_count : int
        The number of sources.
    variations : list of dict
        For each source, its values of the varying keys.
    """
    metadata_name: str
    schema_version: str
    common_values: Dict[str, Union[str,int,float]]
    varying_keys: List[str]
    metadata_count: int
    variations: List[Dict[str, Union[str,int,float]]]


def create_compact_metadata(metadata_dicts: List[Dict[str, Union[str,int,float]]],
                            metadata_name: str = "") -> _CompactMetadata:
    """
    Give the metadata of several sources in a compact form that avoids what is redundant.

    Parameters
    ----------
    metadata_dicts : list of dict
        The metadata of each source, by key.
    metadata_name : str, optional
        A name for the metadata. Default to an empty string.

    Returns
    -------
    compact : dict
        The values that all the sources have in common, the other keys, and for each source the values of those: see
        ``_CompactMetadata``. An empty dictionary if no metadata is given.
    """
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
    """
    Give the keys of several dictionaries in the order that most of them have.

    Each key goes at the position where it is most often, and the keys that share a position come in the order that they
    first appear in.

    Parameters
    ----------
    metadata_dicts : list of dict
        The dictionaries.

    Returns
    -------
    keys : list of str
        All the keys that appear in any of the dictionaries. Empty if there is no dictionary.
    """
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
    """
    Give back the metadata of one source from the compact form.

    Parameters
    ----------
    compact_metadata : dict
        The compact metadata, as ``create_compact_metadata`` gives it.
    metadata_index : int
        The number of the source.

    Returns
    -------
    metadata : dict
        The common values, with the values of the source that vary.

    Raises
    ------
    IndexError
        If the number is not smaller than the number of sources.
    """
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
