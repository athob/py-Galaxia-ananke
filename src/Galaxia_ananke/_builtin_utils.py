#!/usr/bin/env python
"""
Module utilities using built-in implementation
"""
from typing import Type, TypeVar, Any, Union, List, Dict, Callable
from typing_extensions import Self, ParamSpec
from functools import total_ordering
from itertools import zip_longest
import contextlib as cl
import dataclasses as dc
import subprocess
import pathlib
import json
import re


__all__ = ['Singleton', 'classproperty', 'State', 'execute', 'make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'common_entries', 'get_version_of_command']


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


if __name__ == '__main__':
    pass
