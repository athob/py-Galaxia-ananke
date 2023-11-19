#!/usr/bin/env python
"""
Module miscellaneous utilities
"""
import pathlib

__all__ = ['make_symlink', 'compare_given_and_required', 'confirm_equal_length_arrays_in_dict', 'Singleton']


def make_symlink(file_path, dest_dir):
    file_path = pathlib.Path(file_path).resolve()
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
    for key in set(dictionary.keys()) - {control}:
        assert len(dictionary[key]) == len(dictionary[control]), f"Array representing property {key} in the provided {error_message_dict_name + bool(error_message_dict_name)*' '}input dictionary does not have the same length as property {control}."

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


if __name__ == '__main__':
    pass
