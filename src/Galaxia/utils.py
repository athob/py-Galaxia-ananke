#!/usr/bin/env python
"""
Docstring
"""
__all__ = ['Singleton']


class Singleton(type):
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
