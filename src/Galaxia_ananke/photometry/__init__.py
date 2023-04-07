#!/usr/bin/env python
"""
Docstring
"""
from .Photometry import Photometry

__all__ = ['available_photo_systems']

available_photo_systems = Photometry()


if __name__ == '__main__':
    pass
