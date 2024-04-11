#!/usr/bin/env python
"""
Docstring
"""
import re

from ..constants import *
from ..utils import compare_given_and_required
from .Isochrone import Isochrone
from .InterfaceSvoFpsDriver import InterfaceSvoFpsDriver

__all__ = ['PhotoSystem']


class PhotoSystem:
    """
    PhotoSystem
    """
    _required_cmd_magnames_dictkeys = {'magnitude', 'color_minuend', 'color_subtrahend'}
    def __init__(self, isochrone: Isochrone) -> None:
        self.__isochrone = isochrone
        self.__interface_svofps = InterfaceSvoFpsDriver(self)
    
    def __repr__(self) -> str:
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['category', 'name']])
        return f'{cls}({description})'
    
    def check_cmd_magnames(self, cmd_magnames):
        if isinstance(cmd_magnames, str):
            check = set(re.split('[ ,-]', cmd_magnames))
        elif isinstance(cmd_magnames, dict):
            compare_given_and_required(cmd_magnames.keys(), self._required_cmd_magnames_dictkeys, error_message="Given cmd_magnames dict covers wrong set of keys")
            check = set(cmd_magnames.values())
            cmd_magnames = f"{cmd_magnames['magnitude']},{cmd_magnames['color_minuend']}-{cmd_magnames['color_subtrahend']}"
        else:
            raise ValueError(f"cmd_magnames should be either a string or a dict: given {type(cmd_magnames)} instead")
        if not check.issubset(set(self.mag_names)):
            raise ValueError(f"CMD magnitude names '{cmd_magnames}' don't match isochrone names: choose among {self.mag_names}")
        else:
            return cmd_magnames

    @property
    def _isochrone(self):
        return self.__isochrone

    @property
    def _interface_svofps(self):
        return self.__interface_svofps

    @property
    def category(self):
        return self._isochrone.category
    
    @property
    def name(self):
        return self._isochrone.name

    @property
    def key(self):
        return f"{self.category}/{self.name}"

    @property
    def mag_names(self):
        return self._isochrone.mag_names

    @property
    def to_export_keys(self):
        return [f"{self.name.lower()}_{mag_name.lower()}" for mag_name in self.mag_names]


if __name__ == '__main__':
    pass
