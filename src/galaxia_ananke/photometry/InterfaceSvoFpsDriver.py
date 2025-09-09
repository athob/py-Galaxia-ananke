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
Docstring
"""
from __future__ import annotations
from typing import TYPE_CHECKING, List, Dict
from functools import cached_property

from .._constants import *
from .SvoFpsDriver import SvoFpsDriver

if TYPE_CHECKING:
    from .PhotoSystem import PhotoSystem

__all__ = ['InterfaceSvoFpsDriver']


class InterfaceSvoFpsDriver:
    """
    InterfaceSvoFpsDriver
    """
    _plus_in_isochrone_name: str = '+'
    _special_cases: Dict[str, str]  = {
        'DCMC': 'DENIS+2MASS',
        'Euclid': 'Euclid__NISP0+Euclid__VIS',
        'GAIA__0': 'GAIA__GAIA0',
        'GAIA__DR2': 'GAIA__GAIA2r',
        'Generic__Johnson_UBVRIJHK': 'Generic__Johnson_UBVRIJHKL',
        'Generic__Stroemgren': 'Generic__Stromgren',
        'LSST_DP0': 'LSST',
        'Pan_STARRS1': 'PAN-STARRS',
        'HST__WFC3': 'HST__WFC3_UVIS1+HST__WFC3_UVIS2+HST__WFC3_IR'
        }
    def __init__(self, photosystem: PhotoSystem) -> None:
        self.__photosystem: PhotoSystem = photosystem
        self.__svofpsdrivers: List[SvoFpsDriver] = [SvoFpsDriver(self, i) for i in range(len(self.list_of_instruments))]
    
    @classmethod
    def _special_cases_formatter(cls, name: str) -> str:
        return cls._special_cases.get(name, name)

    @property
    def _photosystem(self) -> PhotoSystem:
        return self.__photosystem
    
    @cached_property
    def list_of_instruments(self) -> List[str]:
        return self._plus_in_isochrone_name.join(
            map(self._special_cases_formatter,
                self._special_cases_formatter(self.name).split(self._plus_in_isochrone_name))
                ).split(self._plus_in_isochrone_name)
    
    @property
    def svofpsdrivers(self) -> List[SvoFpsDriver]:
        return self.__svofpsdrivers

    @property
    def name(self) -> str:
        return self._photosystem.name
    
    @property
    def mag_names(self) -> List[str]:
        return self._photosystem.mag_names


if __name__ == '__main__':
    raise NotImplementedError()
