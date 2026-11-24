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
    The link between a photometric system and the SVO Filter Profile Service: it finds the instruments of the system there.

    The name of a photometric system is the one of the folder of its isochrones. It names one instrument of the service, or
    several, joined by ``+``, and the names that CMD gives to some instruments are not the ones of the service (see the
    special cases). The interface makes one ``SvoFpsDriver`` for each instrument.
    """
    _plus_in_isochrone_name: str = '+'
    _special_cases: Dict[str, str]  = {
        'DCMC': 'DENIS+2MASS',
        'Euclid': 'Euclid__NISP0+Euclid__VIS',
        'GAIA__0': 'GAIA__GAIA0',
        'GAIA__DR2': 'GAIA__GAIA2r',
        'Generic__Johnson_UBVRIJHK': 'Generic__Johnson_UBVRIJHKL',
        'Generic__Stroemgren': 'Generic__Stromgren',
        'LSST': 'LSST__LSST',                       # the facility LSST of SVO holds two instruments, LSST and LSSTv19 (added in October 2026, with two filters, g and i)
        'LSST_DP0': 'LSST__LSST',
        'Pan_STARRS1': 'PAN-STARRS',
        'HST__WFC3': 'HST__WFC3_UVIS1+HST__WFC3_UVIS2+HST__WFC3_IR'
        }
    def __init__(self, photosystem: PhotoSystem) -> None:
        """
        Make the interface and the drivers of the instruments of a photometric system.

        Parameters
        ----------
        photosystem : PhotoSystem
            The photometric system.
        """
        self.__photosystem: PhotoSystem = photosystem
        self.__svofpsdrivers: List[SvoFpsDriver] = [SvoFpsDriver(self, i) for i in range(len(self.list_of_instruments))]
    
    @classmethod
    def _special_cases_formatter(cls, name: str) -> str:
        """
        Give the name that the service has for an instrument that CMD names differently.

        Parameters
        ----------
        name : str
            The name, as the isochrones have it.

        Returns
        -------
        name : str
            The name of the service, or the one given if it is the same.
        """
        return cls._special_cases.get(name, name)

    @property
    def _photosystem(self) -> PhotoSystem:
        """The photometric system that this is the interface of."""
        return self.__photosystem
    
    @cached_property
    def list_of_instruments(self) -> List[str]:
        """
        The names that the service has for the instruments of the system, in the order of its name.

        The special cases are applied to the whole name and to each of its parts, and the result may have more instruments
        than the name has parts.
        """
        return self._plus_in_isochrone_name.join(
            map(self._special_cases_formatter,
                self._special_cases_formatter(self.name).split(self._plus_in_isochrone_name))
                ).split(self._plus_in_isochrone_name)
    
    @property
    def svofpsdrivers(self) -> List[SvoFpsDriver]:
        """The drivers of the instruments of the system, in the order of ``list_of_instruments``."""
        return self.__svofpsdrivers

    @property
    def name(self) -> str:
        """The name of the photometric system."""
        return self._photosystem.name
    
    @property
    def mag_names(self) -> List[str]:
        """The names of the magnitudes of the photometric system."""
        return self._photosystem.mag_names


if __name__ == '__main__':
    raise NotImplementedError()
