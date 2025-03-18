#!/usr/bin/env python
"""
Docstring
"""
from typing import Union, List, Dict
from numpy.typing import ArrayLike
from functools import cached_property
from operator import itemgetter
import re

import numpy as np
from astropy import table, units

from .._constants import *
from .._defaults import *
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
    
    def check_cmd_magnames(self, cmd_magnames: Union[str,Dict[str,str]]) -> str:
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
    
    def _extract_from_svo_filter_list(self, key: str, unit=None, aggregate=lambda stack: np.nanmean(stack, axis=1)):
        stack = np.vstack(
            [
                (
                _colu:= sfd.svo_filter_list[key],  # ignored
                _unit:= getattr(_colu, 'unit', None) if unit is None else unit,  # ignored
                _colu.value if _unit is None or _colu.dtype.hasobject else _colu.to(_unit)  # returned
                )[-1]
            for sfd in self._interface_svofps.svofpsdrivers
            ]).T
        return aggregate(stack)

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
    def categ_and_name(self):
        return [self.category, self.name]

    @property
    def key(self):
        return f"{self.category}/{self.name}"

    @cached_property
    def mag_names(self):
        return self._isochrone.mag_names
    
    @cached_property
    def effective_wavelengths(self):
        return self._extract_from_svo_filter_list('WavelengthEff', unit=DEF_UNIT.wavelength)

    @cached_property
    def zeropoints(self):
        return self._extract_from_svo_filter_list('ZeroPoint').to(
            DEF_UNIT.spectral, equivalencies=units.spectral_density(self.effective_wavelengths)
            )

    @cached_property
    def to_export_keys(self):
        return [f"{self.name.lower()}_{mag_name.lower()}" for mag_name in self.mag_names]

    def get_nearest_isochrone_track_for_age_and_metallicity(self, ages: ArrayLike, metallicities: ArrayLike) -> List[table.QTable]:
        z_keys = np.array(list(self._isochrone.qtables_unique_ages_dictionary.keys()))
        nearest_z_keys = z_keys[np.digitize(metallicities, (z_keys[:-1]+z_keys[1:])/2)]
        nearest_z_sorter = np.argsort(nearest_z_keys)
        unique_z_keys, index_z_keys = np.unique(nearest_z_keys[nearest_z_sorter], return_index=True)
        split_per_nearest_z_keys = itemgetter(*tuple(np.split(nearest_z_sorter, index_z_keys[1:])))
        nearest_ages_per_z_keys = {z_key: age_keys[np.digitize(sub_ages, (age_keys[:-1]+age_keys[1:])/2)]
                                   for z_key, sub_ages in zip(unique_z_keys, split_per_nearest_z_keys(ages))
                                   if (age_keys:=np.array(list(self._isochrone.qtables_unique_ages_dictionary[z_key]))) is not None}
        nearest_ages_keys = np.hstack(list(nearest_ages_per_z_keys.values()))[np.argsort(nearest_z_sorter)]
        return [self._isochrone.qtables_dictionary[z][age] for z,age in zip(nearest_z_keys,nearest_ages_keys)]


if __name__ == '__main__':
    pass
