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
The photometric systems: the isochrones and the magnitudes of a system, with what the SVO Filter Profile Service knows of its
filters.
"""
from typing import Any, Callable, Union, List, Dict, Optional
from numpy.typing import NDArray
from functools import cached_property
import re
import warnings

import numpy as np
from astropy import table, units

from .._constants import *
from .._defaults import *
from ..utils import compare_given_and_required, tuple_itemgetter
from .Isochrone import Isochrone
from .InterfaceSvoFpsDriver import InterfaceSvoFpsDriver

__all__ = ['PhotoSystem']


class PhotoSystem:
    """
    A photometric system: its isochrones and the magnitudes that they give, with what the SVO Filter Profile Service knows
    of its filters.

    A system is told by its category and its name, which make its key ``category/name``. The name of each magnitude in the
    catalogue is made of the name of the system and the one of the magnitude (see ``to_export_keys``).
    """
    _required_cmd_magnames_dictkeys = {'magnitude', 'color_minuend', 'color_subtrahend'}
    _ab_zeropoint = 3631 * units.Jy  # by definition of the AB system
    def __init__(self, isochrone: Isochrone) -> None:
        """
        Make the photometric system of an isochrone.

        Parameters
        ----------
        isochrone : Isochrone
            The isochrone, whose category and name are those of the system.
        """
        self.__isochrone = isochrone
        self.__interface_svofps = InterfaceSvoFpsDriver(self)
    
    def __repr__(self) -> str:
        """
        Give the representation of the system.

        Returns
        -------
        representation : str
            The class, with the category and the name of the system.
        """
        cls = self.__class__.__name__
        description = ', '.join([(f"{prop}={getattr(self, prop)}") for prop in ['category', 'name']])
        return f'{cls}({description})'
    
    def check_cmd_magnames(self, cmd_magnames: Union[str,Dict[str,str]]) -> str:
        """
        Check that the magnitude and the colors of a selection in the color-magnitude diagram are magnitudes of the system.

        Parameters
        ----------
        cmd_magnames : str or dict
            The names, either as ``"band1,band2-band3"`` (the magnitude, and the two filters whose difference is the color),
            or as a dictionary with the keys ``magnitude``, ``color_minuend`` and ``color_subtrahend``.

        Returns
        -------
        cmd_magnames : str
            The names as ``"band1,band2-band3"``.

        Raises
        ------
        ValueError
            If it is neither a string nor a dictionary, if the keys of the dictionary are not the three, or if a name is not
            a magnitude of the system: the message gives those that are.
        """
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
    
    def _extract_from_svo_filter_list(self, key: str, unit: Optional[units.UnitBase] = None,
                                      aggregate: Callable[[Any], Any] = lambda stack: np.nanmean(stack, axis=1)
                                      ) -> units.Quantity:
        """
        Give, for each magnitude of the system, a quantity of the SVO service for its filter.

        The service is asked for each instrument of the system. A filter that no instrument has gives NaN, and a filter that
        several have gives the aggregate of what they say.

        Parameters
        ----------
        key : str
            The column of the table of the service that holds the quantity.
        unit : astropy.units.UnitBase, optional
            The unit to give the quantity in. Default to the unit of the column of the service, if it has one.
        aggregate : callable, optional
            What gives the quantity of a filter from what each instrument says of it: it takes the array of the instruments,
            one row for each magnitude, and reduces the instruments. Default to the mean that leaves out the NaN.

        Returns
        -------
        quantity : astropy.units.Quantity
            The quantity of each magnitude, in the order of ``mag_names``.
        """
        parts = [
            None if sfd.svo_filter_list is None else (  # None: the service has no filter for the instrument (with a warning)
                _colu:= sfd.svo_filter_list[key],  # ignored
                _unit:= getattr(_colu, 'unit', None) if unit is None else unit,  # ignored
                _colu.value if _unit is None or _colu.dtype.hasobject else _colu.to(_unit)  # returned
                )[-1]
            for sfd in self._interface_svofps.svofpsdrivers
            ]
        if unit is None:  # the unit of the filters that were found, for those that were not
            unit = next((getattr(part, 'unit', None) for part in parts if part is not None), None)
        nothing = np.full(len(self.mag_names), np.nan) * (1 if unit is None else unit)
        stack = np.vstack([nothing if part is None else part for part in parts]).T
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)  # the mean of nothing, for a filter that no instrument has: NaN
            return aggregate(stack)

    @property
    def _isochrone(self) -> Isochrone:
        """The isochrone of the system."""
        return self.__isochrone

    @property
    def _interface_svofps(self) -> InterfaceSvoFpsDriver:
        """The interface of the system with the SVO Filter Profile Service."""
        return self.__interface_svofps

    @property
    def category(self) -> str:
        """The category of the system: the one of its isochrone."""
        return self._isochrone.category
    
    @property
    def name(self) -> str:
        """The name of the system: the one of its isochrone."""
        return self._isochrone.name

    @property
    def categ_and_name(self) -> List[str]:
        """The category and the name of the system, as a list."""
        return [self.category, self.name]

    @property
    def key(self) -> str:
        """The key of the system: its category and its name, with a slash between."""
        return f"{self.category}/{self.name}"

    @cached_property
    def mag_names(self) -> List[str]:
        """The names of the magnitudes of the system, as the descriptor of its isochrones gives them."""
        return self._isochrone.mag_names
    
    @cached_property
    def effective_wavelengths(self) -> units.Quantity:
        """The effective wavelength of each filter in the order of the magnitudes (NaN for a filter that SVO has none of)."""
        return self._extract_from_svo_filter_list('WavelengthEff', unit=DEF_UNIT.wavelength)

    @cached_property
    def zeropoints(self) -> units.Quantity:
        """The zero point of each filter in the order of the magnitudes (NaN for a filter that SVO has none of)."""
        return self._extract_from_svo_filter_list('ZeroPoint', unit=units.Jy).to(
            DEF_UNIT.spectral, equivalencies=units.spectral_density(self.effective_wavelengths)
            )

    @property
    def magnitude_system(self) -> Optional[str]:
        """The magnitude system ('Vega', 'AB' or 'ST') of the magnitudes, as the isochrone descriptor says it, or None."""
        return self._isochrone.magnitude_system

    @cached_property
    def _vega_zeropoints(self) -> units.Quantity:
        """
        In Jy, in the order of the magnitudes: the zero point that SVO gives of each filter when it is that of the Vega
        system (the mean over the instruments that have the filter, as for zeropoints), NaN for a filter that SVO has no
        such zero point of.
        """
        columns = []
        for sfd in self._interface_svofps.svofpsdrivers:
            filters = sfd.svo_filter_list
            if filters is None:  # the service has no filter for the instrument (with a warning)
                columns.append(np.full(len(self.mag_names), np.nan))
                continue
            frame = filters.to_pandas()
            to_jansky = (filters['ZeroPoint'].unit or units.Jy).to(units.Jy)
            is_vega = (frame['MagSys'] == 'Vega').to_numpy()
            columns.append(np.where(is_vega, frame['ZeroPoint'].to_numpy(dtype=float) * to_jansky, np.nan))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)  # a filter that no instrument has: the mean of nothing
            return np.nanmean(np.vstack(columns).T, axis=1) * units.Jy

    @cached_property
    def vega_to_ab_offsets(self) -> Dict[str, float]:
        """
        What to add to a magnitude of the Vega system to have it in the AB system, ``m_AB = m_Vega + offset``, for each
        magnitude of the system by its key in the output (``to_export_keys``): ``2.5 log10(3631 Jy / zero point of Vega)``,
        with the zero point that the SVO filter profile service gives of the filter (asked of the service once, which needs
        the network). It is as good as the match between the filter of SVO and the one that the isochrones were made with,
        and it is NaN, with a warning, for a filter that SVO has no Vega zero point of. Whether the magnitudes of the system
        are in the Vega system is said by ``magnitude_system``: see ``offsets_to_magnitude_system`` for what to add
        whatever it is.
        """
        zeropoints = self._vega_zeropoints.to_value(units.Jy)
        with np.errstate(divide='ignore', invalid='ignore'):
            offsets = 2.5 * np.log10(self._ab_zeropoint.to_value(units.Jy) / zeropoints)
        missing = [mag_name for mag_name, offset in zip(self.mag_names, offsets) if np.isnan(offset)]
        if missing:
            warnings.warn(f"SVO has no zero point of the Vega system for the filters {missing} of {self.key}: "
                          f"their offsets are NaN", stacklevel=3)
        return dict(zip(self.to_export_keys, offsets.tolist()))

    def offsets_to_magnitude_system(self, target: str) -> Dict[str, float]:
        """
        Give what to add to each magnitude of the system to have it in a magnitude system.

        It is zero for the magnitudes that are in it already, whichever system they are in (``magnitude_system``). The
        conversion between Vega and AB is the one of ``vega_to_ab_offsets``, which needs the network, but nothing is asked
        of the service when there is nothing to convert.

        Parameters
        ----------
        target : str
            The magnitude system to convert to: 'Vega' or 'AB', whatever the case.

        Returns
        -------
        offsets : dict
            What to add to each magnitude, by its key in the catalogue (``to_export_keys``).

        Raises
        ------
        ValueError
            If the target is not given, is not a magnitude system, or if the descriptor of the isochrones does not say which
            system the magnitudes are in.
        NotImplementedError
            If the conversion is not one between Vega and AB.
        """
        checked_target = Isochrone._check_magnitude_system(target)
        if checked_target is None:
            raise ValueError("A magnitude system to convert to is needed")
        target = checked_target
        native = self.magnitude_system
        if native is None:
            raise ValueError(f"The magnitude system of {self.key} is not stated (its descriptor has no line "
                             f"'# Magnitude system: ...'): nothing can be converted")
        if native == target:
            return dict.fromkeys(self.to_export_keys, 0.0)
        if {native, target} != {'Vega', 'AB'}:
            raise NotImplementedError(f"Only the conversions between the Vega and the AB systems are available, "
                                      f"not from {native} to {target}")
        sign = 1 if native == 'Vega' else -1
        return {key: sign * offset for key, offset in self.vega_to_ab_offsets.items()}

    @cached_property
    def to_export_keys(self) -> List[str]:
        """
        The names of the magnitudes in the catalogue: the name of the system and the one of the magnitude,
        in lower case, joined by an underscore.
        """
        return [f"{self.name.lower()}_{mag_name.lower()}" for mag_name in self.mag_names]
    
    @cached_property
    def mag_names_to_export_keys_mapping(self) -> Dict[str,str]:
        """The name of each magnitude in the catalogue, by its name in the isochrones."""
        return dict(zip(self.mag_names, self.to_export_keys))

    def get_nearest_isochrone_track_for_age_and_metallicity(self, ages: NDArray[np.floating[Any]],
                                                            metallicities: NDArray[np.floating[Any]]) -> List[table.QTable]:
        """
        Give, for each age and metallicity, the table of the isochrone at the nearest age of the nearest metallicity.

        The nearest is the one of the values that the isochrone has, the limit between two of them being halfway.

        Parameters
        ----------
        ages : numpy.ndarray
            The ages, in the decimal logarithm of the years, one for each star.
        metallicities : numpy.ndarray
            The metallicities, as the decimal logarithm of their ratio to the solar one, one for each star.

        Returns
        -------
        tables : list of astropy.table.QTable
            The table of each star, in the order that was given.
        """
        z_keys = np.array(list(self._isochrone.qtables_unique_ages_dictionary.keys()))
        nearest_z_keys = z_keys[np.digitize(metallicities, (z_keys[:-1]+z_keys[1:])/2)]
        nearest_z_sorter = np.argsort(nearest_z_keys)
        unique_z_keys, index_z_keys = np.unique(nearest_z_keys[nearest_z_sorter], return_index=True)
        split_per_nearest_z_keys = tuple_itemgetter(*tuple(np.split(nearest_z_sorter, index_z_keys[1:])))
        nearest_ages_per_z_keys = {z_key: age_keys[np.digitize(sub_ages, (age_keys[:-1]+age_keys[1:])/2)]
                                   for z_key, sub_ages in zip(unique_z_keys, split_per_nearest_z_keys(ages))
                                   if (age_keys:=np.array(list(self._isochrone.qtables_unique_ages_dictionary[z_key]))) is not None}
        nearest_ages_keys = np.hstack(list(nearest_ages_per_z_keys.values()))[np.argsort(nearest_z_sorter)]
        return [self._isochrone.qtables_dictionary[z][age] for z,age in zip(nearest_z_keys,nearest_ages_keys)]


if __name__ == '__main__':
    raise NotImplementedError()
