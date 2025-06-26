#!/usr/bin/env python
"""
Docstring
"""
from __future__ import annotations
from typing import TYPE_CHECKING, List
from functools import cached_property

import numpy as np
from astroquery.svo_fps import SvoFps
from astropy.table import Table

from .._constants import *

if TYPE_CHECKING:
    from .InterfaceSvoFpsDriver import InterfaceSvoFpsDriver

__all__ = ['SvoFpsDriver']


class SvoFpsDriver:
    """
    SvoFpsDriver
    """
    _slash_in_isochrone_name: str = '__'
    def __init__(self, interface: InterfaceSvoFpsDriver, instrument_index: int = 0) -> None:
        self.__interface: InterfaceSvoFpsDriver = interface
        self.__instrument_index: int = instrument_index
    
    @property
    def _interface(self) -> InterfaceSvoFpsDriver:
        return self.__interface

    @cached_property
    def __facility_instrument(self) -> List[str]:
        return self._instrument_name.split(self._slash_in_isochrone_name)
    
    @property
    def _instrument_name(self) -> str:
        return self._interface.list_of_instruments[self.__instrument_index]

    @cached_property
    def facility(self) -> str:
        return self.__facility_instrument[0]
    
    @cached_property
    def instrument(self) -> str:
        if len(self.__facility_instrument) > 1:
            return self.__facility_instrument[1]

    @cached_property
    def svo_filter_list(self) -> Table:
        # TODO ideally the following would be enough, but SVO seems to not be consistent with instrument
        # facility_filter_list = SvoFps.get_filter_list(facility=self.facility,
        #                                               instrument=self.instrument)
        facility_filter_list: Table = SvoFps.get_filter_list(facility=self.facility)
        if self.instrument is not None:
            facility_filter_list = facility_filter_list[(facility_filter_list.to_pandas().filterID.str.split('[/.]').str[1] == self.instrument).tolist()]
        _temp = facility_filter_list.to_pandas().filterID.str.split('[/.]')
        facility_filter_list['Instrument'] = np.ma.array(_temp.str[1], dtype='object')
        facility_filter_list['Band'] = np.ma.array(_temp.str[2], dtype='object')
        lower_case_bands = [band.lower() for band in facility_filter_list['Band']]
        for mag_name in set(self._interface.mag_names) - set(facility_filter_list['Band']):
            if mag_name.lower() in lower_case_bands:
                lower_case_bands.index(mag_name.lower())
                facility_filter_list.add_row(facility_filter_list[lower_case_bands.index(mag_name.lower())])
            else:
                facility_filter_list.add_row(len(facility_filter_list.keys())*[None])
            facility_filter_list[-1][['Instrument','Facility']] = facility_filter_list[0][['Instrument','Facility']]
            facility_filter_list[-1][['Band']] = [mag_name]
        facility_filter_list = facility_filter_list[facility_filter_list.to_pandas().Band.reset_index().set_index('Band').loc[self._interface.mag_names].to_numpy().T[0]]
        return facility_filter_list


if __name__ == '__main__':
    raise NotImplementedError()
