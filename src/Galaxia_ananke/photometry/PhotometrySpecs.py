#!/usr/bin/env python
"""
Docstring
"""
from . import available_photo_systems
from ..constants import *

__all__ = ['PhotometrySpecs']


class PhotometrySpecs:
    def __init__(self, photo_sys=DEFAULT_PSYS, cmd_magnames=DEFAULT_CMD, cmd_box=DEFAULT_CMD_BOX) -> None:
        self._set_isochrones_from_photosys(photo_sys)
    
    def _set_isochrones_from_photosys(self, photo_sys):
        if isinstance(photo_sys, str):
            photo_sys = [photo_sys]
        self.__isochrones = [available_photo_systems[psys] for psys in photo_sys]
    
    def _check_cmd_magnames(cmd_magnames):
        pass
        # cmd_magnames = isochrone.check_cmd_magnames(cmd_magnames)
        
    @property
    def isochrones(self):
        return self.__isochrones

    @property
    def photo_sys(self):
        return {isochrone.key for isochrone in self.__isochrones}
