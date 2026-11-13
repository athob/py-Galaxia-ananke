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
from typing import TYPE_CHECKING, Tuple, List, Dict
from numpy.typing import ArrayLike, NDArray
from functools import cached_property
from operator import itemgetter
import dataclasses as dc

import numpy as np
import pandas as pd
from astropy import table, units, constants

from .._constants import *
from ..utils import classproperty

if TYPE_CHECKING:
    from .Isochrone import Isochrone
    from .IsochroneFile import IsochroneFile

__all__ = ['python_formatting', 'cmd_2_1_formatting', 'cmd_2_8_formatting', 'cmd_3_0_formatting', 'cmd_3_6_formatting']


units.Zsun = units.def_unit('Zsun') # TODO could define PhysicalType "abundance" ?
# CMD2-label: units.integer = units.Quantity(1, dtype='int') # units.def_unit('integer', units.Quantity(1, dtype='int')) # TODO pull issue in astropy for integer dtype unit?


class _BaseFormatting:
    _zini = 'Z_ini'
    _unit_zini = units.LogUnit(units.Zsun)
    _age = 'Age'
    _unit_age = units.LogUnit(units.yr)
    _mini = 'M_ini'
    _mact = 'M_act'
    _unit_mass = units.Msun
    _radi = 'Radius'
    _unit_radi = units.Rsun
    _lum = 'Lum'
    _unit_lum = units.LogUnit(units.Lsun)
    _teff = 'T_eff'
    _unit_temp = units.LogUnit(units.Kelvin)
    _grav = 'Grav'
    _unit_grav = units.LogUnit(units.Gal)
    # CMD2-label: _label = 'Label'
    # CMD2-label: _unit_label = units.integer
    @classproperty
    def units_mapper(cls):
        return {
            cls._zini: cls._unit_zini,
            cls._age: cls._unit_age,
            cls._mini: cls._unit_mass,
            cls._mact: cls._unit_mass,
            cls._radi: cls._unit_radi,
            cls._lum: cls._unit_lum,
            cls._teff: cls._unit_temp,
            cls._grav: cls._unit_grav,
            # CMD2-label: cls._label: cls._unit_label
        }

@dc.dataclass()
class Formatting(_BaseFormatting):
    """
    Formatting object: where the columns that the package needs are in the files of the isochrones of a photometric system. The files of
    the systems were made at different times by different versions of the CMD service, so their layouts differ. The magnitudes are not
    part of it: the descriptor of the system says where they start.

    The label is retired for now: the column that the layouts took for it is a label in the files of CMD 3.x, but the bolometric magnitude
    in the files of CMD 2.x (and so for the default layout), and neither the package nor the backend uses it. Its logic is commented out,
    here, in the layouts below and in the tests, with the mark `CMD2-label`: uncomment it, and drop what replaces it, once the systems of
    CMD 2.x are retired (`grep -rn CMD2-label src/` lists it all).
    """
    i_age: int
    i_mini: int
    i_mact: int
    i_lum: int
    i_teff: int
    i_grav: int
    # CMD2-label: i_label: int
    lin_age: bool = False
    feh_sun: float = 0.0152  # PARSEC (Bressan et al. 2012) and later tracks; 0.019 for Marigo et al. (2008) and previous ones (input form of CMD)

    @cached_property
    def format_mapper(self) -> Dict[int,str]:
        return {
            self.i_age: self._age,
            self.i_mini: self._mini,
            self.i_mact: self._mact,
            self.i_lum: self._lum,
            self.i_teff: self._teff,
            self.i_grav: self._grav,
            # CMD2-label: self.i_label: self._label
            }

    @cached_property
    def nonmagnitude_itemgetter(self) -> itemgetter:
        return itemgetter(*self.format_mapper.keys())

    def metallicity_converter(self, metallicity: ArrayLike) -> NDArray:
        return np.log10(metallicity/np.array(self.feh_sun))

    def qtable_per_age_from_isochronefile(self, iso_file: IsochroneFile) -> table.QTable:
        isochrone: Isochrone = iso_file._isochrone
        units_mapper = self.units_mapper
        column_names: Tuple[str] = self.nonmagnitude_itemgetter(iso_file.data.keys()) + isochrone.magnitude_itemgetter(iso_file.data.keys())
        final_names: List[str] = list(self.format_mapper.values()) + isochrone.mag_names
        column_mapper: Dict[str,str] = dict(zip(column_names, final_names))
        reduced_data: pd.DataFrame = iso_file.data[column_names].to_pandas().rename(columns=column_mapper).reset_index().sort_values([self._age, self._mini, 'index']).drop('index', axis=1)
        reduced_data.insert(0, self._zini, self.metallicity_converter(iso_file.metallicity))
        reduced_data.insert(5, self._radi, np.sqrt(constants.G*reduced_data[self._mact].to_numpy()*self._unit_mass/(reduced_data[self._grav].to_numpy()*self._unit_grav).to(self._unit_grav.physical_unit)).to(self._unit_radi).value)
        if self.lin_age:
            reduced_data[self._age] = np.log10(reduced_data[self._age])
        return {
            age: table.QTable({
                key: units.Magnitude(value)
                if key in isochrone.mag_names
                else value*units_mapper[key]
                for key,value in df.to_dict(orient='list').items()
                })
            for age,df in reduced_data.groupby(self._age)
            }


# The layouts of the files of the isochrones, named after the CMD version that writes the files (the first line of every file says it, as in
# "# File generated by CMD 3.7"). The columns that each layout takes are the ones that the backend reads (IsochroneBase::readfile of
# Galaxia-ananke), and the tests pin them, with the solar metallicity: 0.019 for the models of Marigo et al. (2008) and previous ones, 0.0152
# for PARSEC (Bressan et al. 2012) and later ones, as the input form of CMD says.
# CMD2-label: python_formatting = Formatting(0, 1, 2, 3, 4, 5, 6, False)
# CMD2-label: cmd_2_1_formatting = Formatting(0, 1, 2, 3, 4, 5, 6, False, feh_sun=0.019)
# CMD2-label: cmd_2_8_formatting = Formatting(1, 2, 3, 4, 5, 6, 7, False)
# CMD2-label: cmd_3_0_formatting = Formatting(1, 2, 3, 4, 5, 6, 7, True)
# CMD2-label: cmd_3_6_formatting = Formatting(2, 3, 5, 6, 7, 8, 9, False)

# What is made in Python (`Isochrone(name, data)`), and what a system that is not known is read with: Age (in log), M_ini, M_act, Lum, T_eff,
# Grav, then the magnitudes
python_formatting = Formatting(0, 1, 2, 3, 4, 5, False)
# CMD 2.1 and 2.2, the models of Marigo et al. (2008), for APO__SDSS, DCMC, Generic__Johnson_UBVRIJHK and Generic__Stroemgren: log(age/yr),
# M_ini, M_act, logL, logTe, logG, mbol, the magnitudes, then C/O, M_hec, period, pmode, logMdot and int_IMF. The columns of python_formatting
# and Zsun = 0.019
cmd_2_1_formatting = Formatting(0, 1, 2, 3, 4, 5, False, feh_sun=0.019)
# CMD 2.8, PARSEC v1.2S + COLIBRI PR16, for WFIRST+HST__WFC3: Z, log(age/yr), M_ini, M_act, logL, logTe, logG, mbol, the magnitudes, then C/O,
# M_hec, period, pmode, logMdot, slope and int_IMF
cmd_2_8_formatting = Formatting(1, 2, 3, 4, 5, 6, False)
# CMD 3.0, PARSEC v1.2S + COLIBRI PR16, for WFIRST, LSST, GAIA__0, GAIA__DR2, HST__WFC3 and Pan_STARRS1: Zini, Age (linear, in years), Mini,
# Mass, logL, logTe, logg, label, McoreTP, C_O, period0, period1, pmode, Mloss, tau1m, X, Y, Xc, Xn, Xo, Cexcess, Z, mbolmag, the magnitudes
cmd_3_0_formatting = Formatting(1, 2, 3, 4, 5, 6, True)
# CMD 3.6 and 3.7, PARSEC v1.2S + COLIBRI S_37 + S_35 + PR16, for CTIO__DECam, LSST_DP0, Roman, Euclid, JWST__NIRCam and JWST__MIRI: Zini, MH,
# logAge, Mini, int_IMF, Mass, logL, logTe, logg, label, McoreTP, C_O, period0 to period4, pmode, Mloss, tau1m, X, Y, Xc, Xn, Xo, Cexcess, Z,
# mbolmag, the magnitudes
cmd_3_6_formatting = Formatting(2, 3, 5, 6, 7, 8, False)

if __name__ == '__main__':
    raise NotImplementedError()
