#!/usr/bin/env python
"""
Docstring
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import pathlib
import itertools
import h5py as h5
import ebf
import vaex
from astropy import units, coordinates

from .constants import *

if TYPE_CHECKING:
    from . import Survey

__all__ = ['Output']


def shift_g_lon(lon): # restrict longitude values to be within (-180,180)
    return -((-lon+180)%360-180)


class Output:
    _pos = ['px', 'py', 'pz']  # positions
    _vel = ['vx', 'vy', 'vz']  # velocities
    _cel = ['ra', 'dec']  # celestial coordinates
    _gal = ['glon', 'glat']  # galactic coordinates
    _rad = 'rad'  # distance
    _dmod = 'dmod'  # distance modulus
    _mtip = 'mtip'  # mass at tip of giant branch for star of given age & metallicity (Msun)
    _mact = 'mact'  # current stellar mass (Msun)
    _mini = 'smass'  # ZAMS stellar mass (Msun)
    _age = 'age'  # age (Gyr)
    _grav = 'grav'  # surface gravity (log g)
    _feh = 'feh'  # [fe/h]
    _teff = 'teff'  # temperature (returned log10(teff/K) by Galaxia)
    _lum = 'lum'  # luminosity (returned log10(lum/lsun) by Galaxia)
    _parentid = 'parentid'  # id of parent particle 
    _partid = 'partid'  # flag of particle's central star
    #####
    _pi = 'pi'  # parallax
    _mu = ['mura', 'mudec']  # proper motions
    _mugal = ['mul', 'mub']  # galactic proper motions
    _vr = 'vr'  # radial velocity
    #####
    _export_keys = _pos + _vel + _cel + _gal + [_rad, _dmod, _teff, _lum, _grav, _mtip, _mact, _mini, _age, _feh, _parentid, _partid]
    _postprocess_keys = [_pi] + _mu + _mugal + [_vr]
    _vaex_under_list = ['_repr_html_']
    def __init__(self, survey: Survey, parameters: dict) -> None:
        self.__survey = survey
        self.__parameters = parameters
        self.__vaex = None
        self.__path = None

    def __dir__(self):
        return sorted({i for i in self.__vaex.__dir__() if not i.startswith('_')}.union(
            super(Output, self).__dir__()).union(
            self._vaex_under_list if self.__vaex is not None else []))

    def __repr__(self) -> str:
        return repr(self._vaex)
    
    def __getitem__(self, item):
        return self._vaex[item]
    
    def __setitem__(self, item, value):
        self._vaex[item] = value
    
    def __getattr__(self, item):
        if (item in self.__vaex.__dir__() and not item.startswith('_')) or (item in self._vaex_under_list and self.__vaex is not None):
            return getattr(self.__vaex, item)
        else:
            return self.__getattribute__(item)

    @classmethod
    def _compile_export_mag_names(cls, isochrones):
        return list(itertools.chain.from_iterable([isochrone.to_export_keys for isochrone in isochrones]))
    
    @classmethod
    def _make_export_keys(cls, isochrones, extra_keys=[]):
        return cls._export_keys + extra_keys + cls._compile_export_mag_names(isochrones)

    @classmethod
    def _make_catalogue_keys(cls, isochrones, extra_keys=[]):
        return cls._make_export_keys(isochrones, extra_keys=cls._postprocess_keys+extra_keys)

    def _make_input_optional_keys(self):
        return [k if k != 'id' else 'satid' for k in self.survey.input.optional_keys()]

    def _ebf_to_hdf5(self):
        hdf5_file = self._hdf5
        with h5.File(hdf5_file, 'w') as f5:
            for k in self.export_keys:
                f5.create_dataset(name=k, data=ebf.read(self._ebf, f"/{k}"))
            print(f"Exported the following quantities to {hdf5_file}")
            print(list(f5.keys()))
        self.__vaex = vaex.open(hdf5_file)

    def _post_process(self):
        self._pp_convert_cartesian_to_galactic()
        self._pp_convert_galactic_to_icrs()
        self._vaex[self._pi] = 1.0/self._vaex[self._rad]  # parallax in mas (from distance in kpc)
        self._vaex[self._teff] = 10**self._vaex[self._teff]  #Galaxia returns log10(teff/K)
        self._vaex[self._lum] = 10**self._vaex[self._lum]  #Galaxia returns log10(lum/lsun)
        self.flush_extra_columns_to_hdf5(with_columns=[self._teff, self._lum])

    def _pp_convert_cartesian_to_galactic(self):
        """
        converts positions & velocities from mock catalog Cartesian coordinates (relative to solar position) 
        into Galactic coordinates, assuming Sun is on -x axis (use rotateStars)
        """
        gc = coordinates.Galactic(u = self._vaex[self._pos[0]].to_numpy()*units.kpc,
                                  v = self._vaex[self._pos[1]].to_numpy()*units.kpc,
                                  w = self._vaex[self._pos[2]].to_numpy()*units.kpc,
                                  U = self._vaex[self._vel[0]].to_numpy()*units.km/units.s,
                                  V = self._vaex[self._vel[1]].to_numpy()*units.km/units.s,
                                  W = self._vaex[self._vel[2]].to_numpy()*units.km/units.s,
                                  representation_type = coordinates.CartesianRepresentation,
                                  differential_type   = coordinates.CartesianDifferential)
        self._vaex[self._gal[0]] = shift_g_lon(gc.spherical.lon.value)
        self._vaex[self._gal[1]] = gc.spherical.lat.value
        self._vaex[self._rad]    = gc.spherical.distance.value
        ####################################
        self._vaex[self._mugal[0]] = gc.sphericalcoslat.differentials['s'].d_lon_coslat.value
        self._vaex[self._mugal[1]] = gc.sphericalcoslat.differentials['s'].d_lat.value
        self._vaex[self._vr]       = gc.sphericalcoslat.differentials['s'].d_distance.value
        self.flush_extra_columns_to_hdf5(with_columns=self._gal+[self._rad])

    def _pp_convert_galactic_to_icrs(self):
        """
        converts PMs in galactic coordinates (mulcosb, mub) in arcsec/yr (as output by Galaxia)
        to ra/dec in mas/yr (units of output catalog)
        """
        c = coordinates.Galactic(l               = self._vaex[self._gal[0]].to_numpy()*units.degree,
                                 b               = self._vaex[self._gal[1]].to_numpy()*units.degree,
                                 distance        = self._vaex[self._rad].to_numpy()*units.kpc,
                                 pm_l_cosb       = self._vaex[self._mugal[0]].to_numpy()*units.mas/units.yr,
                                 pm_b            = self._vaex[self._mugal[1]].to_numpy()*units.mas/units.yr,
                                 radial_velocity = self._vaex[self._vr].to_numpy()*units.km/units.s)
        c_icrs = c.transform_to(coordinates.ICRS())
        self._vaex[self._cel[0]] = c_icrs.ra.value
        self._vaex[self._cel[1]] = c_icrs.dec.value
        ####################################
        self._vaex[self._mu[0]]  = c_icrs.pm_ra_cosdec.to(units.mas/units.yr).value
        self._vaex[self._mu[1]]  = c_icrs.pm_dec.to(units.mas/units.yr).value
        self.flush_extra_columns_to_hdf5(with_columns=self._cel)
    
    def _pp_convert_icrs_to_galactic(self):
        """
        converts PMs from ICRS coordinates (muacosd, mudec) to Galactic (mul, mub)
        input and output in mas/yr for PMs and degrees for positions
        also exports the galactic lat and longitude
        """
        c = coordinates.ICRS(ra           = self._vaex[self._cel[0]].to_numpy()*units.degree,
                             dec          = self._vaex[self._cel[1]].to_numpy()*units.degree,
                             pm_ra_cosdec = self._vaex[self._mu[0]].to_numpy()*units.mas/units.yr,
                             pm_dec       = self._vaex[self._mu[1]].to_numpy()*units.mas/units.yr)

        c_gal = c.transform_to(coordinates.Galactic())
        self._vaex[self._gal[0]]   = shift_g_lon(c_gal.l.value)
        self._vaex[self._gal[1]]   = c_gal.b.value
        self._vaex[self._mugal[0]] = c_gal.pm_l_cosb.to(units.mas/units.yr).value
        self._vaex[self._mugal[1]] = c_gal.pm_b.to(units.mas/units.yr).value
        self.flush_extra_columns_to_hdf5(with_columns=self._gal+self._mugal)

    def __name_with_ext(self, ext):
        name_base = self._file_base
        return name_base.parent / f"{name_base.name}{ext}"
    
    def save(self, path):
        """
            Save output to new path
        """
        old_path = self._path
        self.__path = pathlib.Path(path)
        self._vaex.close()
        old_path.rename(self._path)
        self.__vaex = vaex.open(self._path)

    @property
    def survey(self):
        return self.__survey
    
    @property
    def isochrones(self):
        return self.survey.isochrones

    @property
    def export_keys(self):
        return self._make_export_keys(self.isochrones, extra_keys=self._make_input_optional_keys())
    
    @property
    def catalogue_keys(self):
        return self._make_catalogue_keys(self.isochrones, extra_keys=self._make_input_optional_keys())
    
    @property
    def output_name(self):
        return f"{self.survey.surveyname}.{self.survey.inputname}"

    @property
    def rsun_skycoord(self):
        _temp = [self._parameters[k] for k in TTAGS.rSun]
        return coordinates.SkyCoord(u=_temp[0], v=_temp[1], w=_temp[2], unit='kpc', representation_type='cartesian', frame='galactic')

    @property
    def _parameters(self):
        return self.__parameters
    
    @property
    def _vaex(self):
        if self.__vaex is None:
            raise RuntimeError("Don't attempt creating an Output object on your own, those are meant to be returned by Survey")
        else:
            return self.__vaex

    @property
    def _path(self):
        if self.__path is None:
            return self._hdf5
        else:
            return self.__path

    @property
    def _file_base(self):
        return GALAXIA_TMP / self.output_name
    
    @property
    def _ebf(self):
        return self.__name_with_ext('.ebf')

    @property
    def _hdf5(self):
        return self.__name_with_ext('.h5')
    
    def flush_extra_columns_to_hdf5(self, with_columns=[]):  # temporary until vaex supports it
        hdf5_file = self._hdf5
        old_column_names = set(vaex.open(hdf5_file).column_names)
        with h5.File(hdf5_file, 'r+') as f5:
            extra_columns = [k for k in set(self.column_names)-old_column_names if not k.startswith('__')]
            for k in extra_columns:
                f5.create_dataset(name=k, data=self[k].to_numpy())
            if extra_columns:
                print(f"Exported the following quantities to {hdf5_file}")
                print(extra_columns)
            for k in with_columns:
                f5[k][...] = self[k].to_numpy()
            if with_columns:
                print(f"Overwritten the following quantities to {hdf5_file}")
                print(with_columns)
        self.__vaex = vaex.open(hdf5_file)


if __name__ == '__main__':
    pass
