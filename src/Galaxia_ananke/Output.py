#!/usr/bin/env python
"""
Contains the Output class definition

Please note that this module is private. The Output class is
available in the main ``Galaxia`` namespace - use that instead.
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
from . import Input

if TYPE_CHECKING:
    from . import Survey

__all__ = ['Output']


def shift_g_lon(lon): # restrict longitude values to be within (-180,180)
    return -((-lon+180)%360-180)


class Output:
    _position_prop = (('px', 'py', 'pz'), "Position coordinates in kpc")
    _velocity_prop = (('vx', 'vy', 'vz'), "Velocity coordinates in km/s")
    _celestial_prop = (('ra', 'dec'), "Celestial equatorial coordinates in degrees")
    _galactic_prop = (('glon', 'glat'), "Celestial galactic coordinates in degrees")
    _distance_prop = ('rad', "Distance in kpc")
    _modulus_prop = ('dmod', "Distance modulus in magnitude units")
    _trgbmass_prop = ('mtip', "Tip of the Red Giant Branch stellar mass in solar masses")
    _currentmass_prop = ('mact', "Current stellar mass in solar masses")
    _zamsmass_prop = ('smass', "Zero Age Main Sequence stellar mass in solar masses")
    _age_prop = ('age', "Stellar ages in years and decimal logarithmic scale")
    _surfacegravity_prop = ('grav', "Surface gravity in CGS units and decimal logarithmic scale")
    _metallicity_prop = ('feh', "Stellar metallicity [Fe/H] in dex relative to solar")
    _temperature_prop = ('teff', "Surface temperature in Kelvin and decimal logarithmic scale")
    _luminosity_prop = ('lum', "Stellar luminosity in solar luminosities and decimal logarithmic scale")
    _parentindex_prop = Input._parentindex_prop
    _particleflag_prop = ('partid', "Flag = 1 if star not at center of its parent particle")
    _export_properties =  {_position_prop, _velocity_prop, _celestial_prop, _galactic_prop, _distance_prop, _modulus_prop, _trgbmass_prop, _currentmass_prop, _zamsmass_prop, _age_prop, _surfacegravity_prop, _metallicity_prop, _temperature_prop, _luminosity_prop, _parentindex_prop, _particleflag_prop}
    #####
    _parallax_prop = ('pi', "Parallax in milliarcseconds")
    _propermotion_prop = (('mura', 'mudec'), "Equatorial proper motions in milliarcseconds per year")
    _galacticpropermotion_prop = (('mul', 'mub'), "Galactic proper motions in milliarcseconds per year")
    _radialvelocity_prop = ('vr', "Radial velocity in km/s")
    _postprocess_properties = {_parallax_prop, _propermotion_prop, _galacticpropermotion_prop, _radialvelocity_prop}
    #####
    _all_optional_properties = Input._optional_properties - {_parentindex_prop}
    #####
    _export_keys = tuple(sum([(_p[0],) if isinstance(_p[0], str) else _p[0] for _p in _export_properties], ()))
    _postprocess_keys = tuple(sum([(_p[0],) if isinstance(_p[0], str) else _p[0] for _p in _postprocess_properties], ()))
    #####
    _pos = _position_prop[0]  # ['px', 'py', 'pz']  # positions
    _vel = _velocity_prop[0]  # ['vx', 'vy', 'vz']  # velocities
    _cel = _celestial_prop[0]  # ['ra', 'dec']  # celestial coordinates
    _gal = _galactic_prop[0]  # ['glon', 'glat']  # galactic coordinates
    _rad = _distance_prop[0]  # 'rad'  # distance
    _dmod = _modulus_prop[0]  # 'dmod'  # distance modulus
    _mtip = _trgbmass_prop[0]  # 'mtip'  # mass at tip of giant branch for star of given age & metallicity (Msun)
    _mact = _currentmass_prop[0]  # 'mact'  # current stellar mass (Msun)
    _mini = _zamsmass_prop[0]  # 'smass'  # ZAMS stellar mass (Msun)
    _age = _age_prop[0]  # 'age'  # age (Gyr)
    _grav = _surfacegravity_prop[0]  # 'grav'  # surface gravity (log g)
    _feh = _metallicity_prop[0]  # 'feh'  # [fe/h]
    _teff = _surfacegravity_prop[0]  # 'teff'  # temperature (returned log10(teff/K) by Galaxia)
    _lum = _luminosity_prop[0]  # 'lum'  # luminosity (returned log10(lum/lsun) by Galaxia)
    _parentid = _parentindex_prop[0]  # 'parentid'  # id of parent particle 
    _partid = _particleflag_prop[0]  # 'partid'  # flag of particle's central star
    #####
    _pi = _parallax_prop[0]  # 'pi'  # parallax
    _mu = _propermotion_prop[0]  # ('mura', 'mudec')  # proper motions
    _mugal = _galacticpropermotion_prop[0]  # ('mul', 'mub')  # galactic proper motions
    _vr = _radialvelocity_prop[0]  # 'vr'  # radial velocity
    #####
    _vaex_under_list = ['_repr_html_']
    def __init__(self, survey: Survey, parameters: dict) -> None:
        """
            Driver to exploit the output of Galaxia.

            Call signature::

                output = Output(survey, parameters)

            Parameters
            ----------
            survey : :obj:`Survey`
                Survey object that returned this output.
            
            parameters : dict
                Dictionary all of parameters passed by Survey that were used
                to generate this output.
            
            Notes
            -----
            An Output object almost behaves as a vaex DataFrame, also please
            consult vaex online tutorials for more hands-on information:
                https://vaex.io/docs/tutorial.html

            The DataFrame represents the catalogue with columns corresponding
            to properties of the synthetic stars. Those include the photometric
            magnitudes per filter, with each filter identified by a key in the
            lowercase format "photosys_filtername" where photo_sys corresponds
            to the photometric system and filtername corresponds to a filter
            name of that system. With those are also always included the
            following properties: {_output_properties}
            
            Additionally, depending on what optional properties were provided
            with the input particle data, the output can also include the
            following properties: {_optional_properties}
        """
        self.__survey = survey
        self.__parameters = parameters
        self.__vaex = None
        self.__path = None

    __init__.__doc__ = __init__.__doc__.format(_output_properties=''.join(
                                                   [f"\n                 -{desc} via key `{str(key)}`"
                                                     for key, desc in _export_properties.union(_postprocess_properties)]),
                                               _optional_properties=''.join(
                                                   [f"\n                 -{desc} via key `{str(key)}`"
                                                     for key, desc in _all_optional_properties]))
    
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
        return tuple(itertools.chain.from_iterable([isochrone.to_export_keys for isochrone in isochrones]))
    
    @classmethod
    def _make_export_keys(cls, isochrones, extra_keys=()):
        return tuple(set(cls._export_keys).union(extra_keys).union(cls._compile_export_mag_names(isochrones)))

    @classmethod
    def _make_catalogue_keys(cls, isochrones, extra_keys=()):
        return cls._make_export_keys(isochrones, extra_keys=cls._postprocess_keys+extra_keys)

    def _make_input_optional_keys(self):
        return tuple(k if k != 'id' else 'satid' for k in self.survey.input.optional_keys())

    def _ebf_to_hdf5(self):
        hdf5_file = self._hdf5
        with h5.File(hdf5_file, 'w') as f5:
            for k in self.export_keys:
                # print(f"Exporting {k}...")
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
        self.flush_extra_columns_to_hdf5(with_columns=(self._teff, self._lum))

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
        self.flush_extra_columns_to_hdf5(with_columns=self._gal+(self._rad,))

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
    
    def flush_extra_columns_to_hdf5(self, with_columns=()):  # temporary until vaex supports it
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
