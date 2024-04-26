#!/usr/bin/env python
"""
Package defaults
"""
import pathlib
from dataclasses import dataclass

from astropy import coordinates, units

from .utils import Singleton
from .templates import *


__all__ = ['GALAXIA_TMP', 'DEFAULT_PSYS', 'DEFAULT_CMD', 'DEFAULT_CMD_BOX', 'DEFAULT_SIMNAME', 'DEFAULT_SURVEYNAME', 'DEFAULT_PARFILE', 'DEF_UNIT', 'DEFAULTS_FOR_PARFILE']

DEFAULT_PSYS = ['padova/GAIA__DR2']
DEFAULT_CMD = 'G,Gbp-Grp'
# DEFAULT_CMD = {'magnitude': 'G', 'color_minuend': 'Gbp', 'color_subtrahend': 'Grp'}
DEFAULT_CMD_BOX = {'app_mag': [-1000,1000], 'abs_mag': [-1000,20], 'color': [-1000,1000]}

DEFAULT_SIMNAME = 'sim'
DEFAULT_SURVEYNAME = 'survey'
DEFAULT_PARFILE = 'survey_params'

GALAXIA_TMP = pathlib.Path.cwd()  # CACHE / TMP_DIR   # TODO use temporary directory

@dataclass(frozen=True)
class DefaultUnits(metaclass=Singleton):
    position: units.Unit   = units.kpc
    velocity: units.Unit   = units.km/units.s
    wavelength: units.Unit = units.micron
    irradiance: units.Unit = units.Unit('erg/(cm2 s)')
    spectral: units.Unit   = irradiance/wavelength

DEF_UNIT = DefaultUnits()

heliocentric_center = coordinates.SkyCoord(x=0*DEF_UNIT.position, y=0*DEF_UNIT.position, z=0*DEF_UNIT.position, frame='hcrs', representation_type='cartesian')
rSun = heliocentric_center.galactocentric.cartesian.xyz.to(DEF_UNIT.position).value
vSun = heliocentric_center.galactocentric.frame.galcen_v_sun.get_d_xyz().to(DEF_UNIT.velocity).value

DEFAULTS_FOR_PARFILE = {
    # TTAGS.output_file: ,  # TODO use temporary file
    TTAGS.output_dir: GALAXIA_TMP,
    TTAGS.app_mag_lim_lo: DEFAULT_CMD_BOX['app_mag'][0],
    TTAGS.app_mag_lim_hi: DEFAULT_CMD_BOX['app_mag'][1],
    TTAGS.abs_mag_lim_lo: DEFAULT_CMD_BOX['abs_mag'][0],
    TTAGS.abs_mag_lim_hi: DEFAULT_CMD_BOX['abs_mag'][1],
    TTAGS.color_lim_lo: DEFAULT_CMD_BOX['color'][0],
    TTAGS.color_lim_hi: DEFAULT_CMD_BOX['color'][1],
    TTAGS.geometry_opt: 0,  # shouldn't use?
    TTAGS.survey_area: 207.455,
    TTAGS.fsample: 1,
    TTAGS.pop_id: 10,
    TTAGS.warp_flare_on: 0,
    TTAGS.longitude: 76.2730,
    TTAGS.latitude: 13.4725,
    TTAGS.star_type: 0,
    TTAGS.photo_error: 0,
    TTAGS.rand_seed: 17052,  # TODO randomize default?
    TTAGS.r_max: 500,
    TTAGS.r_min: 0,
    TTAGS.nres: 64,
    TTAGS.nstart: 0,
    TTAGS.rSun0: rSun[0],
    TTAGS.rSun1: rSun[1],
    TTAGS.rSun2: rSun[2],
    TTAGS.vSun0: vSun[0],
    TTAGS.vSun1: vSun[1],
    TTAGS.vSun2: vSun[2]}
