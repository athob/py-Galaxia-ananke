#!/usr/bin/env python
"""
Package templates
"""
from string import Template
from dataclasses import dataclass

from .utils import Singleton
from .constants import *


__all__ = ['TTAGS', 'FILENAME_TEMPLATE', 'PARFILE_TEMPLATE']

@dataclass(frozen=True)
class TemplateTags(metaclass=Singleton):
    name: str            = 'name'
    pname: str           = 'pname'
    output_file: str     = 'output_file'
    output_dir: str      = 'output_dir'
    photo_categ: str     = 'photo_categ'
    photo_sys: str       = 'photo_sys'
    mag_color_names: str = 'mag_color_names'
    app_mag_lim_lo: str  = 'app_mag_lim_lo'
    app_mag_lim_hi: str  = 'app_mag_lim_hi'
    abs_mag_lim_lo: str  = 'abs_mag_lim_lo'
    abs_mag_lim_hi: str  = 'abs_mag_lim_hi'
    color_lim_lo: str    = 'color_lim_lo'
    color_lim_hi: str    = 'color_lim_hi'
    geometry_opt: str    = 'geometry_opt'
    survey_area: str     = 'survey_area'
    fsample: str         = 'fsample'
    pop_id: str          = 'pop_id'
    warp_flare_on: str   = 'warp_flare_on'
    longitude: str       = 'longitude'
    latitude: str        = 'latitude'
    star_type: str       = 'star_type'
    photo_error: str     = 'photo_error'
    rand_seed: str       = 'rand_seed'
    r_max: str           = 'r_max'
    r_min: str           = 'r_min'
    nres: str            = 'nres'
    nstart: str          = 'nstart'
    rSun0: str           = 'rSun0'
    rSun1: str           = 'rSun1'
    rSun2: str           = 'rSun2'
    vSun0: str           = 'vSun0'
    vSun1: str           = 'vSun1'
    vSun2: str           = 'vSun2'
    
    @property
    def rSun(self):
        return [self.rSun0, self.rSun1, self.rSun2]
    
    @property
    def vSun(self):
        return [self.vSun0, self.vSun1, self.vSun2]

TTAGS = TemplateTags()

FILENAME_TEMPLATE = Template(NBODY1+f"/${{{TTAGS.name}}}/\n\t1\t1\n${{{TTAGS.pname}}}\n")  # TODO Template can't work for N>1 files
PARFILE_TEMPLATE = Template(f"""outputFile\t${{{TTAGS.output_file}}}\t#don't fiddle
outputDir\t${{{TTAGS.output_dir}}}\t#where to output the survey
photoCateg\t${{{TTAGS.photo_categ}}}\t#name of folder where to select magnitude system
photoSys\t${{{TTAGS.photo_sys}}}\t#magnitude system (see ananke-for-wings/GalaxiaData/Isochrones/padova/ for options)
magcolorNames\t${{{TTAGS.mag_color_names}}}\t#magnitude and color to use for selecting the CMD box
appMagLimits[0]\t${{{TTAGS.app_mag_lim_lo}}}\t#upper and lower limits in apparent mag
appMagLimits[1]\t${{{TTAGS.app_mag_lim_hi}}}
absMagLimits[0]\t${{{TTAGS.abs_mag_lim_lo}}}\t#upper and lower limits in absolute mag
absMagLimits[1]\t${{{TTAGS.abs_mag_lim_hi}}}
colorLimits[0]\t${{{TTAGS.color_lim_lo}}}\t#upper and lower limits in color defined on line 4
colorLimits[1]\t${{{TTAGS.color_lim_hi}}}
geometryOption\t${{{TTAGS.geometry_opt}}}\t#don't fiddle
surveyArea\t${{{TTAGS.survey_area}}}\t#not used
fSample\t${{{TTAGS.fsample}}}\t#don't fiddle
popID\t${{{TTAGS.pop_id}}}\t#don't fiddle
warpFlareOn\t${{{TTAGS.warp_flare_on}}}\t#not used
longitude\t${{{TTAGS.longitude}}}\t#not used
latitude\t${{{TTAGS.latitude}}}\t#not used
starType\t${{{TTAGS.star_type}}}\t#don't fiddle
photoError\t${{{TTAGS.photo_error}}}\t#not used
seed\t${{{TTAGS.rand_seed}}}\t#change if you want a different random sample
r_max\t${{{TTAGS.r_max}}}\t#max distance from galactic center to include
r_min\t${{{TTAGS.r_min}}}\t#min distance from galactic center to include
nres\t${{{TTAGS.nres}}}\t#nres
nstart\t${{{TTAGS.nstart}}}\t#integer at which to start numbering synthetic stars
rSun[0]\t${{{TTAGS.rSun0}}}\t#location of survey viewpoint relative to galactic center
rSun[1]\t${{{TTAGS.rSun1}}}
rSun[2]\t${{{TTAGS.rSun2}}}
vSun[0]\t${{{TTAGS.vSun0}}}\t#velocity of survey viewpoint relative to galactic center (not used)
vSun[1]\t${{{TTAGS.vSun1}}}
vSun[2]\t${{{TTAGS.vSun2}}}
""")
