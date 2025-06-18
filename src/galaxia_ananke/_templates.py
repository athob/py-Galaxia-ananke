#!/usr/bin/env python
"""
Package templates
"""
from typing import List
from string import Template
from dataclasses import dataclass

from .utils import Singleton
from ._constants import *


__all__ = ['FTTAGS', 'FILENAME_TEMPLATE', 'PARFILENAME_TEMPLATE', 'PARFILE_TEMPLATE', 'CTTAGS', 'RUN_TEMPLATE', 'HDIMBLOCK_TEMPLATE', 'APPEND_TEMPLATE']

@dataclass(frozen=True)
class FileTemplateTags(metaclass=Singleton):
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
    def rSun(self) -> List[str]:
        return [self.rSun0, self.rSun1, self.rSun2]
    
    @property
    def vSun(self) -> List[str]:
        return [self.vSun0, self.vSun1, self.vSun2]
    
    def append_photo_categ(self, n: int) -> str:
        return f"{self.photo_categ}_{n}"
    
    def append_photo_sys(self, n: int) -> str:
        return f"{self.photo_sys}_{n}"
    
    def append_photo(self, n: int) -> List[str]:
        return [self.append_photo_categ(n), self.append_photo_sys(n)]


FTTAGS = FileTemplateTags()

FILENAME_TEMPLATE = Template(NBODY1+f"/${{{FTTAGS.name}}}/\n\t1\t1\n${{{FTTAGS.pname}}}\n")  # TODO Template can't work for N>1 files
PARFILENAME_TEMPLATE = Template(f"${{{FTTAGS.name}}}_params")
PARFILE_TEMPLATE = Template(f"""outputFile\t${{{FTTAGS.output_file}}}\t#don't fiddle
outputDir\t${{{FTTAGS.output_dir}}}\t#where to output the survey
photoCateg\t${{{FTTAGS.photo_categ}}}\t#name of folder where to select magnitude system
photoSys\t${{{FTTAGS.photo_sys}}}\t#magnitude system (see ananke-for-wings/GalaxiaData/Isochrones/padova/ for options)
magcolorNames\t${{{FTTAGS.mag_color_names}}}\t#magnitude and color to use for selecting the CMD box
appMagLimits[0]\t${{{FTTAGS.app_mag_lim_lo}}}\t#upper and lower limits in apparent mag
appMagLimits[1]\t${{{FTTAGS.app_mag_lim_hi}}}
absMagLimits[0]\t${{{FTTAGS.abs_mag_lim_lo}}}\t#upper and lower limits in absolute mag
absMagLimits[1]\t${{{FTTAGS.abs_mag_lim_hi}}}
colorLimits[0]\t${{{FTTAGS.color_lim_lo}}}\t#upper and lower limits in color defined on line 4
colorLimits[1]\t${{{FTTAGS.color_lim_hi}}}
geometryOption\t${{{FTTAGS.geometry_opt}}}\t#don't fiddle
surveyArea\t${{{FTTAGS.survey_area}}}\t#not used
fSample\t${{{FTTAGS.fsample}}}\t#don't fiddle
popID\t${{{FTTAGS.pop_id}}}\t#don't fiddle
warpFlareOn\t${{{FTTAGS.warp_flare_on}}}\t#not used
longitude\t${{{FTTAGS.longitude}}}\t#not used
latitude\t${{{FTTAGS.latitude}}}\t#not used
starType\t${{{FTTAGS.star_type}}}\t#don't fiddle
photoError\t${{{FTTAGS.photo_error}}}\t#not used
seed\t${{{FTTAGS.rand_seed}}}\t#change if you want a different random sample
r_max\t${{{FTTAGS.r_max}}}\t#max distance from galactic center to include
r_min\t${{{FTTAGS.r_min}}}\t#min distance from galactic center to include
nres\t${{{FTTAGS.nres}}}\t#nres
nstart\t${{{FTTAGS.nstart}}}\t#integer at which to start numbering synthetic stars
rSun[0]\t${{{FTTAGS.rSun0}}}\t#location of survey viewpoint relative to galactic center
rSun[1]\t${{{FTTAGS.rSun1}}}
rSun[2]\t${{{FTTAGS.rSun2}}}
vSun[0]\t${{{FTTAGS.vSun0}}}\t#velocity of survey viewpoint relative to galactic center (not used)
vSun[1]\t${{{FTTAGS.vSun1}}}
vSun[2]\t${{{FTTAGS.vSun2}}}
""")

@dataclass(frozen=True)
class CommandTemplateTags(metaclass=Singleton):
    hdim_block: str = 'hdim_block'
    hdim: str       = 'hdim'
    nfile: str      = 'nfile'
    ngen: str       = 'ngen'
    parfile: str    = 'parfile'
    pcat: str       = 'pcat'
    psys: str       = 'psys'
    filename: str   = 'filename'

CTTAGS = CommandTemplateTags()

RUN_TEMPLATE = Template(f"{GALAXIA} -r${{{CTTAGS.hdim_block}}} --nfile=${{{CTTAGS.nfile}}} --ngen=${{{CTTAGS.ngen}}} ${{{CTTAGS.parfile}}}")
HDIMBLOCK_TEMPLATE = Template(f" --hdim=${{{CTTAGS.hdim}}}")
APPEND_TEMPLATE = Template(f"{GALAXIA} -a --pcat=${{{CTTAGS.pcat}}} --psys=${{{CTTAGS.psys}}} ${{{CTTAGS.filename}}}")


if __name__ == '__main__':
    raise NotImplementedError()
