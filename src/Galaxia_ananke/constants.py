#!/usr/bin/env python
"""
Package parameters
"""
import sys
import pathlib
from string import Template
from dataclasses import dataclass

from .utils import Singleton

__all__ = ['NAME', 'GALAXIA_SUBMODULE_NAME', 'GALAXIA_URL', 'GALAXIA_EXEC', 'SRC_DIR', 'LOG_DIR', 'LEGACY_PHOTOCAT', 'CUSTOM_PHOTOCAT', 'DEFAULT_PSYS', 'DEFAULT_CMD', 'DEFAULT_CMD_BOX', 'DEFAULT_SIMNAME', 'DEFAULT_SURVEYNAME', 'DEFAULT_PARFILE', 'GALAXIA_DATA', 'GALAXIA_NBODY1', 'GALAXIA_FILENAMES', 'GALAXIA_LOG', 'GALAXIA_TMP', 'GALAXIA', 'ISOCHRONES_PATH', 'CACHE', 'TTAGS', 'FILENAME_TEMPLATE', 'PARFILE_TEMPLATE', 'DEFAULTS_FOR_PARFILE']

NAME = 'Galaxia_ananke'
GALAXIA_SUBMODULE_NAME = 'galaxia-ananke'
GALAXIA_URL = 'https://sourceforge.net/projects/galaxia/files/galaxia-0.7.2.tar.gz/download'
GALAXIA_EXEC = 'galaxia'
SRC_DIR = 'src'
BIN_DIR = 'bin'
LOG_DIR = 'log'
TMP_DIR = 'tmp'
GALDATA_DIR = 'GalaxiaData'

NBODY1 = 'nbody1'
FILENAMES = 'filenames'
LEGACY_PHOTOCAT = 'padova'
CUSTOM_PHOTOCAT = 'py_custom'

DEFAULT_PSYS = ['padova/GAIADR2']
DEFAULT_CMD = 'Gmag,G_BPmag-G_RPmag'
# DEFAULT_CMD = {'magnitude': 'Gmag', 'color_minuend': 'G_BPmag', 'color_subtrahend': 'G_RPmag'}
DEFAULT_CMD_BOX = {'app_mag': [-1000,1000], 'abs_mag': [-1000,20], 'color': [-1000,1000]}

DEFAULT_SIMNAME = 'sim'
DEFAULT_SURVEYNAME = 'survey'
DEFAULT_PARFILE = 'survey_params'

PREFIX_ENV_VAR = "ANANKE_SYSTEM_PREFIX"
if PREFIX_ENV_VAR in pathlib.os.environ:
    PREFIX = pathlib.Path(pathlib.os.environ[PREFIX_ENV_VAR])
else:
    PREFIX = pathlib.Path(sys.prefix)
if not pathlib.os.access(PREFIX, pathlib.os.W_OK):
    PermissionError(f"Installation cannot complete: to proceed, please give write permission to directory {PREFIX} or define a custom system prefix via the environment variable {PREFIX_ENV_VAR}")  # TODO make the variable conda-dependent?

GLOBAL_CACHE = PREFIX / '.cache'

CACHE = pathlib.Path(GLOBAL_CACHE).expanduser().resolve() / NAME
GALAXIA_DATA = CACHE / GALDATA_DIR
GALAXIA_NBODY1 = GALAXIA_DATA / NBODY1
GALAXIA_FILENAMES = GALAXIA_NBODY1 / FILENAMES
GALAXIA_LOG = CACHE / LOG_DIR
GALAXIA_TMP = pathlib.Path.cwd()  # CACHE / TMP_DIR   # TODO use temporary directory
GALAXIA = CACHE / BIN_DIR / GALAXIA_EXEC

ISOCHRONES_PATH = GALAXIA_DATA / 'Isochrones'

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
    TTAGS.rSun0: 0.0,
    TTAGS.rSun1: 0.0,
    TTAGS.rSun2: 0.0,
    TTAGS.vSun0: 0.0,
    TTAGS.vSun1: 0.0,
    TTAGS.vSun2: 0.0}
