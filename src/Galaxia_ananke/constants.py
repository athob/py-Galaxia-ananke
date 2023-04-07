#!/usr/bin/env python
"""
Package parameters
"""
import sys
import pathlib
from string import Template

__all__ = ['NAME', 'GALAXIA_SUBMODULE_NAME', 'GALAXIA_URL', 'GALAXIA_EXEC', 'SRC_DIR', 'LOG_DIR', 'LEGACY_PHOTOCAT', 'CUSTOM_PHOTOCAT', 'DEFAULT_PSYS', 'DEFAULT_CMD', 'GALAXIA_DATA', 'GALAXIA_NBODY1', 'GALAXIA_FILENAMES', 'GALAXIA_LOG', 'GALAXIA_TMP', 'GALAXIA', 'ISOCHRONES_PATH', 'CACHE', 'FILENAME_TEMPLATE', 'PARFILE_TEMPLATE', 'DEFAULTS_FOR_PARFILE']

NAME = 'Galaxia_ananke'
GALAXIA_SUBMODULE_NAME = 'galaxia_ananke'
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

DEFAULT_PSYS = ['padova/WFIRST-HST']
DEFAULT_CMD = 'F814W,F555W-F814W'
# DEFAULT_CMD = {'magnitude': 'F814W', 'color_minuend': 'F555W', 'color_subtrahend': 'F814W'}
DEFAULT_CMD_BOX = {'app_mag': [-1000,1000], 'abs_mag': [-1000,5], 'color': [-1000,1000]}

GLOBAL_CACHE = pathlib.Path(sys.prefix) / '.cache'

CACHE = pathlib.Path(GLOBAL_CACHE).expanduser().resolve() / NAME
GALAXIA_DATA = CACHE / GALDATA_DIR
GALAXIA_NBODY1 = GALAXIA_DATA / NBODY1
GALAXIA_FILENAMES = GALAXIA_NBODY1 / FILENAMES
GALAXIA_LOG = CACHE / LOG_DIR
GALAXIA_TMP = pathlib.Path.cwd()  # CACHE / TMP_DIR   # TODO use temporary directory
GALAXIA = CACHE / BIN_DIR / GALAXIA_EXEC

ISOCHRONES_PATH = GALAXIA_DATA / 'Isochrones'

FILENAME_TEMPLATE = Template(NBODY1+"/${name}/\n\t1\t1\n${pname}\n")  # TODO Template can't work for N>1 files
PARFILE_TEMPLATE = Template("""outputFile\t${output_file}\t#don't fiddle
outputDir\t${output_dir}\t#where to output the survey
photoCateg\t${photo_categ}\t#name of folder where to select magnitude system
photoSys\t${photo_sys}\t#magnitude system (see ananke-for-wings/GalaxiaData/Isochrones/padova/ for options)
magcolorNames\t${mag_color_names}\t#magnitude and color to use for selecting the CMD box
appMagLimits[0]\t${app_mag_lim_lo}\t#upper and lower limits in apparent mag
appMagLimits[1]\t${app_mag_lim_hi}
absMagLimits[0]\t${abs_mag_lim_lo}\t#upper and lower limits in absolute mag
absMagLimits[1]\t${abs_mag_lim_hi}
colorLimits[0]\t${color_lim_lo}\t#upper and lower limits in color defined on line 4
colorLimits[1]\t${color_lim_hi}
geometryOption\t${geometry_opt}\t#don't fiddle
surveyArea\t${survey_area}\t#not used
fSample\t${fsample}\t#don't fiddle
popID\t${pop_id}\t#don't fiddle
warpFlareOn\t${warp_flare_on}\t#not used
longitude\t${longitude}\t#not used
latitude\t${latitude}\t#not used
starType\t${star_type}\t#don't fiddle
photoError\t${photo_error}\t#not used
seed\t${rand_seed}\t#change if you want a different random sample
r_max\t${r_max}\t#max distance from galactic center to include
r_min\t${r_min}\t#min distance from galactic center to include
nres\t${nres}\t#nres
nstart\t${nstart}\t#integer at which to start numbering synthetic stars
rSun[0]\t${rSun0}\t#location of survey viewpoint relative to galactic center
rSun[1]\t${rSun1}
rSun[2]\t${rSun2}
vSun[0]\t${vSun0}\t#velocity of survey viewpoint relative to galactic center (not used)
vSun[1]\t${vSun1}
vSun[2]\t${vSun2}
""")
DEFAULTS_FOR_PARFILE = {
    # 'output_file': ,  # TODO use temporary file
    'output_dir': GALAXIA_TMP,
    'app_mag_lim_lo': DEFAULT_CMD_BOX['app_mag'][0],
    'app_mag_lim_hi': DEFAULT_CMD_BOX['app_mag'][1],
    'abs_mag_lim_lo': DEFAULT_CMD_BOX['abs_mag'][0],
    'abs_mag_lim_hi': DEFAULT_CMD_BOX['abs_mag'][1],
    'color_lim_lo': DEFAULT_CMD_BOX['color'][0],
    'color_lim_hi': DEFAULT_CMD_BOX['color'][1],
    'geometry_opt': 0,  # shouldn't use?
    'survey_area': 207.455,
    'fsample': 1,
    'pop_id': 10,
    'warp_flare_on': 0,
    'longitude': 76.2730,
    'latitude': 13.4725,
    'star_type': 0,
    'photo_error': 0,
    'rand_seed': 17052,  # TODO randomize default?
    'r_max': 500,
    'r_min': 0,
    'nres': 64,
    'nstart': 0,
    'rSun0': 0.0,
    'rSun1': 0.0,
    'rSun2': 0.0,
    'vSun0': 0.0,
    'vSun1': 0.0,
    'vSun2': 0.0}
