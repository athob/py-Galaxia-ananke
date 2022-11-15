#!/usr/bin/env python
"""
Package parameters
"""
import pathlib
from string import Template

__all__ = ['NAME', 'GALAXIA_SUBMODULE_NAME', 'GALAXIA_URL', 'GALAXIA_EXEC', 'SRC_DIR', 'LOG_DIR', 'GALAXIA_DATA', 'GALAXIA_NBODY1', 'GALAXIA_FILENAMES', 'GALAXIA_LOG', 'GALAXIA_TMP', 'GALAXIA', 'CACHE', 'EXPORT_KEYS', 'FILENAME_TEMPLATE', 'PARFILE_TEMPLATE', 'DEFAULTS_FOR_PARFILE']

NAME = 'Galaxia'
GALAXIA_SUBMODULE_NAME = 'galaxia'
GALAXIA_URL = 'https://sourceforge.net/projects/galaxia/files/galaxia-0.7.2.tar.gz/download'
GALAXIA_EXEC = 'galaxia'
SRC_DIR = 'src'
BIN_DIR = 'bin'
LOG_DIR = 'log'
TMP_DIR = 'tmp'
GALDATA_DIR = 'GalaxiaData'

NBODY1 = 'nbody1'
FILENAMES = 'filenames'

GLOBAL_CACHE = '~/.cache'  # TODO must adapt to OS

CACHE = pathlib.Path(GLOBAL_CACHE).expanduser().resolve() / NAME
GALAXIA_DATA = CACHE / GALDATA_DIR
GALAXIA_NBODY1 = GALAXIA_DATA / NBODY1
GALAXIA_FILENAMES = GALAXIA_NBODY1 / FILENAMES
GALAXIA_LOG = CACHE / LOG_DIR
GALAXIA_TMP = pathlib.Path.cwd()  # CACHE / TMP_DIR   # TODO use temporary directory
GALAXIA = CACHE / BIN_DIR / GALAXIA_EXEC

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
    'photo_categ': 'padova',
    'photo_sys': 'WFIRST-HST',  # TODO should it have a default value?
    'mag_color_names': 'F814W,F555W-F814W',  # TODO as above?
    'app_mag_lim_lo': -1000,
    'app_mag_lim_hi': 1000,
    'abs_mag_lim_lo': -7.0,
    'abs_mag_lim_hi': -2.5,
    'color_lim_lo': -1000,
    'color_lim_hi': 1000,
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
EXPORT_KEYS = ['rad', 'wfirst_r062', 'teff', 'alpha', 'wfirst-hst_f606w', 'wfirst-hst_w149', 'wfirst_h158', 'mtip',
               'pz', 'px', 'py', 'wfirst-hst_f160w', 'feh', 'wfirst-hst_y106', 'wfirst_j129', 'wfirst-hst_f110w',
               'wfirst-hst_h158', 'wfirst-hst_f184', 'lum', 'wfirst_w149', 'mact', 'wfirst-hst_f814w',
               'wfirst-hst_f475w', 'parentid', 'dmod', 'wfirst-hst_z087', 'partid', 'wfirst-hst_f555w', 'wfirst_y106',
               'wfirst_f184', 'age', 'grav', 'wfirst-hst_j129', 'wfirst_z087', 'smass']  # TODO this differs with photo_sys and mag_color_names
