#!/usr/bin/env python
"""
Package parameters
"""
import pathlib
from string import Template

__all__ = ['NAME', 'GALAXIA_SUBMODULE_NAME', 'GALAXIA_URL', 'GALAXIA_EXEC', 'SRC_DIR', 'LOG_DIR', 'GALAXIA_DATA', 'GALAXIA_NBODY1', 'GALAXIA_FILENAMES', 'GALAXIA_LOG', 'GALAXIA', 'CACHE', 'EXPORT_KEYS', 'FILENAME_TEMPLATE', 'PARFILE_TEMPLATE', 'DEFAULTS_FOR_PARFILE']

NAME = 'Galaxia'
GALAXIA_SUBMODULE_NAME = 'galaxia'
GALAXIA_URL = 'https://sourceforge.net/projects/galaxia/files/galaxia-0.7.2.tar.gz/download'
GALAXIA_EXEC = 'galaxia'
SRC_DIR = 'src'
BIN_DIR = 'bin'
LOG_DIR = 'log'
GALDATA_DIR = 'GalaxiaData'

NBODY1 = 'nbody1'
FILENAMES = 'filenames'
FILENAME_TEMPLATE = Template(NBODY1+"/${name}/\n\t1\t1\n${pname}\n")  # TODO Template can't work for N>1 files
PARFILE_TEMPLATE = Template("""outputFile\t${surveyname}\t#don't fiddle
outputDir\t./\t#where to output the survey
photoSys\tWFIRST-HST\t#magnitude system (see ananke-for-wings/GalaxiaData/Isochrones/padova/ for options)
magcolorNames\tF814W,F555W-F814W\t#magnitude and color to use for selecting the CMD box
appMagLimits[0]\t-1000\t#upper and lower limits in apparent mag
appMagLimits[1]\t1000
absMagLimits[0]\t-7.0\t#upper and lower limits in absolute mag
absMagLimits[1]\t2.5
colorLimits[0]\t-1000\t#upper and lower limits in color defined on line 4
colorLimits[1]\t1000
geometryOption\t0\t#don't fiddle
surveyArea\t207.455\t#not used
fSample\t${fsample}\t#don't fiddle
popID\t10\t#don't fiddle
warpFlareOn\t0\t#not used
longitude\t76.2730\t#not used
latitude\t13.4725\t#not used
starType\t0\t#don't fiddle
photoError\t0\t#not used
seed\t17052\t#change if you want a different random sample
r_max\t500\t#max distance from galactic center to include
r_min\t0\t#min distance from galactic center to include
nres\t${nres}\t#nres
nstart\t0\t#integer at which to start numbering synthetic stars
rSun[0]\t${rSun0}\t#location of survey viewpoint relative to galactic center
rSun[1]\t${rSun1}
rSun[2]\t${rSun2}
vSun[0]\t${vSun0}\t#velocity of survey viewpoint relative to galactic center (not used)
vSun[1]\t${vSun1}
vSun[2]\t${vSun2}
""")
DEFAULTS_FOR_PARFILE = {
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
               'wfirst_f184', 'age', 'grav', 'wfirst-hst_j129', 'wfirst_z087', 'smass']

GLOBAL_CACHE = '~/.cache'  # TODO must adapt to OS

CACHE = pathlib.Path(GLOBAL_CACHE).expanduser().resolve() / NAME
GALAXIA_DATA = CACHE / GALDATA_DIR
GALAXIA_NBODY1 = GALAXIA_DATA / NBODY1
GALAXIA_FILENAMES = GALAXIA_NBODY1 / FILENAMES
GALAXIA_LOG = CACHE / LOG_DIR
GALAXIA = CACHE / BIN_DIR / GALAXIA_EXEC
