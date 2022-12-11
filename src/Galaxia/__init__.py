#!/usr/bin/env python
"""
Docstring
"""
import re
import pathlib
import subprocess
import itertools
import numpy as np
import h5py as h5
import ebf

from .constants import *
from . import photometry

__all__ = []


def write_kernel_for_galaxia(masses, rho_pos, rho_vel, name, knorm, ngb=64):
    kernels = np.vstack([np.sqrt(ngb) * knorm / np.cbrt(rho_pos),
                         np.sqrt(ngb) * knorm / np.cbrt(rho_vel)]).T  # TODO what if rho_vel is None
    kname = GALAXIA_TMP / f"{name}_d{[6, 3][rho_vel is None]}n{ngb}_den.ebf"
    ebf.initialize(kname)
    ebf.write(kname, '/density', rho_pos, "a")
    ebf.write(kname, '/h_cubic', kernels, "a")
    ebf.write(kname, '/mass', masses, "a")
    return pathlib.Path(kname)


def write_particles_for_galaxia(particles, name):
    # TODO check here if particles has all required keys, add a list of required keys to constants.py
    pname = GALAXIA_TMP / f"{name}.ebf"
    ebf.initialize(pname)
    for key in particles.keys():
        ebf.write(pname, '/' + key, particles[key], 'a')
    # and a placeholder field (for multi-component systems this can be used to e.g. keep track of which component)
    ebf.write(pname, '/id', np.zeros(particles[key].shape[0]), 'a')
    return pathlib.Path(pname)


def make_symlink(file_path, dest_dir):
    file_path = pathlib.Path(file_path).resolve()
    symlink_name = dest_dir / file_path.name
    try:
        symlink_name.unlink()
    except FileNotFoundError:
        pass
    symlink_name.symlink_to(file_path)


def prepare_galaxia_input(particles, rho_pos, rho_vel, simname, isochrone, cmd_magnames, ngb, **kwargs):
    check = set(re.split('[ ,-]', cmd_magnames))
    if not check.issubset(set(isochrone.mag_names)):
        raise ValueError(f"CMD magnitude names '{cmd_magnames}' don't match isochrone names: choose among {isochrone.mag_names}")
    masses = particles['mass']  # TODO not really a good design
    knorm = kwargs.pop('knorm', 0.596831)
    kname = write_kernel_for_galaxia(masses, rho_pos, rho_vel, simname, knorm, ngb)
    pname = write_particles_for_galaxia(particles, simname)
    temp_dir = GALAXIA_NBODY1 / simname
    temp_dir.mkdir(parents=True, exist_ok=True)
    make_symlink(kname, temp_dir)
    make_symlink(pname, temp_dir)
    temp_filename = (GALAXIA_FILENAMES / simname).with_suffix('.txt')
    temp_filename.write_text(FILENAME_TEMPLATE.substitute(name=simname, pname=pname.name))
    parfile = pathlib.Path(kwargs.pop('parfile', 'survey_params'))  # TODO make temporary? create a global record of temporary files?
    if not parfile.is_absolute():
        parfile = GALAXIA_TMP / parfile
    parfile.write_text(PARFILE_TEMPLATE.substitute(DEFAULTS_FOR_PARFILE, photo_categ=isochrone.category, photo_sys=isochrone.name, mag_color_names=cmd_magnames, nres=ngb, **kwargs))  # output_file=surveyname, fsample=fsample))
    return simname, parfile


def run_galaxia_on_input(simname, parfile, surveyname, hdim=None):
    surveyname_base = GALAXIA_TMP / f"{surveyname}.{simname}"
    ebf_output_file = surveyname_base.parent / (surveyname_base.name + '.ebf')
    cmd = f"{GALAXIA} -r{(' --hdim=' + str(hdim) if hdim is not None else '')} --nfile={simname} {parfile}"
    print(cmd)
    subprocess.call(cmd.split(' '))
    return ebf_output_file


def append_galaxia_output(ebf_output_file, isochrone):
    cmd = f"{GALAXIA} -a --pcat={isochrone.category} --psys={isochrone.name} {ebf_output_file}"
    print(cmd)
    subprocess.call(cmd.split(' '))


def convert_ebf_to_hdf5(ebf_file, isochrones):
    export_keys = list(itertools.chain.from_iterable([EXPORT_KEYS]+[isochrone.to_export_keys for isochrone in isochrones]))
    hdf5_file = ebf_file.with_suffix('.hdf5')
    data = ebf.read(ebf_file)
    with h5.File(hdf5_file, 'w') as f5:
        for k in export_keys:
            f5.create_dataset(name=k, data=data[k])
        print(f"Exported the following quantities to {hdf5_file}")
        print(list(f5.keys()))
    return hdf5_file


def make_survey_from_particles(particles, rho_pos, rho_vel, photo_sys=DEFAULT_PSYS, cmd_magnames=DEFAULT_CMD, simname='sim', surveyname='survey', fsample=1, hdim=None, ngb=64, knorm=0.596831, **kwargs):
    if isinstance(photo_sys, str):
        photo_sys = [photo_sys]
    isochrones = [photometry.available_photo_systems[psys] for psys in photo_sys]
    simname, parfile = prepare_galaxia_input(particles, rho_pos, rho_vel, simname, isochrones[0], cmd_magnames, ngb, knorm=knorm, output_file=surveyname, fsample=fsample, **kwargs)
    ebf_output_file = run_galaxia_on_input(simname, parfile, surveyname, hdim=hdim)
    for isochrone in isochrones[1:]:
        append_galaxia_output(ebf_output_file, isochrone)
    hdf5_output_file = convert_ebf_to_hdf5(ebf_output_file, isochrones)
    return hdf5_output_file


if __name__ == '__main__':
    pass
