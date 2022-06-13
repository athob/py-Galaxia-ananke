#!/usr/bin/env python
"""
Docstring
"""
import pathlib
import subprocess
import numpy as np
import h5py as h5
import ebf

from .constants import *

__all__ = []


def write_kernel_for_galaxia(masses, rho_pos, rho_vel=None, name='.', knorm=0.596831, ngb=64):
    kernels = np.vstack([np.sqrt(ngb) * knorm / np.cbrt(rho_pos),
                         np.sqrt(ngb) * knorm / np.cbrt(rho_vel)]).T  # TODO what if rho_vel = None
    kname = f"{name}_d{[6, 3][rho_vel is None]}n{ngb}_den.ebf"
    ebf.initialize(kname)
    ebf.write(kname, '/density', rho_pos, "a")
    ebf.write(kname, '/h_cubic', kernels, "a")
    ebf.write(kname, '/mass', masses, "a")
    return pathlib.Path(kname)


def write_particles_for_galaxia(particles, name='.'):
    pname = f"{name}.ebf"
    ebf.initialize(pname)
    for key in particles.keys():
        ebf.write(pname, '/' + key, particles[key], 'a')
    # and a placeholder field (for multi-component systems this can be used to e.g. keep track of which component)
    ebf.write(pname, '/id', np.zeros(particles[key].shape[0]), 'a')
    return pathlib.Path(pname)


def make_symlink(file_path, dest_dir):
    file_path = pathlib.Path(file_path).resolve()
    symlink_pname = dest_dir / file_path.name
    try:
        symlink_pname.unlink()
    except FileNotFoundError:
        pass
    symlink_pname.symlink_to(file_path)


def prepare_galaxia_input_files(kname, pname):
    temp_name = 'temp'
    temp_dir = GALAXIA_NBODY1 / temp_name
    temp_dir.mkdir(parents=True, exist_ok=True)
    make_symlink(kname, temp_dir)
    make_symlink(pname, temp_dir)
    (GALAXIA_FILENAMES / temp_name).with_suffix('.txt').write_text(FILENAME_TEMPLATE.substitute(name=temp_name,
                                                                                                pname=pname.name))
    return temp_name


def ebf_to_hdf5(name):
    # list of keys to export to hdf5 - don't fiddle with these unless you know what you're doing
    data = ebf.read(f"{name}.ebf")
    f = h5.File(f"{name}.hdf5", 'w')
    for k in EXPORT_KEYS:
        f.create_dataset(name=k, data=data[k])
    print(f"Exported the following quantities to {name}.hdf5")
    print(list(f.keys()))
    f.close()


def make_survey(simname, surveyname, fsample=1, hdim=None, ngb=64, photo_sys='WFIRST', parfile='survey_params'):
    with open(parfile, 'w') as f:
        f.write(PARFILE_TEMPLATE.substitute(DEFAULTS_FOR_PARFILE,
                                            nres=ngb,
                                            surveyname=surveyname, fsample=fsample))
    cmd = f"{GALAXIA} -r{(' --hdim=' + str(hdim) if hdim is not None else '')} --nfile={simname} {parfile}"
    print(cmd)
    subprocess.call(cmd.split(' '))
    surveyname_full = f"{surveyname}.{simname}.ebf"
    cmd = f"{GALAXIA} -a --psys={photo_sys} {surveyname_full}"
    print(cmd)
    subprocess.call(cmd.split(' '))
    surveyname_stem = f"{surveyname}.{simname}"
    ebf_to_hdf5(surveyname_stem)


if __name__ == '__main__':
    pass
