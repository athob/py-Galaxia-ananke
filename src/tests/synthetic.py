#!/usr/bin/env python
#
# Author: Adrien CR Thob
# Copyright (C) 2026  Adrien CR Thob
#
# This file is part of the py-Galaxia-ananke project,
# <https://github.com/athob/py-Galaxia-ananke>, which is licensed
# under the GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later).
#
# The full copyright notice, including terms governing use, modification,
# and redistribution, is contained in the files LICENSE and COPYRIGHT,
# which can be found at the root of the source code distribution tree:
# - LICENSE <https://github.com/athob/py-Galaxia-ananke/blob/main/LICENSE>
# - COPYRIGHT <https://github.com/athob/py-Galaxia-ananke/blob/main/COPYRIGHT>
#
"""
Synthetic output of the backend, for the tests of Output: surveys that are prepared (so that they have an Output) without running
anything, and EBF files laid out as Galaxia writes them for such a survey: a dataset for each key that Output exports (floats in
single precision, the indices in 64 bits), the partition of each star, and the log of the run. What the stars are is simple on
purpose, so that what the tests expect can be worked out by hand.
"""
import types

import ebf
import numpy as np

from ..galaxia_ananke import Input, Survey, photometry
from ..galaxia_ananke._defaults import DEFAULT_CMD
from ..galaxia_ananke.photometry.PhotoSystem import PhotoSystem
from ..galaxia_ananke.photometry.Photometry import nested_dict
from .utils import patch_nbody1

INTEGER_KEYS = ('parentid', 'partid', 'partitionid', 'satid')


def photo_system(name, mag_names):
    return PhotoSystem(types.SimpleNamespace(category='padova', name=name, mag_names=list(mag_names)))


def use_photometry(monkeypatch):
    """The registry of the photometric systems is made of two systems that need nothing of the installation."""
    monkeypatch.setattr(photometry, 'available_photo_systems', nested_dict(
        {'padova': {'GAIA__DR2': photo_system('GAIA__DR2', ['G', 'Gbp', 'Grp']), 'LSST': photo_system('LSST', ['u', 'g'])}}))


def particles(**extra):
    """Three particles that the stars are the children of."""
    return {'age': np.array([9.5, 10.0, 8.0]),
            'feh': np.array([-0.5, -1.0, 0.0]),
            'massinit': np.array([1000.0, 2000.0, 1500.0]),
            'pos3': np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]),
            'vel3': np.array([[0.5, 1.5, 2.5], [3.5, 4.5, 5.5], [6.5, 7.5, 8.5]]),
            **extra}


def kernels():
    return np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])


def prepared_survey(tmp_path, monkeypatch, n_gens=(0,), photo_sys='padova/GAIA__DR2', input_particles=None, **input_kwargs):
    """A survey that is prepared: it has its parameters and its Output, that has nothing to read yet."""
    use_photometry(monkeypatch)
    patch_nbody1(monkeypatch, tmp_path)
    (tmp_path / 'in').mkdir()
    (tmp_path / 'out').mkdir()
    inp = Input(particles() if input_particles is None else input_particles, kernels(),
                **{'input_dir': tmp_path / 'in', 'name': 'sim', 'ngb': 64, **input_kwargs})
    inp.input_sorter = None
    survey = Survey(inp, photo_sys=photo_sys, verbose=False)
    survey._prepare_survey_parameters_and_output(DEFAULT_CMD, n_gens, output_dir=tmp_path / 'out', fsample=1)
    return survey


def log_text(index=0, parameters=(('outputFile', 'survey'), ('fSample', '1'), ('seed', '17052'))):
    """A log as Output reads them: a header, the parameter file between its tags, a footer. The first parameter comes without
    the comment sign that the others have, which is how the parser takes the first line that follows the opening tag."""
    body = "\n".join(f"{key}\t{value}" if i == 0 else f"# {key}\t{value}" for i, (key, value) in enumerate(parameters))
    return (f"# Galaxia (a log made for the tests) number {index}\n# <parameterfile>\n{body}\n# </parameterfile>\n"
            f"# End of the log number {index}\n")


def star_columns(keys, n, partitions=None, parents=None, first=0, **overrides):
    """The columns of n stars for the keys that Output exports: the numbers are simple functions of the index of the star,
    which starts at `first`."""
    i = np.arange(n) + first
    x, y, z = 1.0 + i, 0.5 * i, 0.25 * i
    values = {'px': x, 'py': y, 'pz': z, 'vx': 10.0 * i, 'vy': -5.0 * i, 'vz': 2.0 * i,
              'rad': np.sqrt(x**2 + y**2 + z**2), 'dmod': 10.0 + 0.1 * i, 'mtip': 1.5 + 0.01 * i, 'mact': 0.8 + 0.05 * i,
              'minit': 1.0 + i, 'age': 9.0 + 0.1 * i, 'grav': 4.0 + 0.05 * i, 'feh': -0.5 + 0.01 * i, 'teff': 3.7 + 0.01 * i,
              'lum': 0.1 * i, 'parentid': i % 3 if parents is None else np.asarray(parents),
              'partitionid': 0 * i if partitions is None else np.asarray(partitions), 'partid': 0 * i, 'satid': 0 * i}
    columns = {}
    for k, key in enumerate(keys):
        data = overrides[key] if key in overrides else values.get(key, 15.0 + 0.1 * i + k)       # a magnitude, or another optional key
        columns[key] = np.asarray(data, dtype='int64' if key in INTEGER_KEYS else 'float32')
    return columns


def write_ebf(path, columns, log=None):
    """An EBF file with a dataset for each column and a log."""
    ebf.initialize(str(path))
    for key, array in columns.items():
        ebf.write(str(path), f"/{key}", array, 'a')
    ebf.write(str(path), '/log', log_text() if log is None else log, 'a')
    return path


def write_backend_output(survey, n_stars=5, partitions=None, parents=None, **overrides):
    """
    The files that the backend makes for the survey, one for each generation: {name of the output}.{generation}.ebf, with the same
    number of stars each. The index of the stars goes on from a file to the next, so that no two stars have the same position
    (the partitions and the parents are those of the arguments, in each file). Returns the paths and the columns that were
    written, in the order of the generations.
    """
    output, written = survey.output, []
    for index, generation in enumerate(survey.n_gens):
        path = output.output_dir / f"{output.output_name}.{generation}.ebf"
        columns = star_columns(output.export_keys, n_stars, partitions=partitions, parents=parents, first=n_stars * index, **overrides)
        written.append((write_ebf(path, columns, log=log_text(index)), columns))
    return written


if __name__ == '__main__':
    pass
