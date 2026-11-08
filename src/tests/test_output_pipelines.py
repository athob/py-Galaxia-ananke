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
Tests of the post-processing pipelines of Output, which are functions of a table of stars: the conversion of Cartesian coordinates
relative to the Sun into Galactic ones, of Galactic into equatorial ones (and back), and the last conversions (parallax, and the
linear temperature and luminosity). They are checked against what physics says (a star at 1 kpc whose tangential velocity is
4.74 km/s moves by 1 milliarcsecond per year), against formulas that do not use astropy, and against the known place of the
center and of the pole of the Galaxy.
"""
import numpy as np
import pandas as pd
import pytest

from ..galaxia_ananke import Output

K = 4.740470463533348             # km/s of tangential velocity for 1 mas/yr at 1 kpc
cartesian_to_galactic = Output._Output__pp_convert_cartesian_to_galactic
galactic_to_icrs = Output._Output__pp_convert_galactic_to_icrs
icrs_to_galactic = Output._Output__pp_convert_icrs_to_galactic
last_conversions = Output._Output__pp_last_conversions


def frame(**columns):
    return pd.DataFrame({key: np.asarray(value, dtype=float) for key, value in columns.items()})


def stars(x, y, z, vx=None, vy=None, vz=None):
    zeros = np.zeros(len(x))
    return frame(px=x, py=y, pz=z, vx=zeros if vx is None else vx, vy=zeros if vy is None else vy, vz=zeros if vz is None else vz)


def expected_motion(df):
    """The proper motions and the radial velocity of stars from their Cartesian positions and velocities, by hand."""
    x, y, z, vx, vy, vz = (df[key].to_numpy() for key in ('px', 'py', 'pz', 'vx', 'vy', 'vz'))
    r, rho = np.sqrt(x**2 + y**2 + z**2), np.sqrt(x**2 + y**2)
    v_l = (-y * vx + x * vy) / rho
    v_b = (-x * z * vx - y * z * vy + rho**2 * vz) / (r * rho)
    return v_l / (K * r), v_b / (K * r), (x * vx + y * vy + z * vz) / r


def random_stars(n=12, seed=11):
    rng = np.random.default_rng(seed)
    positions = rng.uniform(1.0, 6.0, (n, 3)) * rng.choice([-1, 1], (n, 3))
    velocities = rng.normal(0.0, 120.0, (n, 3))
    return stars(*positions.T, *velocities.T)


# -- from Cartesian to Galactic -----------------------------------------------------

def test_stars_on_the_axes_have_the_longitude_and_the_latitude_of_their_axis():
    df = stars([1, 0, -1, 0, 1], [0, 1, 0, -1, 1], [0, 0, 0, 0, 0])
    cartesian_to_galactic(df)
    assert df['glon'].tolist() == pytest.approx([0, 90, 180, -90, 45])
    assert df['glat'].tolist() == pytest.approx([0, 0, 0, 0, 0])
    assert df['rad'].tolist() == pytest.approx([1, 1, 1, 1, 2 ** 0.5])


def test_a_star_above_the_plane_has_a_positive_latitude_and_the_pole_is_at_90():
    df = stars([1, 0], [0, 0], [1, 2])
    cartesian_to_galactic(df)
    assert df['glat'].tolist() == pytest.approx([45, 90]) and df['rad'].tolist() == pytest.approx([2 ** 0.5, 2])


def test_the_longitude_is_between_minus_180_and_180():
    df = stars([-1, -1, 1], [1, -1, -1], [0, 0, 0])
    cartesian_to_galactic(df)
    assert df['glon'].tolist() == pytest.approx([135, -135, -45])


def test_a_star_at_one_kpc_that_moves_across_the_sky_by_4_74_km_per_s_moves_by_1_mas_per_year():
    df = stars([1, 1], [0, 0], [0, 0], vy=[K, 0], vz=[0, K])
    cartesian_to_galactic(df)
    assert df['mul'].tolist() == pytest.approx([1.0, 0.0], abs=1e-6)             # along the longitude
    assert df['mub'].tolist() == pytest.approx([0.0, 1.0], abs=1e-6)             # along the latitude
    assert df['vr'].tolist() == pytest.approx([0.0, 0.0], abs=1e-9)


def test_the_proper_motion_goes_down_with_the_distance():
    df = stars([1, 2, 4], [0, 0, 0], [0, 0, 0], vy=[K, K, K])
    cartesian_to_galactic(df)
    assert df['mul'].tolist() == pytest.approx([1.0, 0.5, 0.25], rel=1e-6)


def test_the_velocity_along_the_line_of_sight_is_the_radial_velocity_and_moves_nothing_across():
    df = stars([1, 0], [0, 3], [0, 0], vx=[5, 0], vy=[0, -7])
    cartesian_to_galactic(df)
    assert df['vr'].tolist() == pytest.approx([5.0, -7.0]) and df['mul'].tolist() == pytest.approx([0.0, 0.0], abs=1e-6)
    assert df['mub'].tolist() == pytest.approx([0.0, 0.0], abs=1e-6)


def test_the_angles_the_distances_and_the_motions_of_any_star_are_those_of_the_formulas():
    df = random_stars()
    expected = expected_motion(df)
    x, y, z = (df[key].to_numpy() for key in ('px', 'py', 'pz'))
    cartesian_to_galactic(df)
    assert df['glon'].to_numpy() == pytest.approx(np.degrees(np.arctan2(y, x)), abs=1e-9)
    assert df['glat'].to_numpy() == pytest.approx(np.degrees(np.arcsin(z / np.sqrt(x**2 + y**2 + z**2))), abs=1e-9)
    assert df['rad'].to_numpy() == pytest.approx(np.sqrt(x**2 + y**2 + z**2))
    assert df['mul'].to_numpy() == pytest.approx(expected[0], rel=1e-6) and df['mub'].to_numpy() == pytest.approx(expected[1], rel=1e-6)
    assert df['vr'].to_numpy() == pytest.approx(expected[2], rel=1e-6)


def test_a_table_without_stars_gets_the_columns_all_the_same():
    df = stars([], [], [])
    cartesian_to_galactic(df)
    assert {'glon', 'glat', 'rad', 'mul', 'mub', 'vr'} <= set(df.columns) and len(df) == 0


# -- from Galactic to equatorial ----------------------------------------------------

def test_the_center_and_the_pole_of_the_galaxy_are_where_they_are_on_the_sky():
    df = frame(glon=[0, 0], glat=[0, 90], rad=[8, 8], mul=[0, 0], mub=[0, 0], vr=[0, 0])
    galactic_to_icrs(df)
    assert df['ra'].tolist() == pytest.approx([266.405, 192.859], abs=2e-3)
    assert df['dec'].tolist() == pytest.approx([-28.936, 27.128], abs=2e-3)
    assert df['mura'].tolist() == pytest.approx([0, 0], abs=1e-9) and df['mudec'].tolist() == pytest.approx([0, 0], abs=1e-9)


def test_the_proper_motion_is_turned_but_not_changed_in_size():
    rng = np.random.default_rng(5)
    df = frame(glon=rng.uniform(-170, 170, 15), glat=rng.uniform(-80, 80, 15), rad=rng.uniform(0.5, 9, 15),
               mul=rng.normal(0, 5, 15), mub=rng.normal(0, 5, 15), vr=rng.normal(0, 90, 15))
    galactic_to_icrs(df)
    assert np.hypot(df['mura'], df['mudec']).to_numpy() == pytest.approx(np.hypot(df['mul'], df['mub']).to_numpy(), rel=1e-6)


def test_the_equatorial_coordinates_of_a_table_without_stars_are_there_all_the_same():
    df = frame(glon=[], glat=[], rad=[], mul=[], mub=[], vr=[])
    galactic_to_icrs(df)
    assert {'ra', 'dec', 'mura', 'mudec'} <= set(df.columns) and len(df) == 0


# -- from equatorial to Galactic ----------------------------------------------------

def test_going_to_equatorial_and_back_gives_the_same_galactic_coordinates_and_motions():
    rng = np.random.default_rng(9)
    df = frame(glon=rng.uniform(-170, 170, 15), glat=rng.uniform(-80, 80, 15), rad=rng.uniform(0.5, 9, 15),
               mul=rng.normal(0, 5, 15), mub=rng.normal(0, 5, 15), vr=rng.normal(0, 90, 15))
    original = df.copy()
    galactic_to_icrs(df)
    for key in ('glon', 'glat', 'mul', 'mub'):
        del df[key]
    icrs_to_galactic(df)
    assert df['glon'].to_numpy() == pytest.approx(original['glon'].to_numpy(), abs=1e-6)
    assert df['glat'].to_numpy() == pytest.approx(original['glat'].to_numpy(), abs=1e-6)
    assert df['mul'].to_numpy() == pytest.approx(original['mul'].to_numpy(), abs=1e-6)
    assert df['mub'].to_numpy() == pytest.approx(original['mub'].to_numpy(), abs=1e-6)


def test_the_galactic_center_in_equatorial_coordinates_is_at_the_origin_of_the_galactic_ones():
    df = frame(ra=[266.405], dec=[-28.936], mura=[0], mudec=[0])
    icrs_to_galactic(df)
    assert df['glon'].tolist() == pytest.approx([0], abs=2e-3) and df['glat'].tolist() == pytest.approx([0], abs=2e-3)


def test_the_galactic_coordinates_of_a_table_without_stars_are_there_all_the_same():
    df = frame(ra=[], dec=[], mura=[], mudec=[])
    icrs_to_galactic(df)
    assert {'glon', 'glat', 'mul', 'mub'} <= set(df.columns) and len(df) == 0


# -- the last conversions -----------------------------------------------------------

def test_the_parallax_is_the_inverse_of_the_distance_and_the_logs_are_made_linear():
    df = frame(rad=[0.5, 2.0, 8.0], teff=[3.7, 4.0, 3.5], lum=[0.0, 1.0, -2.0])
    last_conversions(df)
    assert df['pi'].tolist() == pytest.approx([2.0, 0.5, 0.125])                       # mas, for kpc
    assert df['teff'].tolist() == pytest.approx([10 ** 3.7, 10 ** 4.0, 10 ** 3.5])
    assert df['lum'].tolist() == pytest.approx([1.0, 10.0, 0.01])
    assert df['rad'].tolist() == [0.5, 2.0, 8.0]                                         # what was not asked for is left


if __name__ == '__main__':
    pass
