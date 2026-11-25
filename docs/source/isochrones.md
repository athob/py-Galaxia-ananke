# The isochrone library

The photometric systems of `py-Galaxia-ananke` are made of isochrone tables, which its backend, [`Galaxia-ananke`](https://github.com/athob/Galaxia-ananke), carries with it. The directory `GalaxiaData/Isochrones/padova` of the backend has one folder for each system, and the name of that folder is the name of the system: `padova/GAIA__DR2`, for instance. A folder has one file for each of the 34 metallicities, named `output_<Z>.dat` with the mass fraction `Z` to 6 decimals, and a descriptor, `IsoFileDescriptor.txt`. The first line of the descriptor gives the name of the system, the number of columns of the files, the number of columns before the magnitudes, and the number and the names of the magnitudes. Its comment lines say the magnitude system of the magnitudes (Vega, AB or ST) and, for some of the systems, the settings of the tables that the headers of the files do not give. The systems that you add yourself are in the category `py_custom` (see `Photometry.add_isochrone`).

## The systems

There are 16 systems, with 544 files, which take 7.83 GB in all. The CMD column is the version of the [CMD service](http://stev.oapd.inaf.it/cgi-bin/cmd) of the Padova group that made the files.

| System | CMD | Generated | Magnitudes | Size (MB) |
| --- | --- | --- | --- | --- |
| `APO__SDSS` | 2.1 | 11 Aug 2008 | 5 | 118 |
| `Generic__Johnson_UBVRIJHK` | 2.1 | 16 Jun 2008 | 8 | 137 |
| `DCMC` | 2.2 | 29 Jul 2010 | 4 | 113 |
| `Generic__Stroemgren` | 2.2 | 24 May 2010 | 7 | 131 |
| `GAIA__0` | 3.0 | 11 Jul 2017 | 3 | 140 |
| `CTIO__DECam` | 3.6 | 28 Feb 2022 | 6 | 604 |
| `LSST_DP0` | 3.6 | 29 May 2022 | 6 | 604 |
| `Euclid` | 3.7 | 30 Jun 2023 | 6 | 605 |
| `JWST__MIRI` | 3.7 | 1 Jul 2023 | 9 | 651 |
| `JWST__NIRCam` | 3.7 | 1 Jul 2023 | 22 | 839 |
| `2MASS+Spitzer+WISE` | 3.9 | 16 Nov 2026 | 14 | 725 |
| `GAIA__DR2` | 3.9 | 15 Nov 2026 | 3 | 562 |
| `HST__WFC3` | 3.9 | 16 Nov 2026 | 17 | 761 |
| `LSST` | 3.9 | 14 Nov 2026 | 6 | 604 |
| `Pan_STARRS1` | 3.9 | 14 Nov 2026 | 6 | 605 |
| `Roman` | 3.9 | 14 Nov 2026 | 8 | 635 |

## Three generations

* `APO__SDSS`, `Generic__Johnson_UBVRIJHK`, `DCMC` and `Generic__Stroemgren` (CMD 2.1 and 2.2) are the models of Marigo et al. (2008) of the original Galaxia. Their four sets of files have 140 rows (0.08% of them), in 6 files each, with a current mass far below the initial one: it is in the tracks, and the same in the four.
* `GAIA__0` (CMD 3.0, July 2017) is made with PARSEC v1.2S and COLIBRI PR16, and its header does not state the bolometric corrections. Its Gaia passbands are those of the first data release, which the page of the CMD service describes as curves from before the launch. It is kept as it is, as it has no known defect.
* The eleven other systems (CMD 3.6, 3.7 and 3.9) are made with PARSEC v1.2S and the COLIBRI sets S_37 (for 0.008 <= Z <= 0.02), S_35 (for 0.0005 <= Z <= 0.006) and PR16 (for Z <= 0.0002 and Z >= 0.03), and with the bolometric corrections of YBC. The ranges that the form of the CMD service quotes for the three sets leave four metallicities of the grid (Z = 0.0004, 0.007, 0.024 and 0.028) in a gap that its page does not explain: the library has a file for each of them.

## How a system is made

`CTIO__DECam`, `LSST_DP0`, `Euclid`, `JWST__NIRCam` and `JWST__MIRI` were made before, with what became the recipe below. `Pan_STARRS1`, `Roman`, `LSST`, `GAIA__DR2`, `HST__WFC3` and `2MASS+Spitzer+WISE` were made again in November 2026 with it, and only the photometric system changes from one to the other:

| Setting | Value |
| --- | --- |
| Tracks | PARSEC v1.2S, with COLIBRI S_37, S_35 and PR16 as above |
| Thermal pulse cycles | resolution `ninTPC` of 10, the default of the form |
| Long-period variables | the periods and the transition of Trabucchi et al. (2021) |
| Mass loss on the red giant branch | Reimers, with eta = 0.2 |
| Bolometric corrections | YBC with the new Vega (the library of Chen et al. 2019, with the Vega SED of Bohlin et al. 2020) |
| Circumstellar dust | M stars 60% silicate and 40% AlOx; C stars 85% AMC and 15% SiC (Groenewegen 2006) |
| Interstellar extinction | A_V = 0: the tables are intrinsic |
| Initial mass function | Kroupa (2001, 2002) and Kroupa et al. (2013), canonical, corrected for unresolved binaries |
| Ages | 177 for each metallicity, from log age 6.60000 to 10.12005 in steps of 0.02 |
| Metallicities | the 34 values of the mass fraction Z of the grid, from 0.0001 to 0.03 (not of [M/H]) |
| Files | tables that are not compressed, one for each metallicity |

The metallicity `[M/H]` of the form is `log(Z/Zsun)`, with Zsun = 0.019 for the models of Marigo et al. (2008) and 0.0152 for PARSEC. The backend does the same: it takes 0.019 for the four systems of CMD 2.1 and 2.2, and 0.0152 for all the others and for the isochrones that are made in Python.

## What changed in the systems that were made again

* `Pan_STARRS1` had 20 files, from CMD 3.0, and no descriptor: it has 34 now.
* `Roman` has the filters of March 2024: 8 bands, without the prism and the grism that the filters of 2021 had.
* `LSST` has the filters of Release 1.9 (September 2023), and not those of March 2012. The difference, 2012 minus 1.9 for the same stars, is +0.568 mag in `u`, +0.036 in `i`, +0.023 in `g`, +0.018 in `y`, +0.012 in `r` and none in `z`.
* `GAIA__DR2` came from CMD 3.0 (April 2018), whose files had a placeholder in place of the current mass of the stars of 8 solar masses or more (52% of their rows): their gravity, computed from it, was wrong by 1.2 to 1.6 dex. That is gone. The magnitudes of most stars, which are dwarfs, changed with the bolometric corrections: for the stars below 0.2 solar masses (39% of the stars), the old file minus the new one is -0.23 mag in `G`, -0.38 in `G_BP` and -0.18 in `G_RP`.
* `HST__WFC3` has the wide filters of the 2021 update of the calibration of the telescope. For the stars below 0.7 solar masses (86% of the stars), the old file minus the new one goes from -0.36 mag in `F218W` and -0.27 in `F555W` to +0.23 in `F336W` and +0.11 in `F160W`.
* `2MASS+Spitzer+WISE` replaces `GAIA__0+TYCHO+2MASS`, and has 14 bands. The CMD service gives references for the passbands of 2MASS and of Spitzer, and none for those of WISE: where the passbands and the zero points of WISE come from is not documented there.
* `WFIRST`, `WFIRST+HST__WFC3` and `GAIA__0+TYCHO+2MASS` are gone.

## The end of an isochrone

In the systems of CMD 3.0 and later, the isochrone of most ages ends with a row that is not a star (a luminosity of -9.999 and magnitudes of 30 or more), above the last star. A star whose mass is between the two is given the values of that row, a remnant that no magnitude limit lets through. Surveys made before the backend did so contain stars that were neither alive nor dead: the changelog says which.

## Size, and how the data reaches an installation

The data is part of the backend repository, which takes 7.4 GB (as `du -h` prints it), most of it the isochrones. Installing the package fetches the repository, as a shallow clone of the commit that the version pins (about 1 GB of git data next to the files), and copies the data into a cache directory under `sys.prefix`: the README says what to expect in disk space and in network access.
