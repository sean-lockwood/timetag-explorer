# Changelog

All notable changes to the science code. The version in `timetag_explorer/__init__.py` is the version of the
package; `SELECTION_SCRIPT_VERSION` (written into generated scripts and into the `SELVER` FITS keyword) is raised
separately, whenever a generated script would compute or write something different.

How the numbers move:

* **patch** (1.0.x): no result changes (fixes to wording, tests, documentation, performance).
* **minor** (1.x.0): results can change (a corrected equation, a new step, new reference data), or features are added.
* **major** (x.0.0): the Python interface changes incompatibly.

Every change that can alter a computed number says so here, in plain words, so that a published result can be traced
to the version that made it.

## [1.2.0] - Unreleased

**Results change: the detector-edge trim is no longer the same on all four sides** (`SELECTION_SCRIPT_VERSION` 1.2).

* The border discarded with the exposure-time correction was 30 hi-res pixels on every side. It is now, in hi-res pixels
  (two per low-res pixel): **left 30, right 6, bottom 28, top 6** (15, 3, 14 and 3 low-res pixels), as the image frame
  shows the detector (left and right are AXIS1, bottom and top AXIS2). The values come from exposure ODDS22A8Q, which
  fills the field of view with external light; they are `EDGE_TRIM` in `constants.py` and may change again if vignetting
  is still visible.
* Events within the extra 24 pixels on the right and on the top and the extra 2 on the bottom (compared with 30 on
  every side) are now kept, and the exposure map counts those pixels as on the detector, so image and light-curve
  numbers change near those edges. The left edge and the middle of the detector are unchanged.
* Selection scripts: `EDGE_OK` follows the new trim (constants `EDGE_X_LOW`, `EDGE_X_HIGH`, `EDGE_Y_LOW`,
  `EDGE_Y_HIGH`), and the header keyword `SELEDGE` is replaced by `SELEDGEL`, `SELEDGER`, `SELEDGEB`, `SELEDGET`.
* `constants.EDGE_TRIM_PIXELS` and `DETECTOR_LIMITS` are replaced by `EDGE_TRIM` (a named tuple) and
  `DETECTOR_LIMITS_X` / `DETECTOR_LIMITS_Y`.

## [1.1.0] - 2026-10-29

**Selection scripts (`SELECTION_SCRIPT_VERSION` 1.1) write a different file: no computed number changes.**

* The event table's per-column coordinate keywords (`TCTYPn`, `TCUNIn`, `TCRPXn`, `TCRVLn`, `TCDLTn`, `TRPOSn`: the
  sky coordinates of the AXIS1 / AXIS2 columns) were dropped from the `_tagsel.fits` file that a script writes, because
  astropy writes them from attributes of the column and the script rebuilt the columns without them. They are now
  kept. Flags, columns and values are unchanged; files written by a 1.0 script lack these keywords.
* `validation/check_selection_scripts.py` now compares every keyword of the event table, not three, so a keyword lost
  in this way cannot pass unnoticed. (The weekly check found this.)

## [1.0.1] - 2026-10-29

Packaging and documentation only; **no computed number changes** (`SELECTION_SCRIPT_VERSION` stays 1.0).

* First release on PyPI (`pip install timetag-explorer`), published from GitHub Actions.
* Sphinx documentation on Read the Docs: <https://timetag-explorer.readthedocs.io>.
* Project links, classifiers and absolute README links, so the PyPI page works; repository URL in `CITATION.cff`.

## [1.0.0] - 2026-10-28

First public release: the science code that the TIME-TAG Explorer web app uses, split out of the app.

* Event lists (`TagDataset`) read from TIME-TAG files; light curves, count-rate images, exposure maps.
* Processing steps: flat field, random sub-pixel offset, exposure-time correction with detector-edge trim.
* Reference frames: image, target rest, planet rest (positions from JPL Horizons), north-up variants.
* Planet maps (Jupiter, Saturn): System III longitude / latitude and polar projections, with pixel-second exposure.
* Selection scripts (`build_selection_script`, script version 1.0).
* **Geometric distortion correction is disabled** in this release: `gdcorr.py` holds placeholder (zero) coefficients
  and the switch stays off. Real coefficients, from a public reference file, will come in a later release.
