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

## [1.4.1] - Unreleased

**Fix: no result changes** (`SELECTION_SCRIPT_VERSION` stays 2.0).

* A selection script for a region drawn in the Explorer failed at its `REGION = ...` line unless numpy had been imported
  first, because the region held numpy floats, which NumPy 2 writes as `np.float64(...)`. Regions are now written with
  plain floats (`regions.region_source`), in the script and in the file's HISTORY. 1.4.0 scripts that use a box, circle
  or free-form region from the app need this version.

## [1.4.0] - Unreleased

**Selection scripts are rebuilt on the package** (`SELECTION_SCRIPT_VERSION` 2.0). No number the Explorer shows changes;
the script and the file it writes do.

* The generated script is now short: the user's choices (time window, region, frame, corrections) and a call to the new
  `timetag_explorer.tagsel.make_tagsel`, which does the work with the package's own code instead of a re-implementation
  of the arithmetic. It names the exact `timetag-explorer` version the Explorer ran, gives the `pip install` line for
  it (a PyPI release, or the git commit or archive address of a branch the Explorer is testing; read from pip's
  `direct_url.json`), and stops if another version is installed unless `ALLOW_OTHER_VERSION = True`.
* **The region is tested at each event's own position**, not on the Explorer's 4-pixel bins. Events within about two
  pixels of a region's edge can differ from the highlighted bins. (Regions drawn on planet maps, which are experimental,
  are not carried over; the script says so.)
* The output (`<rootname>_tagsel.fits`) keeps every original column, `AXIS1` and `AXIS2` included, unchanged, and adds:
  `XCORR`, `YCORR` (the image-frame position with the random sub-pixel offset and/or the geometric distortion
  correction applied, when either is on); a pair per reference frame from those positions (`XREST`/`YREST`,
  `XRESTN`/`YRESTN`, and for moon observations `XPLAN`/`YPLAN`, `XPLANN`/`YPLANN`), counted from 1 like AXIS1/AXIS2;
  the weights `FLATW`, `VIGW` (vignetting model) and `WEIGHT` (their product); and the flags `SEL_TIME`, `SEL_REGION`,
  `EDGE_OK` as before.
* **Primary-header keywords** (ext 0) record what was done: `RANDCORR`, `GEOCORR`, `FLATCORR`, `VIGNCORR`, `EXPTCORR`
  (`PERFORM` or `OMIT`; `PERFORM` only when the step took effect: the flat field was found, the vignetting model applies to
  the aperture), `TIMESEL` and `SPCSEL` (True when a time range / a region was selected), `SELTMIN`/`SELTMAX`,
  `SELREG`/`SELFRAME`, `RANDSEED`, the vignetting ramps (`VIGBSTRT` ... `VIGTLOSS`), the edge trim (`SELEDGEL` ...),
  and `SELVER` and the new `TTEVER` (package version) and `SELURL`. These replace the extension-1 keywords of the 1.x
  scripts (`SELRAND`, `SELGEO` and the others). Where the calibration pipeline had already put one of the five step
  keywords in the file, its original value is kept in the HISTORY.
* **A selection is optional**: with no time range and no region the script still applies the corrections and writes the
  positions and weights (`TIMESEL` and `SPCSEL` are False and there are no `SEL_TIME` / `SEL_REGION` columns).
* The script downloads the SPT file (for the rest frames) and the flat field as well as the TAG file.
* New `TagDataset.vignette_gain()`. The self-contained generator, with its own copy of the geometric distortion
  polynomial, is gone; scripts written by version 1.x keep working as they are.
* `validation/check_selection_scripts.py` now checks the positions and weights as well as the flags.

## [1.3.0] - Unreleased

**New, experimental, off by default: a vignetting model for the F25SRF2 aperture.** No result changes unless it is switched on
(`SELECTION_SCRIPT_VERSION` is unchanged: selection scripts select events and apply no weights).

* New `vignette` module and switch `Corrections.vigcorr` (default `False`). With it on, each event is weighted by
  1 / (bottom sensitivity x top sensitivity) at its row, where each sensitivity is a linear ramp: 1 from the inside up
  to the ramp's start row, falling to 1 - loss at its end row, and on at the same slope beyond it (never below 0.2).
  The numbers are `Corrections.vigbot_start / _end / _loss` (default 430, 0, 0.25) and `vigtop_start / _end / _loss`
  (default 1550, 2047, 0.23), rows being 0-indexed hi-res detector pixels as the image frame shows them. The defaults
  are empirical (exposure ODDS22A8Q) and they are still under study.
* The model is applied only to exposures through **F25SRF2** (the primary-header `APERTURE`, new attribute
  `TagDataset.aperture`; `TagDataset.vignette_applies()`). For any other aperture the switch does nothing.
* It is applied **after the geometric distortion step**: with `geocorr` on, the row is the recorded row plus the
  distortion table's row offset at the recorded pixel (the same lookup the frame coordinates use). The weights multiply
  the flat-field weights and apply to images and to light curves; the exposure-time (edge) trim is unchanged.

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
