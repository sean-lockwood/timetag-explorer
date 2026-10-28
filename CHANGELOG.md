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

## [1.0.0] - Unreleased

First public release: the science code that the TIME-TAG Explorer web app uses, split out of the app.

* Event lists (`TagDataset`) read from TIME-TAG files; light curves, count-rate images, exposure maps.
* Processing steps: flat field, random sub-pixel offset, exposure-time correction with detector-edge trim.
* Reference frames: image, target rest, planet rest (positions from JPL Horizons), north-up variants.
* Planet maps (Jupiter, Saturn): System III longitude / latitude and polar projections, with pixel-second exposure.
* Selection scripts (`build_selection_script`, script version 1.0).
* **Geometric distortion correction is disabled** in this release: `gdcorr.py` holds placeholder (zero) coefficients
  and the switch stays off. Real coefficients, from a public reference file, will come in a later release.
