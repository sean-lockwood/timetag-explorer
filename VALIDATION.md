# Validation

Every number the web app shows about the data comes from this package. This page lists each step, what it is meant to
agree with, how it is tested today, and whether a person has reviewed it. **No step has been independently reviewed
yet**: the "Review" column says so until someone has, and who, when and for which version is then recorded here.

How to read the columns:

* **Intended reference**: the document, equation or file the step is meant to follow. These are the author's intentions,
  to be confirmed or corrected by the reviewer.
* **Automated tests**: what runs on every pull request (`tests/`) or every week (`validation/`). Tests on synthetic data
  show that the code does what it says; they do not show that what it says is right.

| # | Step | Where | Intended reference | Automated tests | Review |
|---|---|---|---|---|---|
| 1 | Event positions and binning: 1-indexed hi-res pixels, 4x4 bins, count rate = counts / good exposure | `dataset.TagDataset.image`, `transforms.make_frame_grid` | STIS Instrument Handbook (TIME-TAG mode); the FITS event table | `test_dataset` | not reviewed |
| 2 | Good exposure from the GTI table; light-curve bins with gaps | `dataset.TagDataset.exposure`, `.light_curve` | the `GTI` extension of the TIME-TAG file | `test_dataset` | not reviewed |
| 3 | Flat-field weighting (custom grid-wire flat, sampled on the detector grid) | `dataset.load_flat_weights` | [How the flat was made](https://www.planetarylightshow.com/aurora_processing.html#Flat-Fielding) | `test_dataset` (weights multiply counts) | not reviewed |
| 4 | Random sub-pixel offset (removes the integer-position grid pattern) | `transforms.frame_coordinates` | the usual "randomise within the pixel" practice | `test_dataset` (total kept, repeatable) | not reviewed |
| 5 | **Geometric distortion** | `gdcorr`, `transforms.load_geo_lut` | a public STIS FUV reference file (to come) | `test_gdcorr` | **disabled in 1.0.0** (placeholder coefficients) |
| 6 | Exposure-time correction and the detector-edge trim (left 30, right 6, bottom 28, top 6 pixels) | `dataset.TagDataset.exposure_map`, `transforms.compute_exposure_map` | an empirical choice for the vignetted border (from exposure ODDS22A8Q, which fills the field of view with external light) | `test_dataset` (trim on and off, each side, the exposure map), `test_selection_script` (EDGE_OK equals the Explorer's mask) | not reviewed |
| 7 | Target rest frames: removal of the commanded slew, rotation by ORIENTAT for north-up | `transforms.frame_shift`, `frames.load_frame_info` | the SPT file's motion-tracking (`MT_LV`) levels; HST pointing keywords | `test_frames` (header parsing) | not reviewed |
| 8 | Planet rest frame: planet (and moon) positions from JPL Horizons, converted to detector offsets | `horizons`, `transforms` | JPL Horizons ephemerides via `astroquery`; the header WCS | `check_selection_scripts` (frames only) | not reviewed |
| 9 | Planet orientation: IAU pole and rotation elements, light-time correction, TDB, System III | `planet.planet_orientation` | IAU WGCCRE report (Archinal et al. 2018) for Jupiter and Saturn | `test_planet` (rotation rate, apparent size, leap seconds) | not reviewed |
| 10 | Sky to planet-surface projection on an oblate spheroid; limb outline | `planet.screen_to_surface`, `.surface_to_screen`, `.limb_outline` | orthographic projection of a spheroid; planetographic latitude, west longitude | `test_planet` (round trips, sub-observer point) | not reviewed |
| 11 | Map frames (longitude/latitude, polar) and surface area of a map bin | `planet.map_plane_coordinates`, `.map_bin_area` | spherical and spheroidal geometry | `test_planet` (whole-sphere and polar-cap areas) | not reviewed |
| 12 | Map brightness: counts over the pixel-seconds that saw each bin | `transforms.map_exposure_from_knots`, `dataset.TagDataset.map_exposure` | derived in the code comments | none yet | not reviewed |
| 13 | Region membership (box, circle, polygon) | `regions` | half-open boxes; closed circles; even-odd polygon test | `test_regions` | not reviewed |
| 14 | Selection scripts: the flags they write equal the app's own masks, event by event | `selection_script` | the app's masks (steps 1, 6, 7, 13) | `test_selection_script` (script is valid, versioned); `validation/check_selection_scripts.py` (weekly, real data) | not reviewed |

## Known gaps

* Step 5 is disabled, so nothing about the distortion correction is validated or used.
* Steps 7, 8 and 12 have no known-answer test yet; they are only exercised end to end.
* The tests run on synthetic data or on one real exposure; no comparison with an independent code (calstis, for
  example) has been made.

## Recording a review

A reviewer adds a row to the table below and, if they found problems, opens issues. A review applies to the version
written beside it, and a later change to that step makes it stale until it is repeated.

| Step(s) | Reviewer | Date | Version reviewed | Outcome |
|---|---|---|---|---|
| (none yet) | | | | |
