"""Numbers that the science code shares: bin sizes, limits and version numbers.

Nothing here depends on a web server or on where files are kept (see `config` for that)."""
from __future__ import annotations

from typing import NamedTuple

# Event-table columns (raw _tag files). For _corrtag files you might use
# 'XCORR' / 'YCORR' or 'XFULL' / 'YFULL' instead.
TIME_COL, X_COL, Y_COL = 'TIME', 'AXIS1', 'AXIS2'

# Flat field: a 2048 x 2048 image in ext 0, indexed [AXIS2, AXIS1] in hi-res
# pixels. The file already holds the multiplicative weight, so each event gets
# weight = file value at its pixel (no inversion). TAG AXIS1/AXIS2 are 1-indexed.
# The flat is optional: without the file, counts are unweighted and its switch is disabled.
FLAT_URL = 'https://www.planetarylightshow.com/static/data/combined_flat.fits'   # where the flat field can be downloaded
FLAT_SHAPE = (2048, 2048)

SPATIAL_BIN = 4          # detector pixels per image bin
DEFAULT_TIME_BIN = 3     # seconds per light-curve bin
MAX_LIGHT_CURVE_BINS = 100_000   # a smaller bin would only exhaust memory
MIN_TIME_BIN = 0.001             # seconds

TIME_DECIMALS = 1        # time selections are rounded to 0.1 s

RANDOM_SEED = 0          # fixed, so repeated renders of one dataset agree

# Version of the script that the "Selection code" button writes. People keep these scripts and the FITS
# files they produce, and the number is written into both (docstring, SCRIPT_VERSION, the SELVER header
# keyword, HISTORY), so raise it whenever build_selection_script changes what the script computes or
# writes (not for wording changes in the dialog).
SELECTION_SCRIPT_VERSION = '2.2'   # 2.2: own header keywords (TTERAND ... TTEEDGE: DONE/OMIT, TTEFFILE, TTETSEL, TTERSEL) in one labelled block; 2.1: north-up frame columns only with INCLUDE_NORTH_UP, no WEIGHT column, light-curve and count-rate examples; 2.0: the script calls the package (tagsel); exact region test; XCORR/YCORR, frame and weight columns; 1.2: EDGE_OK uses a different trim on each side; 1.1: column coordinate keywords kept
DEFAULT_TOOL_URL = 'https://explorer.planetarylightshow.com'   # TAG_EXPLORER_URL, when set, takes its place
MAX_CACHED_FRAMES = 2    # per dataset; each cached frame costs 4 bytes per event
EVENT_CHUNK = 1 << 20    # events handled at a time, so temporary arrays stay small whatever the file size


def _chunks(n: int, size: int | None = None):
    """(lo, hi) slices of range(n), EVENT_CHUNK events each."""
    size = size or EVENT_CHUNK
    for lo in range(0, n, size):
        yield lo, min(lo + size, n)

# Exposure-time correction: a frame pixel seen for less than this fraction of the
# good exposure (the edges of a drifting field of view) is left blank, not amplified.
TIME_CORR_MIN_FRACTION = 0.1
MAP_BIN_DEGREES = 0.25   # bin size of the planet maps (degrees of longitude / latitude / colatitude)
MAP_SEGMENTS = 8         # the planet's rotation is treated as this many straight segments when mapping exposure
MAP_SLICE_SECONDS = 60.0 # events are deprojected in time slices this long
EXPOSURE_SUBSAMPLE = 2   # sub-samples per bin and axis for a moving footprint (static ones use twice as many)
PLANET_SEGMENTS = 8      # the planet's path is treated as this many straight segments


class EdgeTrim(NamedTuple):
    """Detector border discarded, in hi-res pixels (two per low-res pixel), as the image frame shows the detector:
    `left` and `right` are the low and high AXIS1 (x) edges, `bottom` and `top` the low and high AXIS2 (y) edges."""
    left: int
    right: int
    bottom: int
    top: int


# The detector's edges are vignetted, by different amounts on each side. With the exposure-time correction on,
# events within these distances of an edge are discarded, and the exposure map treats the detector as the rectangle
# with that border removed. Set a side to 0 to keep it, all four to 0 to keep the whole detector.
# (From exposure ODDS22A8Q, which fills the field of view with external light: 15, 3, 14 and 3 low-res pixels.)
EDGE_TRIM = EdgeTrim(left=30, right=6, bottom=28, top=6)

# The vignetting model (vignette.py) was measured through this aperture only, so it is applied to no other.
VIGNETTE_APERTURE = 'F25SRF2'
# The rectangle that is kept, in 0-indexed hi-res pixel coordinates (a pixel p covers [p - 0.5, p + 0.5]).
DETECTOR_LIMITS_X = (EDGE_TRIM.left - 0.5, 2047.5 - EDGE_TRIM.right)
DETECTOR_LIMITS_Y = (EDGE_TRIM.bottom - 0.5, 2047.5 - EDGE_TRIM.top)

TimeSel = tuple[float, float]
