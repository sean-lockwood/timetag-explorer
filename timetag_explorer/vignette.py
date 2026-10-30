"""Vignetting model (experimental): a sensitivity drop towards the bottom and the top of the detector.

The combined flat field leaves the rows near the bottom and top edges under-corrected for the F25SRF2 aperture, a
vignetting component the flat does not model. This module is the simplest model of it: two linear ramps, one per edge.

A ramp is given by three numbers: the row where it starts (the sensitivity is 1 there and on the inside), the row where
the loss is reached, and the loss (a fraction). From the start row towards the end row the sensitivity falls linearly
from 1 to 1 - loss, and keeps falling at that slope beyond the end row (but never below MIN_SENSITIVITY). Rows are
0-indexed hi-res detector pixels, as the image frame's y axis shows them. The bottom ramp has start > end (the loss
grows downwards), the top ramp start < end (it grows upwards).

An event is weighted by 1 / (bottom sensitivity x top sensitivity) at its row. The row is where the event lies after
the geometric distortion correction (when that is on), the step before this one in the processing order. The numbers
are Corrections.vigbot_* and Corrections.vigtop_*, applied when Corrections.vigcorr is on and the exposure was taken
through VIGNETTE_APERTURE.
"""
from __future__ import annotations

import numpy as np

from .constants import VIGNETTE_APERTURE  # noqa: F401  (re-exported: the one aperture the model describes)

MIN_SENSITIVITY = 0.2   # a ramp is continued beyond its end row, but never trusted below this (weights stay finite)


def sensitivity(row, start: float, end: float, loss: float):
    """Relative sensitivity of one ramp at detector row(s) `row`: 1 on the inside of `start`, falling linearly to
    1 - loss at `end`, and on at the same slope (not below MIN_SENSITIVITY)."""
    if start == end:
        raise ValueError(f'a ramp needs a start row different from its end row (both are {start})')
    along = np.maximum((start - np.asarray(row, dtype=np.float64)) / (start - end), 0.0)
    return np.maximum(1.0 - loss * along, MIN_SENSITIVITY)


def event_gain(x, y, bottom: tuple[float, float, float], top: tuple[float, float, float],
               y_offset=None) -> np.ndarray:
    """Weight for events at AXIS1 / AXIS2 values `x`, `y` (1-indexed): 1 / sensitivity at their row, as float32.

    `bottom` and `top` are (start, end, loss) of each ramp. `y_offset` is the geometric-distortion table of row
    offsets (hi-res pixels, indexed [row, column]) or None when that correction is off; the offset at the event's
    recorded pixel is added to its row, as the frame coordinates do.
    """
    row = np.asarray(y, dtype=np.float64) - 1.0
    if y_offset is not None:
        ix = np.clip(np.asarray(x, dtype=np.int32) - 1, 0, 2047)
        iy = np.clip(np.asarray(y, dtype=np.int32) - 1, 0, 2047)
        row = row + y_offset[iy, ix]
    s = sensitivity(row, *bottom) * sensitivity(row, *top)
    return (1.0 / s).astype(np.float32)
