"""Shared by the tests: make the package importable from a checkout (no install needed) and build small datasets."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from timetag_explorer.dataset import TagDataset  # noqa: E402
from timetag_explorer.frames import Corrections, FrameInfo  # noqa: E402

OFF = Corrections(randcorr=False, geocorr=False, flatcorr=False, timecorr=False)   # every processing step off


def dataset(time, x, y, gti=((0.0, 100.0),), weight=None, name='odtest01q_tag.fits'):
    """A TagDataset made from arrays (no FITS file): events at times `time` (s) and 1-indexed pixels (x, y)."""
    time = np.asarray(time, dtype=np.float32)
    order = np.argsort(time, kind='stable')
    return TagDataset(
        name=name, time=time[order], x=np.asarray(x, dtype=np.int16)[order], y=np.asarray(y, dtype=np.int16)[order],
        gti=np.asarray(gti, dtype=float), rootname=name.split('_')[0],
        weight=None if weight is None else np.asarray(weight, dtype=np.float32)[order],
        frame_info=FrameInfo(orientat=0.0, slew_angle=0.0, slew_rate=0.0))
