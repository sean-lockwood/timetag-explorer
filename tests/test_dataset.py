"""TagDataset on small synthetic exposures: counts, rates, exposure, regions and the processing steps."""
import support  # noqa: F401
import numpy as np
from support import OFF, dataset

from timetag_explorer.constants import EDGE_TRIM_PIXELS, SPATIAL_BIN
from timetag_explorer.frames import Corrections
from timetag_explorer.regions import Box

# 1000 events, 10 per second for 100 s, all in one pixel
n = 1000
t = np.linspace(0, 99.9, n)
ds = dataset(t, np.full(n, 100), np.full(n, 100))
assert ds.exposure() == 100.0 and abs(ds.tmax - 99.9) < 1e-3
img = ds.image(None, 'image', OFF)
assert abs(img.sum() - 10.0) < 1e-9, img.sum()                 # counts per second
img10 = ds.image((0.0, 10.0), 'image', OFF)
assert abs(img10.sum() - 10.0) < 1e-9                          # a time window gives the rate in that window
edges, rate = ds.light_curve(10.0, None, 'image', OFF)
assert np.allclose(rate, 10.0) and len(rate) == 10

# a gap in the good-time intervals: its bins have no exposure (NaN) and the others are still right
gap = dataset(t[(t < 40) | (t >= 60)], np.full(800, 100), np.full(800, 100), gti=((0, 40), (60, 100)))
assert gap.exposure() == 80.0 and abs(gap.exposure(30, 70) - 20.0) < 1e-9
_, r = gap.light_curve(10.0, None, 'image', OFF)
assert np.isnan(r[4]) and np.isnan(r[5]) and np.allclose(r[[0, 1, 2, 3, 6, 7, 8, 9]], 10.0)

# a region: counts inside it only
x = np.where(np.arange(n) % 2 == 0, 100, 300)      # every other event in each of two places
two = dataset(t, x, np.full(n, 100))
inside = two.light_curve(10.0, Box(90, 110, 90, 110), 'image', OFF)[1]
assert np.allclose(inside, 5.0)
assert abs(two.image(None, 'image', OFF).sum() - 10.0) < 1e-9
mask = two.region_mask(Box(90, 110, 90, 110), 'image', OFF)
assert mask.sum() == 500 and mask.dtype == bool

# flat-field weights multiply the counts
weighted = dataset(t, np.full(n, 100), np.full(n, 100), weight=np.full(n, 2.0))
assert abs(weighted.image(None, 'image', Corrections(randcorr=False, geocorr=False, flatcorr=True, timecorr=False)).sum() - 20.0) < 1e-6
assert abs(weighted.image(None, 'image', OFF).sum() - 10.0) < 1e-9     # switched off: unweighted

# the vignetted border is discarded only with the exposure-time correction
edge = dataset(t, np.full(n, EDGE_TRIM_PIXELS - 20), np.full(n, 100))
assert abs(edge.image(None, 'image', OFF).sum() - 10.0) < 1e-9
assert edge.image(None, 'image', Corrections(randcorr=False, geocorr=False, flatcorr=False, timecorr=True)).sum() == 0.0

# the random sub-pixel offset moves events by less than a pixel, keeps the total, and repeats exactly
rand = Corrections(randcorr=True, geocorr=False, flatcorr=False, timecorr=False)
a = dataset(t, np.full(n, 100), np.full(n, 100)).image(None, 'image', rand)
b = dataset(t, np.full(n, 100), np.full(n, 100)).image(None, 'image', rand)
assert np.array_equal(a, b) and abs(a.sum() - 10.0) < 1e-9 and a.shape[0] * SPATIAL_BIN <= 2048

# a window outside the data
assert ds.image((200.0, 300.0), 'image', OFF).sum() == 0.0
print('dataset ok')
