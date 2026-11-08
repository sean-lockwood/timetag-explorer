"""Several time windows in an image, event masks (a region filter, in any frame), selection masks for rows of a region
table (drawn or made from others), light curves from masks, and the curve normalization the Explorer plots."""
import support  # noqa: F401
import numpy as np
from support import OFF, dataset

from timetag_explorer.analysis import curve_peak, normalized_curve
from timetag_explorer.dataset import time_windows
from timetag_explorer.frames import Corrections
from timetag_explorer.regions import Box
from timetag_explorer.selections import Intervals, RegionSelection, RegionTable, TimeRange, TimeTable


def refused(call, text):
    try:
        call()
    except ValueError as err:
        assert text in str(err), (text, str(err))
    else:
        raise AssertionError(f'expected an error containing {text!r}')


# ---- time windows ---------------------------------------------------------------------------------------------------
inf = np.inf
assert time_windows(None) == [(-inf, inf)]
assert time_windows((10, 20)) == [(10.0, 20.0)]
assert time_windows([(30, 40), (10, 20)]) == [(10.0, 20.0), (30.0, 40.0)]                 # sorted
assert time_windows([(10, 25), (20, 40)]) == [(10.0, 40.0)]                               # overlapping ones merge
assert time_windows([(10, 10), (30, 20)]) == []                                           # empty ones go
assert time_windows(Intervals([[1, 2], [5, 6]])) == [(1.0, 2.0), (5.0, 6.0)]

# ---- images ---------------------------------------------------------------------------------------------------------
n = 400
rng = np.random.default_rng(3)
t = np.sort(rng.uniform(0, 99.9, n))
x = rng.integers(100, 1900, n)
y = rng.integers(100, 1900, n)
ds = dataset(t, x, y)
total = ds.image(None, 'image', OFF).sum()
exposure_all = ds.exposure(0, 100)
assert abs(total * exposure_all - n) < 1e-6                                                # counts / total time, summed

a, b = ds.image((10, 30), 'image', OFF), ds.image((50, 90), 'image', OFF)
both = ds.image([(10, 30), (50, 90)], 'image', OFF)
assert abs(both.sum() * 60.0 - (a.sum() * 20.0 + b.sum() * 40.0)) < 1e-6                   # counts over the total time of both
assert np.allclose(both * 60.0, a * 20.0 + b * 40.0)
assert np.array_equal(ds.image((10, 30), 'image', OFF), ds.image([(10, 30)], 'image', OFF))
assert not ds.image([], 'image', OFF).any()                                               # no window, no events

# an event mask keeps only those events; the rest of the image is zero
grid = ds.grid('image', OFF)
box = Box(500.0, 1000.0, 500.0, 1000.0)                                                     # array coordinates
mask = ds.region_mask(box, 'image', OFF)
assert 0 < mask.sum() < n
filtered = ds.image(None, 'image', OFF, event_mask=mask)
assert abs(filtered.sum() * exposure_all - mask.sum()) < 1e-6
assert not filtered[~ds.region_grid(box, 'image', OFF)].any()
assert (filtered <= ds.image(None, 'image', OFF)).all()
# a mask and several windows together
sel = ds.image([(10, 30), (50, 90)], 'image', OFF, event_mask=mask)
assert abs(sel.sum() * 60.0 - np.sum(mask & ((ds.time >= 10) & (ds.time < 30) | (ds.time >= 50) & (ds.time < 90)))) < 1e-6
# with the exposure-time correction the stub leaves the exposure as it is (see TagDataset.filter_exposure_map)
on = Corrections(randcorr=False, geocorr=False, flatcorr=False, timecorr=True)
emap = ds.exposure_map(None, 'image', on).ravel()
assert ds.filter_exposure_map(emap, mask, None, 'image', on) is emap
assert ds.image(None, 'image', on, event_mask=mask).sum() <= ds.image(None, 'image', on).sum() + 1e-12

# ---- selection masks ------------------------------------------------------------------------------------------------
rows = RegionTable([
    RegionSelection.from_array(Box(500.0, 1000.0, 500.0, 1000.0), 'image', name='big'),
    RegionSelection.from_array(Box(600.0, 900.0, 600.0, 900.0), 'image', name='small'),
    RegionSelection.from_array(Box(100.0, 1500.0, 100.0, 1500.0), 'rest', name='wide'),     # another frame (no slew here: the same pixels)
    RegionSelection(name='ring', expr='big - small'),
    RegionSelection(name='both', expr='big & wide'),
])
big, small = ds.selection_mask(rows, 'big', OFF), ds.selection_mask(rows, 'small', OFF)
assert np.array_equal(big, mask) and ds.selection_mask(rows, 0, OFF).sum() == big.sum()   # by name or by bit
assert np.array_equal(ds.selection_mask(rows, 'ring', OFF), big & ~small)
wide = ds.selection_mask(rows, 'wide', OFF)
assert np.array_equal(ds.selection_mask(rows, 'both', OFF), big & wide)                   # rows in two frames combine per event
assert big.sum() > small.sum() > 0 and wide.sum() > big.sum()
# an image in the 'image' frame filtered by a region drawn in the 'rest' frame
cross = ds.image(None, 'image', OFF, event_mask=wide)
assert abs(cross.sum() * exposure_all - wide.sum()) < 1e-6
rest_image = ds.image(None, 'rest', OFF, event_mask=big)
assert abs(rest_image.sum() * exposure_all - big.sum()) < 1e-6                              # the same events, shown in the other frame
refused(lambda: ds.selection_mask(RegionTable([RegionSelection(name='m', region=Box(0, 1, 0, 1), frame='map_rect')]), 'm', OFF),
        'unavailable')                                                                      # a frame that cannot be built

# ---- light curves from masks ----------------------------------------------------------------------------------------
edges, whole = ds.light_curve(5.0, None, 'image', OFF)
edges, by_region = ds.light_curve(5.0, box, 'image', OFF)
edges, by_mask = ds.light_curve(5.0, None, 'image', OFF, event_mask=mask)
assert np.allclose(by_region, by_mask, equal_nan=True)
both_masks = ds.light_curve(5.0, box, 'image', OFF, event_mask=ds.selection_mask(rows, 'small', OFF))[1]
assert np.nansum(both_masks) <= np.nansum(by_mask) + 1e-12
assert np.nansum(by_mask) <= np.nansum(whole) + 1e-12

# ---- the table: several times, display state --------------------------------------------------------------------------
times = TimeTable([TimeRange(name='a', start=10.0, stop=30.0), TimeRange(name='b', start=50.0, stop=90.0),
                   TimeRange(name='ab', expr='a | b')])
assert times.union(['a', 'b']).rows.tolist() == [[10.0, 30.0], [50.0, 90.0]] and not len(times.union([]).rows)
assert time_windows(times.union(['ab'])) == [(10.0, 30.0), (50.0, 90.0)]
assert np.allclose(ds.image(times.union(['a', 'b']), 'image', OFF), both)
display = RegionTable([RegionSelection(name='x', region=Box(1, 2, 1, 2), frame='image'),
                       RegionSelection(name='y', region=Box(1, 2, 1, 2), frame='image', curve=False, visible=False)])
cols = display.to_columns()
assert cols['CURVE'].tolist() == [True, False] and cols['VISIBLE'].tolist() == [True, False]
assert [(r.visible, r.curve) for r in RegionTable.from_columns(cols)] == [(True, True), (False, False)]
old = {k: v for k, v in cols.items() if k != 'CURVE'}                                      # a file from before selection script 3.1
assert [r.curve for r in RegionTable.from_columns(old)] == [True, True]

# ---- curve normalization ----------------------------------------------------------------------------------------------
curve = np.array([0.0, 2.0, np.nan, 8.0, 4.0])
scaled, peak = normalized_curve(curve)
assert peak == 8.0 and np.nanmax(scaled) == 1.0 and np.isnan(scaled[2]) and scaled[1] == 0.25
assert curve_peak([0.0, 0.0]) == 1.0 and curve_peak([np.nan]) == 1.0 and curve_peak([]) == 1.0
print('filters ok')

# memory: which frames have a grid built already
fresh = dataset(t, x, y)
assert not fresh.has_grid('image') and not fresh.has_grid('rest')
fresh.grid('image', OFF)
assert fresh.has_grid('image') and not fresh.has_grid('rest')
