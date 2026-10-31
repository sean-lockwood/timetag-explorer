"""Vignetting model (vignette.py): linear sensitivity ramps at the bottom and the top of the detector (F25SRF2 only)."""
import support  # noqa: F401
import numpy as np
from support import OFF, dataset

from dataclasses import replace
import timetag_explorer.dataset as dataset_module
from timetag_explorer.frames import Corrections
from timetag_explorer.regions import Box
from timetag_explorer.vignette import MIN_SENSITIVITY, VIGNETTE_APERTURE, event_gain, sensitivity

# the bottom ramp: 1 at and above `start`, 1 - loss at `end`, the same slope below, never under MIN_SENSITIVITY
s = sensitivity(np.array([500.0, 128.0, 78.0, 28.0, 0.0]), start=128.0, end=28.0, loss=0.05)
assert np.allclose(s, [1.0, 1.0, 0.975, 0.95, 1 - 0.05 * 128 / 100]), s
assert sensitivity(-1e6, 128.0, 28.0, 0.05) == MIN_SENSITIVITY
# the top ramp is the same thing the other way up: start below end, the loss grows with the row
s = sensitivity(np.array([0.0, 1900.0, 1950.0, 2000.0, 2047.0]), start=1900.0, end=2000.0, loss=0.10)
assert np.allclose(s, [1.0, 1.0, 0.95, 0.90, 1 - 0.10 * 147 / 100]), s
assert sensitivity(1e6, 1900.0, 2000.0, 0.10) == MIN_SENSITIVITY
try:
    sensitivity(0.0, 28.0, 28.0, 0.05)
    raise SystemExit('a ramp with no length must be refused')
except ValueError:
    pass
# events carry 1-indexed AXIS2; the ramps are in 0-indexed rows (the image frame's y axis); both multiply
bottom, top = (128.0, 28.0, 0.05), (1900.0, 2000.0, 0.10)
x = np.full(3, 1000)
assert np.allclose(event_gain(x, [29, 1000, 2001], bottom, top), [1 / 0.95, 1.0, 1 / 0.90])
assert event_gain(x, [29, 129, 1000], bottom, top).dtype == np.float32
# the geometric distortion moves the row the ramps are read at (by the table's offset at the recorded pixel)
lut = np.zeros((2048, 2048), np.float32)
lut[28, :] = 50.0                                   # recorded row 28 (AXIS2 = 29) is really at row 78
assert np.allclose(event_gain(x[:2], [29, 30], bottom, top, lut), [1 / 0.975, 1 / (1 - 0.05 * 99 / 100)])

# the defaults are the measured ramps; the model is off
c = Corrections()
assert not c.vigcorr and c.vignette_ramps() is None
assert (c.vigbot_start, c.vigbot_end, c.vigbot_loss) == (430.0, 0.0, 0.25)
assert (c.vigtop_start, c.vigtop_end, c.vigtop_loss) == (1400.0, 2047.0, 0.23)
assert replace(c, vigcorr=True).vignette_ramps() == ((430.0, 0.0, 0.25), (1400.0, 2047.0, 0.23))
assert VIGNETTE_APERTURE == 'F25SRF2'

# counts: ten events in the row of the bottom ramp's end, ten in the top ramp's, ten in the middle
n = 30
t = np.linspace(0, 99, n)
xs = np.full(n, 1000)
y = np.concatenate([np.full(10, 29), np.full(10, 1000), np.full(10, 2001)])
ramps = dict(vigbot_start=128.0, vigbot_end=28.0, vigbot_loss=0.05, vigtop_start=1900.0, vigtop_end=2000.0, vigtop_loss=0.10)
on = replace(OFF, **ramps, vigcorr=True)
want = 10 / 0.95 + 10 / 0.90 + 10

ds = dataset(t, xs, y)                                    # no flat field: weights come from the ramps alone
ds.aperture = 'F25SRF2'
assert ds.vignette_applies()
total = lambda d, corr: d.image(None, 'image', corr).sum()
assert abs(total(ds, OFF) - 30 / 100.0) < 1e-9
assert abs(total(ds, on) - want / 100.0) < 1e-7           # (the gain is stored as float32)
assert total(ds, replace(OFF, **ramps)) == total(ds, OFF)  # switch off: the numbers do nothing
# the light curve (running-sum path, and the region path) is the same weighted count
for region in (None, Box(900, 1100, 0, 2047)):
    edges, rate = ds.light_curve(100.0, region, 'image', on)
    assert abs(np.nansum(rate * 100.0) - want) < 1e-5, (region, rate)
# changing the numbers, or switching the model, is not answered from the previous running sum
for corr, w in {on: want, replace(on, vigtop_loss=0.20): 10 / 0.95 + 10 / 0.80 + 10, OFF: 30.0}.items():
    edges, rate = ds.light_curve(100.0, None, 'image', corr)
    assert abs(np.nansum(rate) * 100.0 - w) < 1e-5, (corr, w)
# with a flat field the weights multiply
weighted = dataset(t, xs, y, weight=np.full(n, 2.0))
weighted.aperture = 'F25SRF2'
assert abs(weighted.image(None, 'image', replace(on, flatcorr=True)).sum() - 2 * want / 100.0) < 1e-6
edges, rate = weighted.light_curve(100.0, None, 'image', replace(on, flatcorr=True))
assert abs(np.nansum(rate) * 100.0 - 2 * want) < 1e-5

# only F25SRF2: any other aperture (or none recorded) is left exactly as it was
for other in ('', '52X0.2', '52X2', 'F28X50LP', 'f25srf2 '):
    other_ds = dataset(t, xs, y)
    other_ds.aperture = other.strip()
    if other.strip().upper() == 'F25SRF2':
        assert other_ds.vignette_applies()           # (the header value is compared without regard to case)
        continue
    assert not other_ds.vignette_applies(), other
    assert total(other_ds, on) == total(other_ds, OFF), other
    edges, rate = other_ds.light_curve(100.0, None, 'image', on)
    assert abs(np.nansum(rate) * 100.0 - 30.0) < 1e-6, other

# the geometric step comes first: with the correction on, the rows are those the distortion moved the events to
real_lut = dataset_module.load_geo_lut
try:
    shifted = (np.zeros((2048, 2048), np.float32), lut)   # (dx, dy): rows recorded at 28 land on 78
    dataset_module.load_geo_lut = lambda: shifted
    geo_on = replace(on, geocorr=True)
    edges, rate = ds.light_curve(100.0, None, 'image', geo_on)
    assert abs(np.nansum(rate) * 100.0 - (10 / 0.975 + 10 / 0.90 + 10)) < 1e-5, rate
    edges, rate = ds.light_curve(100.0, None, 'image', on)           # geometric step off: the recorded rows
    assert abs(np.nansum(rate) * 100.0 - want) < 1e-5, rate
    assert abs(ds.image(None, 'image', geo_on).sum() * 100.0 - (10 / 0.975 + 10 / 0.90 + 10)) < 1e-4
finally:
    dataset_module.load_geo_lut = real_lut
print('vignette ok')
