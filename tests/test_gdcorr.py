"""The geometric-distortion polynomial, and the placeholder state of this release."""
import support  # noqa: F401
import numpy as np

from timetag_explorer import gdcorr
from timetag_explorer.frames import DEFAULT_CORRECTIONS, geocorr_available

# The polynomial itself: each coefficient multiplies one monomial.
x, y = 3.0, 5.0
monomials = [y, x, y * y, x * y, x * x, y ** 3, y * y * x, y * x * x, x ** 3, y ** 4, y ** 3 * x, y * y * x * x,
             y * x ** 3, x ** 4]
for k, expected in enumerate(monomials):
    c = [0.0] * 14
    c[k] = 1.0
    assert gdcorr.polyd4(x, y, c) == expected, k

# Placeholder coefficients: nothing is corrected, and the Explorer knows it.
assert gdcorr.PLACEHOLDER == (not (any(gdcorr.a) or any(gdcorr.b)))
if gdcorr.PLACEHOLDER:
    gx, gy = np.meshgrid(np.arange(0, 2048, 64, dtype=float), np.arange(0, 2048, 64, dtype=float))
    xc, yc = gdcorr.gdcorr(gx, gy)
    assert np.array_equal(xc, gx) and np.array_equal(yc, gy)
    assert not geocorr_available() and not DEFAULT_CORRECTIONS.geocorr
else:   # real coefficients: the switch is available and on by default
    assert geocorr_available()
print('gdcorr ok')
