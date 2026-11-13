"""The geometric-distortion correction (GEOCORR) from the STIS FUV-MAMA IDC table: the table is read correctly, the
polynomial is applied as the table's header and STIS ISR 2018-02 define it, the forward and inverse rows agree, and the
hi-res lookup table follows the exposure's optical element. The IDC table is the bundled copy of 1771408ro_idc.fits."""
import support  # noqa: F401
import numpy as np

from timetag_explorer import config, gdcorr
from timetag_explorer.frames import ARCSEC_PER_LOWRES_PIXEL, DEFAULT_CORRECTIONS, geocorr_available

PATH = gdcorr.BUNDLED_IDC
assert PATH.is_file() and gdcorr.idc_file() == PATH and geocorr_available() and DEFAULT_CORRECTIONS.geocorr
rows = gdcorr.idc_rows()

# ---- the table: two FORWARD rows (MIRCUV, MIRFUV) with the same coefficients and different SCALE ---------------------------
assert [r.opt_elem for r in rows] == ['MIRCUV', 'MIRFUV'] and all(r.norder == 4 for r in rows)
cuv, fuv = rows
assert (cuv.scale, fuv.scale) == (0.024395, 0.024468) and fuv.scale == ARCSEC_PER_LOWRES_PIXEL
assert cuv.cx == fuv.cx and cuv.cy == fuv.cy and (cuv.xref, cuv.yref) == (512.0, 512.0)
assert cuv.cx[0] == (0.0,) and cuv.cy[0] == (0.0,)          # the reference point maps to itself

# ---- the coefficients are those of ISR 2018-02, Table 3 (forward: cx, cy), to the 10 digits it prints --------------------------
ISR_FORWARD = {   # (i, j): (cx, cy)
    (1, 0): (-2.346813458e-05, 2.451026177e-02), (1, 1): (2.427007264e-02, -4.484187512e-05),
    (2, 0): (1.584981859e-08, -5.427108105e-08), (2, 1): (1.238586312e-07, 5.584123466e-09),
    (2, 2): (-5.686955680e-09, 1.845072272e-07), (3, 0): (9.821534833e-12, -4.343839277e-11),
    (3, 1): (-1.120214702e-10, 2.636735916e-11), (3, 2): (3.464251745e-11, -1.567986544e-12),
    (3, 3): (3.020217823e-11, 1.396986362e-11), (4, 0): (1.819768264e-13, -7.639717075e-14),
    (4, 1): (2.268814640e-13, 1.277725915e-13), (4, 2): (7.586840087e-14, 1.167041365e-13),
    (4, 3): (6.678498410e-14, -9.002936012e-14), (4, 4): (1.685402081e-14, -9.190396791e-13)}
for (i, j), (isr_x, isr_y) in ISR_FORWARD.items():
    assert abs(fuv.cx[i][j] - isr_x) <= 1e-9 * abs(isr_x), (i, j)
    assert abs(fuv.cy[i][j] - isr_y) <= 1e-9 * abs(isr_y), (i, j)

# ---- the polynomial: each coefficient multiplies one monomial, j = the power of x ------------------------------------------
x, y = 3.0, 5.0
for i in range(5):
    for j in range(i + 1):
        c = tuple(tuple(1.0 if (a, b) == (i, j) else 0.0 for b in range(a + 1)) for a in range(5))
        assert gdcorr.polynomial(c, x, y) == x ** j * y ** (i - j), (i, j)
assert np.array_equal(gdcorr.polynomial(fuv.cx, np.array([[1.0, 2.0]]), np.array([[3.0], [4.0]])).shape, (2, 2))   # broadcasts

# ---- the correction, against a direct evaluation of the formulas in the table's header ------------------------------------
def direct(row, xr, yr):
    """XT = CX10*Y + CX11*X + ... written out term by term, in arcsec; then to pixels and back to hi-res pixels."""
    X, Y = xr / 2.0 + 0.75 - 512.0, yr / 2.0 + 0.75 - 512.0
    out = []
    for c in (row.cx, row.cy):
        t = (c[1][0] * Y + c[1][1] * X + c[2][0] * Y ** 2 + c[2][1] * X * Y + c[2][2] * X ** 2
             + c[3][0] * Y ** 3 + c[3][1] * X * Y ** 2 + c[3][2] * X ** 2 * Y + c[3][3] * X ** 3
             + c[4][0] * Y ** 4 + c[4][1] * X * Y ** 3 + c[4][2] * X ** 2 * Y ** 2 + c[4][3] * X ** 3 * Y + c[4][4] * X ** 4)
        out.append(t / row.scale)
    return 2.0 * (out[0] - X), 2.0 * (out[1] - Y)
for xr, yr in ((0.0, 0.0), (2047.0, 2047.0), (1023.5, 1023.5), (300.0, 1800.0), (1700.0, 90.0)):
    for row in rows:
        dx, dy = gdcorr.forward_offsets(row, xr, yr)
        ex, ey = direct(row, xr, yr)
        assert abs(dx - ex) < 1e-9 and abs(dy - ey) < 1e-9, (xr, yr, dx, ex)
        xc, yc = gdcorr.gdcorr(xr, yr, row)
        assert abs(xc - (xr + dx)) < 1e-12 and abs(yc - (yr + dy)) < 1e-12

# sizes: the correction is the change of scale that makes the FUV-MAMA's slightly rectangular pixels square (a linear term:
# -0.8% in x for MIRFUV) plus the distortion proper, which ISR 2018-02 (Fig. 4) puts at about 4 low-res pixels at most
grid = np.arange(0.0, 2048.0, 16.0)
gx, gy = np.meshgrid(grid, grid)
for row in rows:
    dx, dy = gdcorr.forward_offsets(row, gx, gy)
    X, Y = gx / 2.0 + 0.75 - 512.0, gy / 2.0 + 0.75 - 512.0
    linear_x = 2.0 * ((row.cx[1][1] * X + row.cx[1][0] * Y) / row.scale - X)       # only the terms of the first order
    linear_y = 2.0 * ((row.cy[1][1] * X + row.cy[1][0] * Y) / row.scale - Y)
    assert 2.0 < np.abs(dx - linear_x).max() < 8.0 and 2.0 < np.abs(dy - linear_y).max() < 8.0     # hi-res pixels (1 to 4 low-res)
    assert np.abs(dx).max() < 20.0 and np.abs(dy).max() < 20.0
corner = [gdcorr.forward_offsets(row, 2047.0, 2047.0)[0] for row in (cuv, fuv)]
assert abs(corner[0] - corner[1]) > 2.0, 'the two optical elements use different SCALE'    # 0.3% of ~1000 hi-res pixels

# ---- the inverse rows: the true -> observed mapping undoes the forward one --------------------------------------------------
inverse = {}
if support.ASTROPY_IS_REAL:
    from astropy.io import fits
    with fits.open(PATH) as hdul:
        table = hdul[1].data
        for r, direction in enumerate(table['DIRECTION']):
            if str(direction).strip() == 'INVERSE' and str(table['OPT_ELEM'][r]).strip() == 'MIRFUV':
                inverse = {a: tuple(tuple(float(table[f'C{a}{i}{j}'][r]) for j in range(i + 1)) for i in range(5)) for a in 'XY'}
else:   # (the same rows, read by the helper that stands in for astropy)
    import struct
    raw = PATH.read_bytes()
    start = raw.index(b'INVERSE ANY     MIRFUV')
    values = struct.unpack('>d2i8d2i2d30d', raw[start + 42: start + 42 + 344])[15:]
    names = [f'{k}{i}{j}' for k in 'XY' for i in range(5) for j in range(i + 1)]
    table = dict(zip(names, values))
    inverse = {a: tuple(tuple(table[f'{a}{i}{j}'] for j in range(i + 1)) for i in range(5)) for a in 'XY'}
grid = np.arange(0.0, 2048.0, 64.0)
gx, gy = np.meshgrid(grid, grid)
X, Y = gx / 2.0 + 0.75 - 512.0, gy / 2.0 + 0.75 - 512.0                  # observed, low-res pixels from the reference
xt = gdcorr.polynomial(fuv.cx, X, Y)                                      # true, arcsec
yt = gdcorr.polynomial(fuv.cy, X, Y)
# The inverse rows store the terms of order n divided by 41**(n-1) (41 is about 1/SCALE, the arcseconds-per-pixel factor), so they
# are evaluated here on the true position in units of 1/41 arcsec.  Read literally, in arcseconds, they would not undo the forward
# polynomial (3 low-res pixels off at the corners); read this way they do, to 0.02 pixel, which also checks the forward rows.
# This package only uses the FORWARD rows (observed -> true), so the inverse rows matter here only for that cross-check.
UNIT = 41.0
back_x = gdcorr.polynomial(inverse['X'], xt * UNIT, yt * UNIT) / UNIT      # observed again, low-res pixels
back_y = gdcorr.polynomial(inverse['Y'], xt * UNIT, yt * UNIT) / UNIT
assert np.abs(back_x - X).max() < 0.05 and np.abs(back_y - Y).max() < 0.05, (np.abs(back_x - X).max(), np.abs(back_y - Y).max())
literal_x = gdcorr.polynomial(inverse['X'], xt, yt)
assert np.abs(literal_x - X).max() > 1.0, 'the inverse rows are not meant to be read in arcseconds'

# ---- the reference convention does not matter at the table's accuracy ---------------------------------------------------------
from dataclasses import replace
moved = replace(fuv, xref=fuv.xref + 1.0, yref=fuv.yref + 1.0)              # the reference one whole low-res pixel away
ref = np.array(gdcorr.forward_offsets(fuv, gx, gy)) / 2.0
mov = np.array(gdcorr.forward_offsets(moved, gx, gy)) / 2.0
assert np.abs(ref - mov).max() < 0.05, np.abs(ref - mov).max()               # low-res pixels; the table's RMS is 0.16

# ---- the lookup table: hi-res offsets at every detector pixel, per optical element --------------------------------------------
lut_x, lut_y = gdcorr.load_geo_lut('MIRFUV')
assert lut_x.shape == lut_y.shape == (2048, 2048) and lut_x.dtype == np.float32 and not lut_x.flags.writeable
for r_, c_ in ((0, 0), (2047, 2047), (1023, 1024), (100, 1900)):                    # [row = y, column = x]
    dx, dy = gdcorr.forward_offsets(fuv, float(c_), float(r_))
    assert abs(lut_x[r_, c_] - dx) < 1e-4 and abs(lut_y[r_, c_] - dy) < 1e-4
assert gdcorr.load_geo_lut('mirfuv')[0] is lut_x and gdcorr.load_geo_lut('')[0] is lut_x      # (case-insensitive; the default row)
assert gdcorr.load_geo_lut('G140L')[0] is lut_x                                    # an optical element the table lacks: MIRFUV
assert gdcorr.load_geo_lut('MIRCUV')[0] is not lut_x and not np.array_equal(gdcorr.load_geo_lut('MIRCUV')[0], lut_x)
assert abs(float(lut_x[1023:1025, 1023:1025].mean())) < 0.1                         # near the reference point: almost nothing
assert gdcorr.select_row(rows, 'mircuv') is cuv and gdcorr.select_row(rows, '') is fuv

# ---- another table, and a table that is not an IDC one -------------------------------------------------------------------------
import tempfile
import pathlib
with tempfile.TemporaryDirectory() as folder:
    missing = pathlib.Path(folder) / 'none.fits'
    config.IDC_FILE = missing
    try:
        assert not geocorr_available()
        try:
            gdcorr.idc_rows()
            raise SystemExit('a missing IDC table must be an error')
        except FileNotFoundError as err:
            assert 'none.fits' in str(err)
    finally:
        config.IDC_FILE = None
assert gdcorr.idc_file() == PATH and geocorr_available()

# ---- the package's own reader (astropy) agrees with the byte-level one, when astropy is available -----------------------------
if support.ASTROPY_IS_REAL:
    assert gdcorr.read_idc(PATH) == support.read_idc_without_astropy(PATH)
print('gdcorr ok')
