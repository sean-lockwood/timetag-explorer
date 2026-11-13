"""STIS FUV-MAMA geometric distortion correction (the pipeline's GEOCORR), from the instrument's IDC table.

The reference file is `1771408ro_idc.fits` (the STIS FUV-MAMA IDCTAB in use since 2017-07-07; Sohn 2018, STIS ISR
2018-02), a public STScI file that ships with the package (`timetag_explorer/data`). To use another IDCTAB, point
`config.IDC_FILE` (environment variable `TAG_EXPLORER_IDC`) at it.

What the table holds. Each FORWARD row (observed -> true positions) has, per optical element, a fourth-order polynomial

    XT = sum_ij CXij * X^j * Y^(i-j)          YT = sum_ij CYij * X^j * Y^(i-j)          (0 <= j <= i <= 4)

where (X, Y) is the observed position in LOW-RES pixels relative to the detector reference point (XREF, YREF) = (512, 512),
counted from 1 as in FITS, and (XT, YT) is the true position in arcseconds from that point. Dividing by the row's SCALE
(arcsec per pixel) gives the true position in pixels of a rectified image, whose pixels are square. That is what the
pipeline's x2d image has, and it includes the small change of scale that makes the FUV-MAMA's slightly rectangular pixels
square, not just the higher-order distortion. The table's INVERSE rows give the reverse mapping (true -> observed); they
are not used here, except to check the forward rows in the tests.

What this module does with it. TIME-TAG events are counted in hi-res pixels (2048 x 2048, two to a low-res pixel), so a
hi-res position is first expressed in the low-res frame, the polynomial is evaluated there, and the offset (true minus
observed) is doubled back to hi-res pixels. `load_geo_lut` returns that offset for every detector pixel: the
polynomial is smooth (the offset changes by less than 0.01 pixel across a pixel), so one table lookup per event gives the
same answer as evaluating it per event.

Choices worth knowing:
  * The polynomial's reference is XREF = YREF = 512 in a frame counted from 1, so hi-res pixel i (counted from 0) is at
    i / 2 + 0.75 (a low-res pixel is two hi-res pixels, and pixel 0 is centred a quarter of a low-res pixel below
    the first low-res centre). Moving the reference by a whole low-res pixel changes the offset by less than
    0.01 pixel (tested), so the exact half-pixel convention does not matter at the table's accuracy (RMS 4 mas, 0.16 pixel).
  * The two optical elements' rows have identical coefficients but different SCALE (MIRCUV 0.024395, MIRFUV 0.024468
    arcsec per pixel, 0.3% apart: about 1.5 low-res pixels at the detector's edge), so the row is chosen by the
    exposure's OPT_ELEM; an exposure whose OPT_ELEM the table does not list uses the MIRFUV row. The mean pixel scale the
    reference frames use (`frames.ARCSEC_PER_LOWRES_PIXEL`) is the MIRFUV one.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from . import config

BUNDLED_IDC = Path(__file__).resolve().parent / 'data' / '1771408ro_idc.fits'
DEFAULT_OPT_ELEM = 'MIRFUV'     # the row used when the table has none for an exposure's optical element
LOWRES_FROM_HIRES_OFFSET = 0.75 # low-res pixel counted from 1 = hi-res pixel counted from 0, over 2, plus this
DETECTOR_PIXELS = 2048


def idc_file() -> Path:
    """The IDC table in use: `config.IDC_FILE` if set, otherwise the copy that ships with the package."""
    return Path(config.IDC_FILE) if config.IDC_FILE else BUNDLED_IDC


@dataclass(frozen=True)
class IdcRow:
    """One FORWARD row of the IDC table: the polynomials, as coefficient tuples `cx[i][j]` (i = 0..norder, j = 0..i)."""
    opt_elem: str
    filter: str
    scale: float        # arcsec per pixel of the rectified image
    xref: float         # reference point of the observed positions, low-res pixels counted from 1
    yref: float
    cx: tuple
    cy: tuple

    @property
    def norder(self) -> int:
        return len(self.cx) - 1


def read_idc(path) -> tuple:
    """The FORWARD rows of an IDC table (a FITS binary table in extension 1), as `IdcRow`s."""
    from astropy.io import fits
    with fits.open(path) as hdul:
        table = hdul[1]
        norder = int(table.header.get('NORDER', 4))
        data = table.data
        names = set(data.columns.names)
        rows = []
        for r in range(len(data)):
            if str(data['DIRECTION'][r]).strip().upper() != 'FORWARD':
                continue

            def coefficients(axis):   # a missing column is a zero coefficient
                return tuple(tuple(float(data[f'C{axis}{i}{j}'][r]) if f'C{axis}{i}{j}' in names else 0.0
                                   for j in range(i + 1)) for i in range(norder + 1))
            rows.append(IdcRow(opt_elem=str(data['OPT_ELEM'][r]).strip().upper(), filter=str(data['FILTER'][r]).strip(),
                               scale=float(data['SCALE'][r]), xref=float(data['XREF'][r]), yref=float(data['YREF'][r]),
                               cx=coefficients('X'), cy=coefficients('Y')))
    if not rows:
        raise ValueError(f'{path} has no FORWARD rows: it is not a STIS IDC table')
    return tuple(rows)


@lru_cache(maxsize=2)
def _rows(path: str, modified: int) -> tuple:
    return read_idc(path)       # (modified: a replaced file is read again)


def idc_rows() -> tuple:
    path = idc_file()
    if not path.is_file():
        raise FileNotFoundError(f'IDC table not found: {path}')
    return _rows(str(path), path.stat().st_mtime_ns)


def select_row(rows, opt_elem: str = '') -> IdcRow:
    """The row for an optical element: the one named, else the MIRFUV row, else the first row whose FILTER is ANY."""
    for wanted in (str(opt_elem or '').strip().upper(), DEFAULT_OPT_ELEM):
        for row in rows:
            if row.opt_elem == wanted:
                return row
    for row in rows:
        if row.filter.upper() == 'ANY':
            return row
    raise ValueError(f'the IDC table has no FORWARD row for {opt_elem!r}, {DEFAULT_OPT_ELEM}, or FILTER = ANY')


def polynomial(coeff, x, y):
    """sum_ij coeff[i][j] * x^j * y^(i-j); x and y may be arrays that broadcast against each other."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    order = len(coeff) - 1
    xp = [x ** k for k in range(order + 1)]
    yp = [y ** k for k in range(order + 1)]
    total = 0.0
    for i in range(order + 1):
        for j in range(i + 1):
            if coeff[i][j]:
                total = total + coeff[i][j] * xp[j] * yp[i - j]
    return total


def forward_offsets(row: IdcRow, xr, yr):
    """The correction, in hi-res pixels, that takes an observed position to its true position: (dx, dy) to add to
    (xr, yr), the hi-res pixel position counted from 0 (arrays that broadcast against each other are fine)."""
    xt = np.asarray(xr, dtype=np.float64) / 2.0 + LOWRES_FROM_HIRES_OFFSET - row.xref    # low-res pixels from the reference
    yt = np.asarray(yr, dtype=np.float64) / 2.0 + LOWRES_FROM_HIRES_OFFSET - row.yref
    px = polynomial(row.cx, xt, yt) / row.scale      # true position, in pixels from the reference
    py = polynomial(row.cy, xt, yt) / row.scale
    return 2.0 * (px - xt), 2.0 * (py - yt)


def gdcorr(xr, yr, row: IdcRow | None = None):
    """Distortion-corrected position of the observed hi-res pixel position (xr, yr), counted from 0: (xc, yc).
    `row` is the IDC row to use (default: the MIRFUV row of the table in use)."""
    row = select_row(idc_rows()) if row is None else row
    dx, dy = forward_offsets(row, xr, yr)
    return np.asarray(xr, dtype=np.float64) + dx, np.asarray(yr, dtype=np.float64) + dy


@lru_cache(maxsize=4)
def _lut(row: IdcRow):
    axis = np.arange(DETECTOR_PIXELS, dtype=np.float64)
    dx, dy = forward_offsets(row, axis[None, :], axis[:, None])      # [row = y, column = x]
    lut = (np.ascontiguousarray(dx, dtype=np.float32), np.ascontiguousarray(dy, dtype=np.float32))
    for arr in lut:
        arr.flags.writeable = False
    return lut


def load_geo_lut(opt_elem: str = ''):
    """Distortion offsets (dx, dy) in hi-res pixels at every integer detector pixel (0-indexed, [row = y, column = x]),
    for an exposure with this OPT_ELEM. Built once per optical element and process (about a second)."""
    return _lut(select_row(idc_rows(), opt_elem))
