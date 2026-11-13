"""Shared by the tests: make the package importable from a checkout (no install needed) and build small datasets."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from timetag_explorer.dataset import TagDataset  # noqa: E402
from timetag_explorer.frames import Corrections, FrameInfo  # noqa: E402


# ---- reading the IDC table where astropy is only a stand-in (the sandbox this was developed in has none) ------------------
# The package reads the table with astropy. This reads the same file with the struct module, so that the tests that turn the
# geometric distortion on run without astropy; with the real astropy, the package's own reader is used and
# test_gdcorr.py checks it against this one.

def read_idc_without_astropy(path):
    """The FORWARD rows of an IDC table (a FITS file: primary header, then one binary table), read byte by byte."""
    import struct
    from timetag_explorer.gdcorr import IdcRow
    raw = Path(path).read_bytes()

    def header_end(start):
        pos = start
        while True:
            block = raw[pos:pos + 2880]
            if any(block[i:i + 8] == b'END     ' for i in range(0, 2880, 80)):
                return pos + 2880
            pos += 2880
    table_header = header_end(0)
    data = header_end(table_header)
    cards = [raw[i:i + 80].decode('ascii') for i in range(table_header, data, 80)]
    value = lambda key: next(c[10:].split('/')[0].strip() for c in cards if c[:8].strip() == key)
    width, nrows, nfields = int(value('NAXIS1')), int(value('NAXIS2')), int(value('TFIELDS'))
    names = [value(f'TTYPE{n}').strip("' ") for n in range(1, nfields + 1)]
    forms = [value(f'TFORM{n}').strip("' ") for n in range(1, nfields + 1)]
    layout, offset = {}, 0
    for name, form in zip(names, forms):
        kind = form[-1]
        count = int(form[:-1] or 1)
        size = count * {'A': 1, 'D': 8, 'J': 4}[kind]
        layout[name] = (offset, kind, count)
        offset += size
    assert offset == width, (offset, width)
    norder = int(value('NORDER'))
    rows = []
    for r in range(nrows):
        base = data + r * width

        def field(name):
            start, kind, count = layout[name]
            chunk = raw[base + start: base + start + {'A': count, 'D': 8, 'J': 4}[kind]]
            return chunk.decode('ascii').strip() if kind == 'A' else struct.unpack('>d' if kind == 'D' else '>i', chunk)[0]
        if field('DIRECTION').upper() != 'FORWARD':
            continue
        coefficients = lambda axis: tuple(tuple(field(f'C{axis}{i}{j}') if f'C{axis}{i}{j}' in layout else 0.0
                                                for j in range(i + 1)) for i in range(norder + 1))
        rows.append(IdcRow(opt_elem=field('OPT_ELEM').upper(), filter=field('FILTER'), scale=field('SCALE'),
                           xref=field('XREF'), yref=field('YREF'), cx=coefficients('X'), cy=coefficients('Y')))
    return tuple(rows)


def _astropy_is_a_stand_in() -> bool:
    from unittest.mock import MagicMock
    try:
        from astropy.io import fits
    except Exception:
        return True
    return isinstance(fits, MagicMock) or isinstance(getattr(fits, 'open', None), MagicMock)


ASTROPY_IS_REAL = not _astropy_is_a_stand_in()
if not ASTROPY_IS_REAL:
    from timetag_explorer import gdcorr as _gdcorr    # noqa: E402
    _gdcorr.read_idc = read_idc_without_astropy

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
