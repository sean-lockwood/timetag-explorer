"""The selection-script generator: the script it writes is valid Python, carries its version, and says what it does.
(Running the scripts on real data is done by validation/check_selection_scripts.py.)"""
import support  # noqa: F401
import numpy as np
from support import OFF, dataset

import timetag_explorer as tte
from timetag_explorer.frames import Corrections
from timetag_explorer.regions import Box, Circle, Polygon

n = 200
ds = dataset(np.linspace(0, 99, n), np.full(n, 600), np.full(n, 700))
cases = [((10.0, 20.0), None, 'image', OFF),
         (None, Box(500, 700, 600, 800), 'image', OFF),
         ((10.0, 20.0), Circle(600, 700, 50), 'rest', Corrections(randcorr=True, geocorr=False, timecorr=True)),
         (None, Polygon(((500, 500), (700, 500), (600, 700))), 'image', Corrections(randcorr=False, geocorr=True))]
for time_sel, region, frame, corr in cases:
    try:
        s = tte.build_selection_script(ds, time_sel, region, frame, corr)
    except ValueError as exc:       # a rest frame needs the SPT file's pointing information
        assert frame != 'image', exc
        continue
    assert not s.warnings or frame != 'image', s.warnings
    compile(s.code, s.filename, 'exec')                                  # valid Python
    assert tte.SELECTION_SCRIPT_VERSION in s.code and 'SELVER' in s.code
    assert s.filename.endswith('.py') and ds.rootname in s.filename
    assert isinstance(s.packages, (list, tuple)) and s.notes
    assert f"SCRIPT_VERSION = '{tte.SELECTION_SCRIPT_VERSION}'" in s.code
    if time_sel is not None:
        assert f'TIME_WINDOW = ({time_sel[0]}, {time_sel[1]})' in s.code
    if isinstance(region, Box):
        assert f'X0, X1 = {float(region.x0)}, {float(region.x1)}' in s.code
    if region is not None:
        assert region.describe() in s.code        # the region is written out in words, too
assert tte.SELECTION_SCRIPT_VERSION == '1.2'

# The script must carry each column's coordinate keywords (TCTYPn, TCRPXn ...) into the file it writes: astropy
# writes those from attributes of the column, so a rebuilt column loses them unless they are passed on. (Run just
# the helper out of a generated script, with a stand-in for astropy's Column that records what it is given.)
import types
s = tte.build_selection_script(ds, (10.0, 20.0), Box(400, 600, 400, 600), 'image', Corrections(randcorr=False, geocorr=False))
start, stop = s.code.index('STRUCTURE_KEYWORD ='), s.code.index('def in_file_order')
made = []
fits_stub = types.SimpleNamespace(Column=lambda **kw: made.append(kw) or kw)
namespace = {'np': np, 're': __import__('re'), 'fits': fits_stub}
exec(s.code[start:stop], namespace)
column = types.SimpleNamespace(unit='pixel', coord_type='RA---TAN', coord_unit='deg', coord_ref_point=512.5,
                               coord_ref_value=150.25, coord_inc=-1e-6, time_ref_pos=None)
hdu = types.SimpleNamespace(columns={'AXIS1': column}, data={'AXIS1': np.arange(3, dtype=np.int16)})
namespace['stored_column'](hdu, 'AXIS1')
assert made[0]['coord_ref_point'] == 512.5 and made[0]['coord_type'] == 'RA---TAN' and made[0]['coord_inc'] == -1e-6, made
assert made[0]['format'] == 'I' and made[0]['unit'] == 'pixel'
# EDGE_OK in the script is the Explorer's own edge trim (30, 6, 28, 6 hi-res pixels: left, right, bottom, top), event by
# event: run the script's two lines of constants and its test on a row of events along each axis.
from timetag_explorer.constants import EDGE_TRIM
assert tuple(EDGE_TRIM) == (30, 6, 28, 6)
axis = np.concatenate([np.arange(1, 70), np.arange(1980, 2049)])         # both ends of an axis
xs, ys = np.concatenate([axis, np.full(axis.size, 1000)]), np.concatenate([np.full(axis.size, 1000), axis])
trimmed = dataset(np.linspace(0, 99, xs.size), xs, ys)
on = Corrections(randcorr=False, geocorr=False, flatcorr=False, timecorr=True)
s = tte.build_selection_script(trimmed, (10.0, 20.0), Box(400, 600, 400, 600), 'image', on)
names = {}
exec('\n'.join(line for line in s.code.splitlines() if line.startswith('EDGE_')), names)
assert (names['EDGE_X_LOW'], names['EDGE_X_HIGH'], names['EDGE_Y_LOW'], names['EDGE_Y_HIGH']) == (31, 2042, 29, 2042)
x, y = trimmed.x, trimmed.y                                  # (time-ordered; the dataset sorts by time, which is already so)
in_script = (x >= names['EDGE_X_LOW']) & (x <= names['EDGE_X_HIGH']) & (y >= names['EDGE_Y_LOW']) & (y <= names['EDGE_Y_HIGH'])
assert '(x >= EDGE_X_LOW) & (x <= EDGE_X_HIGH) & (y >= EDGE_Y_LOW) & (y <= EDGE_Y_HIGH)' in s.code
assert np.array_equal(in_script, trimmed.valid_events(on)) and not in_script.all() and in_script.any()
for key, pixels in zip(('SELEDGEL', 'SELEDGER', 'SELEDGEB', 'SELEDGET'), EDGE_TRIM):
    assert f"header['{key}'] = ({pixels}," in s.code, key
assert 'left 30, right 6, bottom 28, top 6' in s.code
print('selection script ok')
