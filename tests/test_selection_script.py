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
assert tte.SELECTION_SCRIPT_VERSION == '1.0'
print('selection script ok')
