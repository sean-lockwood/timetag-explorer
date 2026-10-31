"""The selection-script generator and the columns it adds (tagsel). The script is only the user's choices and a call
to tagsel.make_tagsel, so the script is checked by running it with make_tagsel replaced by a recorder; the columns are
checked on synthetic exposures; writing the FITS file is checked on real data by validation/check_selection_scripts.py."""
import support  # noqa: F401
import json
import types
from dataclasses import replace

import numpy as np
from support import OFF, dataset

import timetag_explorer as tte
from timetag_explorer import selection_script, tagsel
from timetag_explorer.constants import EDGE_TRIM, RANDOM_SEED
from timetag_explorer.frames import Corrections, FrameInfo
from timetag_explorer.regions import Box, Circle, Polygon
from timetag_explorer.transforms import frame_coordinates

assert tte.SELECTION_SCRIPT_VERSION == '2.1'
assert tuple(EDGE_TRIM) == (30, 6, 28, 6)

n = 200
ds = dataset(np.linspace(0, 99, n), np.full(n, 600), np.full(n, 700), weight=np.full(n, 1.25))
ds.aperture = 'F25SRF2'

# ---- the script: valid Python that passes the choices on ----------------------------------------------------------
calls = []
fake = types.ModuleType('timetag_explorer.tagsel')
fake.make_tagsel = lambda *a, **k: calls.append((a, k)) or 'out.fits'


def run_script(code, version=None):
    """Execute a generated script with make_tagsel recorded; returns its global names."""
    import sys
    real = sys.modules['timetag_explorer.tagsel']
    sys.modules['timetag_explorer.tagsel'] = fake
    saved = tte.__version__
    try:
        if version:
            tte.__version__ = version
        names = {'__name__': 'script'}
        exec(compile(code, 'script.py', 'exec'), names)
        return names
    finally:
        sys.modules['timetag_explorer.tagsel'] = real
        tte.__version__ = saved


corr_all = Corrections(randcorr=True, geocorr=False, flatcorr=True, timecorr=True)
cases = [((10.0, 20.0), None, 'image', OFF),
         (None, Box(500, 700, 600, 800), 'image', corr_all),
         ((10.0, 20.0), Circle(600, 700, 50), 'rest', corr_all),
         (None, Polygon(((500.0, 500.0), (700.0, 500.0), (600.0, 700.5))), 'rest_north', replace(corr_all, vigcorr=True))]
for time_sel, region, frame, corr in cases:
    s = tte.build_selection_script(ds, time_sel, region, frame, corr)
    assert not s.warnings, s.warnings
    assert s.filename == f'select_{ds.rootname}.py' and s.notes
    assert tte.SELECTION_SCRIPT_VERSION in s.code and f"SCRIPT_VERSION = '{tte.SELECTION_SCRIPT_VERSION}'" in s.code
    assert f"PACKAGE_VERSION = '{tte.__version__}'" in s.code and s.pip_command in s.code
    assert s.pip_command == f'pip install timetag-explorer=={tte.__version__}'      # (not an installed git copy here)
    assert len(s.code.splitlines()) < 170, 'the script should stay short'           # the work is in the package
    calls.clear()
    run_script(s.code)
    (args, kwargs), = calls
    assert args == (ds.rootname,)
    assert kwargs['time_sel'] == time_sel and kwargs['region'] == region
    assert kwargs['frame'] == (frame if region is not None else 'image')
    got = kwargs['corr']
    assert (got.randcorr, got.geocorr, got.flatcorr, got.timecorr, got.vigcorr) == (
        corr.randcorr, corr.geocorr, corr.flatcorr, corr.timecorr, corr.vigcorr and ds.vignette_applies())
    assert got.vignette_ramps() == (corr.vignette_ramps() if corr.vigcorr else None)
    assert kwargs['script_version'] == tte.SELECTION_SCRIPT_VERSION and kwargs['tool_url'] == tte.selection_script.tool_url()
    assert kwargs['overwrite'] is False and kwargs['north_up'] is False and 'INCLUDE_NORTH_UP = False' in s.code
    if region is not None:
        assert region.describe() in s.code        # the region is written out in words, too

# A region drawn in the app holds numpy floats; the script must not need numpy to run its REGION line
for drawn in (Box(*np.array([508.08, 819.27, 1192.6, 1452.37])), Circle(*np.array([600.5, 700.25, 50.125])),
              Polygon(tuple((np.float64(a), np.float64(b)) for a, b in ((500.5, 500.25), (700.0, 500.0), (600.0, 700.5))))):
    s = tte.build_selection_script(ds, None, drawn, 'image', OFF)
    line = next(l for l in s.code.splitlines() if l.startswith('REGION = '))
    assert 'np.' not in line and 'float64' not in s.code, line
    names = {}
    exec('from timetag_explorer.regions import Box, Circle, Polygon\n' + line, names)   # (no numpy imported)
    assert names['REGION'] == drawn and all(type(v) is float for v in vars(names['REGION']).values() if not isinstance(v, tuple))
    assert tagsel.region_source(drawn) in s.code

# vignetting is written only for an exposure it applies to
other = dataset(np.linspace(0, 99, n), np.full(n, 600), np.full(n, 700)); other.aperture = '52X2'
s = tte.build_selection_script(other, (1.0, 2.0), None, 'image', replace(OFF, vigcorr=True))
assert 'vigcorr' not in s.code and 'vignetting' not in s.code.split('"""')[1]

# a script made for one version of the package refuses to run with another (unless told to)
s = tte.build_selection_script(ds, (10.0, 20.0), None, 'image', OFF)
calls.clear()
try:
    run_script(s.code, version='0.0.1')
    raise SystemExit('a different package version must stop the script')
except SystemExit as exc:
    assert 'needs timetag-explorer' in str(exc) and s.pip_command.replace('pip install ', '') in str(exc), exc
assert not calls
calls.clear(); run_script(s.code.replace('ALLOW_OTHER_VERSION = False', 'ALLOW_OTHER_VERSION = True'), version='0.0.1')
assert len(calls) == 1

# a region drawn on a planet map (experimental) is not carried over, and the script says so
s = tte.build_selection_script(ds, (10.0, 20.0), Box(0, 10, 0, 10), 'map_rect', OFF)
assert len(s.warnings) == 1 and 'map' in s.warnings[0]
calls.clear(); run_script(s.code)
assert calls[0][1]['region'] is None and calls[0][1]['time_sel'] == (10.0, 20.0)

# the instruction to install the package: a release from PyPI, or exactly the branch / commit / archive that is running
class FakeDistribution:
    def __init__(self, direct_url): self.direct_url = direct_url
    def read_text(self, name): return self.direct_url


def spec_for(direct_url):
    real = selection_script.metadata.distribution
    selection_script.metadata.distribution = lambda name: FakeDistribution(direct_url)
    try:
        return selection_script.package_install_spec()
    finally:
        selection_script.metadata.distribution = real

v = tte.__version__
assert spec_for(None) == (f'timetag-explorer=={v}', None)
git = json.dumps({'url': 'https://github.com/sean-lockwood/timetag-explorer', 'vcs_info': {'vcs': 'git', 'commit_id': 'abc123', 'requested_revision': 'branch'}})
assert spec_for(git) == ('timetag-explorer @ git+https://github.com/sean-lockwood/timetag-explorer@abc123', None)
archive = json.dumps({'url': 'https://github.com/sean-lockwood/timetag-explorer/archive/c24f0f78.zip', 'archive_info': {'hash': 'sha256=0'}})
assert spec_for(archive) == ('timetag-explorer @ https://github.com/sean-lockwood/timetag-explorer/archive/c24f0f78.zip', None)
spec, warning = spec_for(json.dumps({'url': 'file:///home/me/timetag-explorer', 'dir_info': {'editable': True}}))
assert spec == f'timetag-explorer=={v}' and 'development copy' in warning
assert selection_script.pip_command('timetag-explorer @ https://x/y.zip') == 'pip install "timetag-explorer @ https://x/y.zip"'

# ---- the columns ---------------------------------------------------------------------------------------------------
# events along a row, some of them in the detector edge
m = 400
t = np.linspace(0, 99, m)
xs = np.linspace(5, 2040, m).astype(int)
ys = np.full(m, 1000)
w = np.linspace(0.5, 2.0, m)
slewing = FrameInfo(orientat=30.0, slew_angle=45.0, slew_rate=500.0)
exp = dataset(t, xs, ys, weight=w); exp.frame_info = slewing; exp.aperture = 'F25SRF2'
on = Corrections(randcorr=True, geocorr=False, flatcorr=True, timecorr=True)
cols = tagsel.selection_columns(exp, (20.0, 60.0), Box(900.0, 1300.0, 900.0, 1100.0), 'image', on)
a = cols.arrays
assert list(a) == ['SEL_TIME', 'SEL_REGION', 'EDGE_OK', 'XCORR', 'YCORR', 'XREST', 'YREST', 'FLATW'], list(a)   # (no north-up, no WEIGHT)
assert cols.flags == ('SEL_TIME', 'SEL_REGION', 'EDGE_OK')
assert [name for name, _ in cols.skipped] == ['Planet rest frame'], cols.skipped      # (not a moon observation)
# north-up columns only on request
full = tagsel.selection_columns(exp, (20.0, 60.0), Box(900.0, 1300.0, 900.0, 1100.0), 'image', on, north_up=True)
assert list(full.arrays) == ['SEL_TIME', 'SEL_REGION', 'EDGE_OK', 'XCORR', 'YCORR', 'XREST', 'YREST', 'XRESTN', 'YRESTN',
                             'FLATW'], list(full.arrays)
assert [name for name, _ in full.skipped] == ['Planet rest frame', 'Planet rest frame, north up'], full.skipped
assert np.array_equal(full.arrays['XREST'], a['XREST'])
assert all(len(v) == m for v in a.values())
# the time flag: the same events as the Explorer's time slice
assert np.array_equal(np.flatnonzero(a['SEL_TIME']), np.arange(*exp.time_slice((20.0, 60.0)).indices(m)))
# the corrected positions are the Explorer's own (jitter, same seed), counted from 1 like AXIS1 / AXIS2
u, v = frame_coordinates(exp.x, exp.y, exp.time, exp.frame_info, 'image', on, seed=RANDOM_SEED)
assert np.allclose(a['XCORR'], u + 1) and np.allclose(a['YCORR'], v + 1) and a['XCORR'].dtype == np.float32
assert np.all(np.abs(a['XCORR'] - exp.x) <= 0.5 + 1e-4) and np.any(a['XCORR'] != exp.x)     # (a half-pixel jitter)
# ... and the region is tested at those positions, event by event (not on bins)
assert np.array_equal(a['SEL_REGION'], Box(900.0, 1300.0, 900.0, 1100.0).contains(u, v)) and a['SEL_REGION'].any() and not a['SEL_REGION'].all()
# the other frames start from the corrected positions
ur, vr = frame_coordinates(exp.x, exp.y, exp.time, exp.frame_info, 'rest', on, seed=RANDOM_SEED)
assert np.allclose(a['XREST'], ur + 1) and np.allclose(a['YREST'], vr + 1)
assert not np.allclose(a['XREST'], a['XCORR']) and not np.allclose(full.arrays['XRESTN'], a['XREST'])      # shifted, then rotated
# edge trim and weights
assert np.array_equal(a['EDGE_OK'], exp.valid_events(on)) and not a['EDGE_OK'].all() and a['EDGE_OK'].any()
assert np.allclose(a['FLATW'], w) and 'WEIGHT' not in a

# a region drawn in the rest frame is tested in that frame
box = Box(900.0, 1300.0, 900.0, 1100.0)
cols = tagsel.selection_columns(exp, None, box, 'rest', on)
assert np.array_equal(cols.arrays['SEL_REGION'], box.contains(ur, vr)) and 'SEL_TIME' not in cols.arrays
# a region in a frame that cannot be built is an error, not a missing column
try:
    tagsel.selection_columns(exp, None, box, 'planet', on)
    raise SystemExit('a region in an unavailable frame must be refused')
except RuntimeError as err:
    assert 'cannot be built' in str(err)

# with every correction off there is no XCORR / YCORR, no weights, no edge flag; AXIS1 / AXIS2 are never rewritten
cols = tagsel.selection_columns(exp, (20.0, 60.0), None, 'image', OFF)
assert list(cols.arrays) == ['SEL_TIME', 'XREST', 'YREST'], list(cols.arrays)
# geometric distortion alone also gives the corrected pair (the jitter is off, so it differs from AXIS only by the table)
cols = tagsel.selection_columns(exp, None, None, 'image', replace(OFF, geocorr=True))
assert 'XCORR' in cols.arrays and np.allclose(cols.arrays['XCORR'], exp.x)                          # (zero offsets for now)
# a rest frame with no slew is the image frame
still = dataset(t, xs, ys); cols = tagsel.selection_columns(still, None, None, 'image', on)
assert np.array_equal(cols.arrays['XREST'], cols.arrays['XCORR']) and np.array_equal(cols.arrays['YREST'], cols.arrays['YCORR'])

# vignetting weights: F25SRF2 only, a separate column from FLATW (there is no combined WEIGHT)
ys2 = np.full(m, 29)                                  # row 28: on the bottom ramp
vig = replace(on, vigcorr=True, vigbot_start=128.0, vigbot_end=28.0, vigbot_loss=0.05)
sens = dataset(t, xs, ys2, weight=w); sens.aperture = 'F25SRF2'
cols = tagsel.selection_columns(sens, None, None, 'image', vig)
assert np.allclose(cols.arrays['VIGW'], 1 / 0.95, rtol=1e-6) and 'WEIGHT' not in cols.arrays
assert np.allclose(cols.arrays['FLATW'], w)
sens.aperture = '52X0.2'
cols = tagsel.selection_columns(sens, None, None, 'image', vig)
assert 'VIGW' not in cols.arrays and 'WEIGHT' not in cols.arrays and np.allclose(cols.arrays['FLATW'], w)
# no flat field: the switch is on but nothing can be written, and the reason is recorded
noflat = dataset(t, xs, ys); noflat.flat_problem = 'Flat field not found'
cols = tagsel.selection_columns(noflat, None, None, 'image', on)
assert 'FLATW' not in cols.arrays and ('flat field', 'Flat field not found') in cols.skipped

# ---- the commented examples: no combined weight, one flag per line, a light curve and a count-rate image --------------
class Cols(dict):
    """A stand-in for the astropy table the examples read: columns by name, a length."""
    def __len__(self):
        return len(self['TIME'])


def run_example(code, number, **extra):
    """Run the commented example `number` (1 or 2) on arrays, with the file reads and the plotting replaced."""
    part = code.split(f'# ---- Example {number}')[1].split('# ---- Example')[0]
    lines = [line[2:] for line in part.splitlines()[1:] if line.startswith('# ')]
    body = '\n'.join(line for line in lines if line and not line.startswith((' ', 'import', 'from', 'with', '#', 'Needs')))
    for old, new in extra.items():
        assert old in body, old
        body = body.replace(old, new)
    nn = 4000
    gen = np.random.default_rng(1)
    time = np.sort(gen.uniform(0, 99, nn))
    table = Cols(TIME=time, SEL_TIME=(time >= 10) & (time < 20), SEL_REGION=gen.random(nn) < 0.5, EDGE_OK=np.ones(nn, bool),
                 FLATW=np.full(nn, 2.0), VIGW=np.full(nn, 0.5), XCORR=gen.uniform(1, 2048, nn), YCORR=gen.uniform(1, 2048, nn))
    noop = lambda *args, **kwargs: None
    names = dict(np=np, OUTPUT='x', t=table, gti={'START': np.array([0.0, 50.0]), 'STOP': np.array([40.0, 99.0])},
                 fits=types.SimpleNamespace(getheader=lambda f, h: {'SELTMIN': 10.0, 'SELTMAX': 20.0}),
                 plt=types.SimpleNamespace(step=noop, xlabel=noop, ylabel=noop, show=noop, imshow=noop, colorbar=noop))
    exec(compile(body, f'example{number}', 'exec'), names)
    return names, table


dsv = dataset(np.linspace(0, 99, n), np.full(n, 600), np.full(n, 700), weight=np.full(n, 1.25)); dsv.aperture = 'F25SRF2'
s = tte.build_selection_script(dsv, (10.0, 20.0), Box(500, 700, 600, 800), 'image', replace(corr_all, vigcorr=True))
code = s.code
for flag in ('SEL_TIME', 'SEL_REGION', 'EDGE_OK'):
    assert code.count(f"# keep &= t['{flag}']") == 2, flag                  # one line per criterion, in both examples
assert code.count('# keep = np.ones(len(t), dtype=bool)') == 2
assert code.count("# w *= t['FLATW']") == 2 and code.count("# w *= t['VIGW']") == 2 and 'WEIGHT' not in code
assert code.count("hdr['SELTMIN']") == 2 and code.count('UnitsWarning') == 4
assert code.index('Example 1') < code.index('Example 2') < code.index('np.histogram2d')     # the light curve comes first
names, table = run_example(code, 1, **{'BIN = 10.0': 'BIN = 4.0'})
keep = table['SEL_TIME'] & table['SEL_REGION']
assert np.allclose(names['edges'], [10, 14, 18, 20]) and np.allclose(names['exposure'], [4, 4, 2])   # (inside the GTI)
assert names['counts'].sum() == keep.sum()                                  # FLATW x VIGW = 2 x 0.5 = 1 per event
assert np.allclose(names['rate'], names['counts'] / names['exposure'])
names, table = run_example(code, 2)
assert names['exposure'] == 10.0 and abs(names['image'].sum() * names['exposure'] - keep.sum()) < 1e-6   # counts per second
assert names['image'].shape == (1024, 1024)
# no time selection: the whole exposure, and no header read
s = tte.build_selection_script(ds, None, Box(500, 700, 600, 800), 'image', corr_all)
assert "keep &= t['SEL_TIME']" not in s.code and 'SELTMIN' not in s.code.split('Example 1')[1]
names, table = run_example(s.code, 2)
assert names['exposure'] == 40.0 + 49.0 and names['t0'] == 0.0
# with every correction off the examples use the recorded pixels and count each event once
s = tte.build_selection_script(ds, (10.0, 20.0), None, 'image', OFF)
assert "t['AXIS1']" in s.code and "w *=" not in s.code and 'w = np.ones(len(t))' in s.code

# ---- no selection: the corrections still apply, and the file has no selection columns --------------------------------
s = tte.build_selection_script(ds, None, None, 'image', corr_all)
assert not s.warnings and 'none (the corrections are still applied' in s.code and 'TIME_WINDOW = None' in s.code
calls.clear(); run_script(s.code)
assert calls[0][1]['time_sel'] is None and calls[0][1]['region'] is None and calls[0][1]['corr'].randcorr
assert "# keep &= t['EDGE_OK']" in s.code and "keep &= t['SEL_TIME']" not in s.code     # the examples keep only the flags that exist
assert 'SELTMIN' not in s.code and 'gti["STOP"].max()' in s.code                         # no time window to read from the header
bare = tagsel.selection_columns(exp, None, None, 'image', on)
assert list(bare.arrays)[:3] == ['EDGE_OK', 'XCORR', 'YCORR'] and bare.flags == ('EDGE_OK',)
assert not {'SEL_TIME', 'SEL_REGION'} & set(bare.arrays)
nothing = tagsel.selection_columns(exp, None, None, 'image', OFF)
assert 'EDGE_OK' not in nothing.arrays and not nothing.flags and 'XCORR' not in nothing.arrays and 'XREST' in nothing.arrays

# ---- the primary header: which steps ran, which selections were made ---------------------------------------------
def header(time_sel, region, frame, corr, dsx=exp):
    cols = tagsel.selection_columns(dsx, time_sel, region, frame, corr)
    return {k: v for k, v, c in tagsel.primary_keywords(dsx, time_sel, region, frame, corr, cols, '2.0', 'https://x')}

h = header((20.0, 60.0), box, 'image', on)
assert [h[k] for k in tagsel.STEP_KEYWORDS] == ['PERFORM', 'OMIT', 'PERFORM', 'OMIT', 'PERFORM'], h     # rand, geo, flat, vign, expt
assert h['TIMESEL'] is True and h['SPCSEL'] is True and (h['SELTMIN'], h['SELTMAX']) == (20.0, 60.0)
assert h['SELREG'] == 'box' and h['SELFRAME'] == 'image' and h['RANDSEED'] == RANDOM_SEED
assert (h['SELEDGEL'], h['SELEDGER'], h['SELEDGEB'], h['SELEDGET']) == tuple(EDGE_TRIM)
assert h['SELVER'] == '2.0' and h['TTEVER'] == tte.__version__ and h['SELURL'] == 'https://x' and 'VIGBSTRT' not in h
h = header((20.0, 60.0), None, 'image', OFF)
assert h['TIMESEL'] is True and h['SPCSEL'] is False and 'SELREG' not in h
assert [h[k] for k in tagsel.STEP_KEYWORDS] == ['OMIT'] * 5 and 'RANDSEED' not in h and 'SELEDGEL' not in h
h = header(None, None, 'image', replace(OFF, geocorr=True))
assert h['TIMESEL'] is False and h['SPCSEL'] is False and h['GEOCORR'] == 'PERFORM' and h['RANDCORR'] == 'OMIT'
# a step is PERFORM only if it took effect: no flat found, or a vignetting model for another aperture
assert header(None, None, 'image', on, noflat)['FLATCORR'] == 'OMIT'
other_ap = dataset(t, xs, ys, weight=w); other_ap.aperture = '52X2'
assert header(None, None, 'image', vig, other_ap)['VIGNCORR'] == 'OMIT'
sens.aperture = 'F25SRF2'
h = header(None, None, 'image', vig, sens)
assert h['VIGNCORR'] == 'PERFORM' and (h['VIGBSTRT'], h['VIGBEND'], h['VIGBLOSS']) == (128.0, 28.0, 0.05)
assert (h['VIGTSTRT'], h['VIGTEND'], h['VIGTLOSS']) == (1550.0, 2047.0, 0.23)
assert all(len(k) <= 8 for k in h)                          # FITS keywords

# ---- helpers for the file ------------------------------------------------------------------------------------------
# flags worked out in time order go back to the file's row order
order = np.array([2, 0, 1])
assert np.array_equal(tagsel.in_file_order(np.array([True, False, False]), order), np.array([False, False, True]))
assert tagsel.in_file_order(np.array([1, 2]), None).tolist() == [1, 2]
# The file keeps each column's coordinate keywords (TCTYPn, TCRPXn ...): astropy writes those from attributes of the
# column, so a rebuilt column loses them unless they are passed on.
made = []
real_fits = tagsel.fits
tagsel.fits = types.SimpleNamespace(Column=lambda **kw: made.append(kw) or kw)
try:
    column = types.SimpleNamespace(unit='pixel', coord_type='RA---TAN', coord_unit='deg', coord_ref_point=512.5,
                                   coord_ref_value=150.25, coord_inc=-1e-6, time_ref_pos=None)
    hdu = types.SimpleNamespace(columns={'AXIS1': column}, data={'AXIS1': np.arange(3, dtype=np.int16)})
    tagsel.stored_column(hdu, 'AXIS1')
    assert made[0]['coord_ref_point'] == 512.5 and made[0]['coord_type'] == 'RA---TAN' and made[0]['coord_inc'] == -1e-6, made
    assert made[0]['format'] == 'I' and made[0]['unit'] == 'pixel'
    made.clear()
    tagsel.column_definition('XCORR', np.zeros(3, np.float32)); tagsel.column_definition('EDGE_OK', np.ones(3, bool))
    tagsel.column_definition('FLATW', np.ones(3, np.float32))
    assert [(c['name'], c['format'], c.get('unit')) for c in made] == [
        ('XCORR', 'E', 'pixel'), ('EDGE_OK', 'L', None), ('FLATW', 'E', None)], made
finally:
    tagsel.fits = real_fits
print('selection script ok')
