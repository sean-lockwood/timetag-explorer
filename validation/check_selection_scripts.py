#!/usr/bin/env python
"""Run the selection scripts the TIME-TAG Explorer generates, and check what they write.

For one exposure this builds a script for each combination of switches, reference frame, region tool and selection
(time, space or both), runs each one in a fresh folder with the real astropy (the scripts call this checkout's
`timetag_explorer.tagsel`), and compares the columns each writes with the Explorer's own numbers, event by event: the
flags, the corrected positions (XCORR, YCORR and one pair per reference frame) and the weights. It also checks that the
TAG file's own columns and tables come through unchanged, and the header keywords and HISTORY.

The data are downloaded once and shared by every run:
  * the first script is run in the empty data folder, so it downloads the TAG file, the SPT file and the flat field
    itself: that tests the download code and is the only request to MAST and the flat-field site (use
    --skip-download-test to fetch the TAG and SPT files from AWS S3 instead; the flat field then comes from its site);
  * every other run gets links to those files, so nothing else is downloaded.

Usage:
    python validation/check_selection_scripts.py OFAP44G5Q                 # a quick pass: corrections all on, all off
    python validation/check_selection_scripts.py OFAP44G5Q --switches all  # every distinct script
    python validation/check_selection_scripts.py OFAP44G5Q --list          # only list the cases
    python validation/check_selection_scripts.py --tag-file /path/to/ofap44g5q_tag.fits --spt-file /path/to/ofap44g5q_spt.fits

Needs the package's own requirements (pip install -e .) from a checkout of the repository; s3fs is optional (without it
the files come from MAST). Run it from the repository's top folder.

Planet frames are off unless --planet-frames is given (they need JPL Horizons, and a moon observation). The region is
tested at each event's own position, so besides the event-by-event comparison this checks that the Explorer's own
binned highlight (4-pixel bins) agrees with it except near the region's edge.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import os
import re
import shutil
import subprocess
import sys
import time
import types
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent.parent      # the repository's top folder
HST_S3_BUCKET, HST_S3_PREFIX = 'stpubdata', 'hst/public'    # the public HST archive on AWS
PIXEL_FRAMES = ('image', 'rest', 'rest_north')
PLANET_FRAMES = ('planet', 'planet_north')
SHAPES = ('box', 'circle', 'polygon')
SELECTIONS = ('time', 'space', 'both', 'none')
# (randcorr, geocorr, flatcorr, timecorr): "quick" is everything on and everything off, "all" is every combination
SWITCH_SETS = {
    'quick': [(True, True, True, True), (False, False, False, False)],
    'all': list(itertools.product((True, False), repeat=4)),
}


# ---- getting the data ---------------------------------------------------------------------------------------------

def fetch_file(rootname, suffix, dest):
    """Copy one file of an exposure from the public HST bucket on AWS, or else through MAST."""
    root = rootname.lower()
    errors = []
    for source in ('AWS S3', 'MAST'):
        part = dest.with_name(dest.name + '.part')
        try:
            if source == 'AWS S3':
                import s3fs
                s3fs.S3FileSystem(anon=True).get(
                    f'{HST_S3_BUCKET}/{HST_S3_PREFIX}/{root[:4]}/{root}/{root}_{suffix}.fits', str(part))
            else:
                from astroquery.mast import Observations
                result = Observations.download_file(
                    f'mast:HST/product/{root}_{suffix}.fits', local_path=str(part), cache=False, verbose=False)
                status = result[0] if isinstance(result, (tuple, list)) else result
                if status == 'ERROR' or not part.is_file():
                    raise FileNotFoundError(f'MAST download failed: {result}')
            os.replace(part, dest)
            return
        except Exception as exc:   # network errors, a missing optional package, no such file ...
            errors.append(f'{source}: {exc}')
        finally:
            part.unlink(missing_ok=True)
    raise OSError(f'{dest.name} could not be downloaded ({"; ".join(errors)})')


def ensure_files(rootname, folder):
    """The exposure's TAG file (required) and SPT file (optional: only the rest frames need it) in `folder`."""
    root = rootname.lower()
    for suffix in ('tag', 'spt'):
        dest = folder / f'{root}_{suffix}.fits'
        if dest.is_file():
            continue
        try:
            fetch_file(root, suffix, dest)
        except OSError as exc:
            if suffix == 'tag':
                raise
            print(f'note: no SPT file ({exc}); the rest frames will be skipped')


# ---- the cases ---------------------------------------------------------------------------------------------

def make_cases(frames, shapes, selections, switches):
    """(frame, shape, selection, (rand, geo, flat, edge)) for every distinct script.

    A time-only (or no-selection) script does not depend on the frame or the shape, so it is listed once, as 'image'."""
    cases = []
    for sw in switches:
        if 'time' in selections:
            cases.append(('image', None, 'time', sw))
        if 'none' in selections:       # no selection: the corrections alone
            cases.append(('image', None, 'none', sw))
        for frame, shape, selection in itertools.product(frames, shapes, [s for s in selections if s not in ('time', 'none')]):
            cases.append((frame, shape, selection, sw))
    return cases


def case_name(case):
    frame, shape, selection, (rand, geo, flat, edge) = case
    flags = ''.join(c if on else '-' for c, on in zip('RGFE', (rand, geo, flat, edge)))
    return f'{frame}/{shape or "-"}/{selection}/{flags}'


class RegionMaker:
    """Regions and a time window that land on real events, so a mismatch cannot hide in an empty selection."""

    def __init__(self, te, ds):
        self.te, self.ds, self._uv = te, ds, {}
        lo, hi = np.percentile(ds.time, [20, 60])
        self.time_sel = (round(float(lo), te.TIME_DECIMALS), round(float(hi), te.TIME_DECIMALS))

    def uv(self, frame, corr):
        """A thinned sample of the events' coordinates in this frame, from the app's own grid."""
        key = (frame, corr.randcorr, corr.geocorr)
        if key not in self._uv:
            g = self.ds.grid(frame, corr)
            step = max(1, len(g.bin_index) // 2_000_000)
            idx = g.bin_index[::step].astype(np.int64)
            if g.keep is not None:
                idx = idx[g.keep[::step]]
            self._uv[key] = (g.xlo + (idx % g.nx + 0.5) * g.binsize, g.ylo + (idx // g.nx + 0.5) * g.binsize)
        return self._uv[key]

    def region(self, shape, frame, corr):
        te = self.te
        u, v = self.uv(frame, corr)
        u25, u50, u75 = np.percentile(u, [25, 50, 75])
        v25, v50, v75 = np.percentile(v, [25, 50, 75])
        r = lambda a: round(float(a), 1)
        if shape == 'box':
            return te.Box(r(u25), r(u75), r(v25), r(v75))
        if shape == 'circle':
            return te.Circle(r(u50), r(v50), r(max(u75 - u25, v75 - v25) / 2))
        if shape == 'polygon':
            radius = max(u75 - u25, v75 - v25) / 2
            points = []
            for k, scale in enumerate((1.0, 0.8, 1.1, 0.9, 1.0)):   # an irregular pentagon
                angle = np.pi / 2 + k * 2 * np.pi / 5
                points.append((r(u50 + scale * radius * np.cos(angle)), r(v50 + scale * radius * np.sin(angle))))
            return te.Polygon(tuple(points))
        raise ValueError(shape)


# ---- what the Explorer says the flags should be ---------------------------------------------------------------

def to_file_order(flag, order):
    """The Explorer works on time-ordered events; the script's columns follow the file's own row order."""
    if order is None:
        return flag
    out = np.empty_like(flag)
    out[order] = flag
    return out


def reference_columns(te, ds, order, case, time_sel, region, north_up=False):
    """The columns the script should write (in the file's row order), from the Explorer's own functions, in the order
    they are written: flags, then positions, then weights. The north-up frames' columns are left out unless
    `north_up` (the scripts' INCLUDE_NORTH_UP, off by default); there is no combined WEIGHT column."""
    frame, shape, selection, (rand, geo, flat, edge) = case
    corr = te.Corrections(randcorr=rand, geocorr=geo, flatcorr=flat, timecorr=edge)
    flags, positions, weights = {}, {}, {}
    if time_sel is not None:
        lo, hi = np.searchsorted(ds.time, time_sel)
        m = np.zeros(len(ds.time), dtype=bool)
        m[lo:hi] = True
        flags['SEL_TIME'] = m
    geo_lut = te.load_geo_lut() if geo else None
    xy = lambda f: te.frame_coordinates(ds.x, ds.y, ds.time, ds.frame_info, f, corr, seed=te.RANDOM_SEED,
                                        geo_lut=geo_lut, track=ds.planet_track() if f in PLANET_FRAMES else None)
    if region is not None:
        flags['SEL_REGION'] = np.asarray(region.contains(*xy(frame)), dtype=bool)   # at each event's own position
    valid = ds.valid_events(corr)
    if valid is not None:
        flags['EDGE_OK'] = np.array(valid, dtype=bool)
    for f, names in (('image', ('XCORR', 'YCORR')),) + tuple(te.POSITION_COLUMNS.items()):
        if f == 'image' and not (rand or geo):
            continue
        if f.endswith('_north') and not north_up:
            continue
        try:
            if f != 'image':
                ds.prepare_frame(f)
        except ValueError:
            continue          # a frame this exposure does not have: no columns
        u, v = xy(f)
        positions[names[0]], positions[names[1]] = u + np.float32(1.0), v + np.float32(1.0)
    flat_w = ds.weight if flat and ds.weight is not None else None
    gain = ds.vignette_gain(corr)
    if flat_w is not None:
        weights['FLATW'] = np.asarray(flat_w, dtype=np.float32)
    if gain is not None:
        weights['VIGW'] = gain
    return {name: to_file_order(values, order) for name, values in {**flags, **positions, **weights}.items()}


# ---- running and checking one script -----------------------------------------------------------------------------

def link_file(source, folder):
    link = folder / source.name
    if link != source and not link.exists():   # (the download test runs in the data folder itself)
        try:
            os.symlink(source, link)
        except OSError:   # e.g. Windows without the privilege
            shutil.copyfile(source, link)


def run_script(code, filename, tag_path, folder, timeout, shared=()):
    """Run the script in a fresh folder holding only links to the shared files (TAG, SPT, flat field). The script
    imports this checkout's package; the planet-frame ephemerides are kept in one shared folder."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / filename).write_text(code)
    for source in (tag_path, *shared):
        link_file(source, folder)
    env = dict(os.environ, PYTHONPATH=str(HERE) + os.pathsep + os.environ.get('PYTHONPATH', ''))
    env.setdefault('TAG_EXPLORER_HORIZONS', str(tag_path.parent / 'horizons'))
    started = time.time()
    proc = subprocess.run([sys.executable, filename], cwd=folder, capture_output=True, text=True, timeout=timeout, env=env)
    (folder / 'log.txt').write_text(f'exit {proc.returncode}\n--- stdout\n{proc.stdout}\n--- stderr\n{proc.stderr}\n')
    return proc, time.time() - started


def same(a, b):
    return a.shape == b.shape and bool(np.array_equal(a, b))


# Keywords that describe the table's layout, which the script rebuilds: column definitions (TTYPEn, TFORMn, TSCALn ...)
# and the sizes. (The column coordinate keywords TCTYPn, TCRPXn ... are NOT in this list: they must survive.)
STRUCTURE_KEYWORD = re.compile(r'T(TYPE|FORM|UNIT|NULL|SCAL|ZERO|DISP|BCOL|DIM)\d+|XTENSION|BITPIX|NAXIS\d*|PCOUNT|GCOUNT'
                               r'|TFIELDS|THEAP')


def same_value(a, b):
    """Header values compared as written (floats to 12 digits, since astropy re-formats them)."""
    if isinstance(a, float) or isinstance(b, float):
        try:
            return bool(np.isclose(float(a), float(b), rtol=1e-12, atol=0.0))
        except (TypeError, ValueError):
            return False
    return a == b


def check_output(te, fits, case, region, expected, tag_path, folder, rootname):
    """Problems found in the script's output file (an empty list means it is right)."""
    frame, shape, selection, (rand, geo, flat, edge) = case
    problems = []
    out_path = folder / f'{rootname}_tagsel.fits'
    if not out_path.is_file():
        return [f'{out_path.name} was not written']
    stray = sorted(p.name for p in folder.iterdir()
                   if p.name not in (out_path.name, 'log.txt') and not p.name.endswith('.py') and not p.is_symlink()
                   and p.name != 'horizons')
    if stray:
        problems.append(f'stray files left behind: {stray}')
    with fits.open(tag_path) as old, fits.open(out_path) as new:
        try:
            new.verify('exception')
        except Exception as exc:   # astropy raises its own VerifyError
            problems.append(f'FITS verify: {exc}')
        if len(new) != len(old) or [h.name for h in new] != [h.name for h in old]:
            return problems + [f'extensions differ: {[h.name for h in old]} -> {[h.name for h in new]}']
        old_names = list(old[1].columns.names)
        want = old_names + list(expected)
        got = list(new[1].columns.names)
        if got != want:
            problems.append(f'columns {got}, expected {want}')
        for name in old_names:   # includes AXIS1 / AXIS2: still the recorded integers
            if name in got and not same(np.asarray(new[1].data[name]), np.asarray(old[1].data[name])):
                problems.append(f'column {name} changed')
        for name, ref in expected.items():
            if name in got:
                values = np.asarray(new[1].data[name])
                if ref.dtype == bool:
                    values = values.astype(bool)
                if not same(values, ref):
                    differing = int(np.count_nonzero(values != ref)) if values.shape == ref.shape else -1
                    problems.append(f'{name}: {differing} of {len(ref):,} events differ from the Explorer')
        for i in range(2, len(old)):   # other tables (good-time intervals ...) are carried over untouched
            if getattr(old[i].data, 'names', None):
                for name in old[i].columns.names:
                    if not same(np.asarray(new[i].data[name]), np.asarray(old[i].data[name])):
                        problems.append(f'extension {old[i].name} column {name} changed')
        header = new[0].header
        host = te.tool_url().split('://')[-1]
        history = '\n'.join(str(c) for c in new[1].header['HISTORY']) + '\n' + '\n'.join(str(c) for c in header['HISTORY'])
        # which steps ran: PERFORM only when the step took effect (its columns were written)
        done = {'RANDCORR': rand, 'GEOCORR': geo, 'FLATCORR': 'FLATW' in expected, 'VIGNCORR': 'VIGW' in expected,
                'EXPTCORR': 'EDGE_OK' in expected}
        checks = [(header.get(key) == ('PERFORM' if on else 'OMIT'), f'{key} {header.get(key)!r}') for key, on in done.items()]
        checks += [(header.get('TIMESEL') is (selection in ('time', 'both')), f'TIMESEL {header.get("TIMESEL")!r}'),
                   (header.get('SPCSEL') is (region is not None), f'SPCSEL {header.get("SPCSEL")!r}'),
                   (('SEL_TIME' in expected) == (selection in ('time', 'both')), 'SEL_TIME column presence'),
                   (('SEL_REGION' in expected) == (region is not None), 'SEL_REGION column presence'),
                   (header.get('SELVER') == te.SELECTION_SCRIPT_VERSION, f'SELVER {header.get("SELVER")!r}'),
                   (header.get('TTEVER') == te.__version__, f'TTEVER {header.get("TTEVER")!r}'),
                   (host in str(header.get('SELURL', '')), f'SELURL {header.get("SELURL")!r}'),
                   (host in history, 'tool URL missing from HISTORY'),
                   ('v' + te.SELECTION_SCRIPT_VERSION in history, 'script version missing from HISTORY')]
        if done['EXPTCORR']:
            checks += [(header.get(key) == pixels, f'{key} {header.get(key)!r}, expected {pixels}')
                       for key, pixels in zip(('SELEDGEL', 'SELEDGER', 'SELEDGEB', 'SELEDGET'), te.EDGE_TRIM)]
        if region is not None:
            checks.append((header.get('SELFRAME') == frame, f'SELFRAME {header.get("SELFRAME")!r}'))
        for ok, text in checks:
            if not ok:
                problems.append('header: ' + text)
        # every other keyword of the event table survives with the same value; only the table's own structure
        # (column definitions, which astropy writes afresh, and the sizes) may differ
        lost, changed = [], []
        for key in dict.fromkeys(old[1].header):
            if STRUCTURE_KEYWORD.fullmatch(key) or key in ('COMMENT', 'HISTORY', ''):
                continue
            if key not in new[1].header:
                lost.append(key)
            elif not same_value(new[1].header[key], old[1].header[key]):
                changed.append(key)
        # ... and the primary header keeps its keywords too, except the five step flags the script sets itself
        steps = ('RANDCORR', 'GEOCORR', 'FLATCORR', 'VIGNCORR', 'EXPTCORR')
        for key in dict.fromkeys(old[0].header):
            if key in steps or key in ('COMMENT', 'HISTORY', ''):
                continue
            if key not in new[0].header:
                lost.append('primary ' + key)
            elif not same_value(new[0].header[key], old[0].header[key]):
                changed.append('primary ' + key)
        if lost:
            problems.append(f'header keywords lost: {", ".join(lost)}')
        if changed:
            problems.append(f'header keywords changed: {", ".join(changed)}')
    return problems


# ---- the driver -----------------------------------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('rootname', nargs='?', help='HST exposure to test, e.g. OFAP44G5Q (not needed with --tag-file)')
    p.add_argument('--tag-file', help='use this TAG file instead of downloading (also stops the download test)')
    p.add_argument('--spt-file', help='this exposure\'s SPT file, needed for the rest frames')
    p.add_argument('--work', default='selection_test', help='working folder (default: ./selection_test)')
    p.add_argument('--switches', choices=sorted(SWITCH_SETS), default='quick',
                   help='quick: corrections all on or all off. all: every combination of jitter, geometric, flat and edge trim')
    p.add_argument('--frames', default=','.join(PIXEL_FRAMES), help='comma-separated (default: image,rest,rest_north)')
    p.add_argument('--planet-frames', action='store_true', help='also try the planet frames (needs Horizons, and a moon observation)')
    p.add_argument('--shapes', default=','.join(SHAPES))
    p.add_argument('--selections', default=','.join(SELECTIONS))
    p.add_argument('--limit', type=int, help='run only this many cases, picked at random')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--jobs', type=int, default=2, help='scripts running at once (each reads the whole file)')
    p.add_argument('--timeout', type=int, default=3600, help='seconds allowed for one script')
    p.add_argument('--skip-download-test', action='store_true', help='fetch the TAG and SPT files from AWS S3 (or MAST) here instead of testing the script\'s own downloads')
    p.add_argument('--keep', action='store_true', help='keep every run folder, not only the failures')
    p.add_argument('--list', action='store_true', help='list the cases and stop')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not args.rootname and not args.tag_file:
        sys.exit('give an exposure name (e.g. OFAP44G5Q) or --tag-file')
    frames = [f for f in args.frames.split(',') if f] + (list(PLANET_FRAMES) if args.planet_frames else [])
    shapes = [s for s in args.shapes.split(',') if s]
    selections = [s for s in args.selections.split(',') if s]
    cases = make_cases(frames, shapes, selections, SWITCH_SETS[args.switches])
    if args.limit and args.limit < len(cases):
        rng = np.random.default_rng(args.seed)
        cases = [cases[i] for i in sorted(rng.choice(len(cases), args.limit, replace=False))]
    print(f'{len(cases)} cases')
    if args.list:
        for c in cases:
            print('  ', case_name(c))
        return 0

    work = Path(args.work).resolve()
    data = work / 'data'
    data.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(HERE))
    import timetag_explorer as tte       # the package's own code: the generator and the reference masks
    from timetag_explorer import constants, dataset as tte_dataset
    from astropy.io import fits

    from timetag_explorer import tagsel, transforms
    te = types.SimpleNamespace(
        Box=tte.Box, Circle=tte.Circle, Polygon=tte.Polygon, Corrections=tte.Corrections, FrameInfo=tte.FrameInfo,
        TagDataset=tte.TagDataset, build_selection_script=tte.build_selection_script, __version__=tte.__version__,
        SELECTION_SCRIPT_VERSION=tte.SELECTION_SCRIPT_VERSION, TIME_DECIMALS=constants.TIME_DECIMALS,
        EDGE_TRIM=constants.EDGE_TRIM, tool_url=tte.selection_script.tool_url, RANDOM_SEED=constants.RANDOM_SEED,
        POSITION_COLUMNS=tagsel.POSITION_COLUMNS, frame_coordinates=transforms.frame_coordinates,
        load_geo_lut=transforms.load_geo_lut,
        _is_sorted=tte_dataset._is_sorted, _load_dataset_uncached=tte_dataset._load_dataset_uncached)

    if args.tag_file:
        tag_path = Path(args.tag_file).resolve()
        rootname = tag_path.name.split('_', 1)[0].lower()
        if args.spt_file:
            spt = Path(args.spt_file).resolve()
            link = tag_path.parent / spt.name
            if not link.exists() and spt.parent != tag_path.parent:
                print(f'note: copy {spt.name} next to the TAG file so the rest frames can use it')
    else:
        rootname = args.rootname.lower()
        tag_path = data / f'{rootname}_tag.fits'

    # 1. The download test: a script run in the empty data folder fetches the TAG file, the SPT file and the flat
    #    field itself.
    results = []
    flat_path = data / 'combined_flat.fits'
    if not args.tag_file and not tag_path.exists() and not args.skip_download_test:
        print(f'download test: running a script in the empty folder {data} (requests to MAST and the flat-field site) ...')
        stub = te.TagDataset(name=tag_path.name, time=np.array([0.0, 1.0], np.float32), x=np.array([1, 2], np.int16),
                             y=np.array([1, 2], np.int16), gti=np.array([[0.0, 1.0]]), rootname=rootname,
                             frame_info=te.FrameInfo(orientat=0.0, slew_angle=0.0, slew_rate=0.0), weight=np.ones(2, np.float32))
        s = te.build_selection_script(stub, (0.0, 1.0), None, 'image', te.Corrections(randcorr=False, geocorr=False))
        proc, secs = run_script(s.code, s.filename, tag_path, data, args.timeout)   # (the link is skipped: no file yet)
        ok = (proc.returncode == 0 and tag_path.is_file() and not tag_path.is_symlink()
              and (data / f'{rootname}_spt.fits').is_file() and flat_path.is_file())
        results.append(('download test', 'PASS' if ok else 'FAIL', secs, '' if ok else (proc.stderr or proc.stdout)[-600:]))
        print(f'  {results[-1][1]} ({secs:.0f} s)')
        if not ok:
            print(proc.stdout + proc.stderr)
            return 1
        out = data / f'{rootname}_tagsel.fits'
        out.unlink(missing_ok=True)
    # 2. The data every run shares (the SPT file comes from S3, or MAST; the flat field from its site).
    if not tag_path.exists():
        print('fetching the data ...')
        ensure_files(rootname, data)
    elif not args.tag_file:
        ensure_files(rootname, data)   # only fetches what is missing: the SPT file
    if not flat_path.exists():
        print('fetching the flat field ...')
        tagsel.download_flat(flat_path)
    os.environ['TAG_EXPLORER_FLAT'] = str(flat_path)
    os.environ['TAG_EXPLORER_HORIZONS'] = str(data / 'horizons')
    tte.config.FLAT_FILE, tte.config.HORIZONS_DIR = flat_path, data / 'horizons'
    spt_path = tag_path.with_name(f'{rootname}_spt.fits')
    shared = [p for p in (spt_path, flat_path) if p.is_file()]
    print(f'loading {tag_path.name} ...')
    ds = te._load_dataset_uncached(str(tag_path))
    if ds.frame_info.problem:
        print(f'note: rest frames unavailable ({ds.frame_info.problem})')
    if ds.weight is None:
        print(f'note: no flat-field weights ({ds.flat_problem})')
    raw_time = np.asarray(fits.getdata(tag_path, ext=1)['TIME'], dtype=np.float32)
    order = None if te._is_sorted(raw_time) else np.argsort(raw_time, kind='stable')
    del raw_time
    print(f'{len(ds.time):,} events, '
          f'{"in time order" if order is None else "NOT in time order: columns are mapped back to file order"}')
    maker = RegionMaker(te, ds)
    print(f'time window {maker.time_sel[0]} to {maker.time_sel[1]} s')

    # 3. Build every script, then run them.
    available, jobs, skipped = {}, [], []
    for case in cases:
        frame, shape, selection, (rand, geo, flat, edge) = case
        if frame not in available:
            try:
                ds.prepare_frame(frame)
                available[frame] = None
            except ValueError as exc:
                available[frame] = str(exc)
        if available[frame]:
            skipped.append((case, available[frame]))
            continue
        corr = te.Corrections(randcorr=rand, geocorr=geo, flatcorr=flat, timecorr=edge)
        region = None if selection in ('time', 'none') else maker.region(shape, frame, corr)
        time_sel = maker.time_sel if selection in ('time', 'both') else None
        try:
            s = te.build_selection_script(ds, time_sel, region, frame, corr)
        except Exception as exc:
            results.append((case_name(case), 'FAIL', 0.0, f'the Explorer could not build the script: {exc!r}'))
            continue
        jobs.append((case, region, time_sel, s))
    for case, why in skipped:
        results.append((case_name(case), 'SKIP', 0.0, why))
    print(f'running {len(jobs)} scripts, {args.jobs} at a time ...')

    runs = work / 'runs'
    shutil.rmtree(runs, ignore_errors=True)

    def run_one(job):
        case, region, time_sel, s = job
        folder = runs / case_name(case).replace('/', '_')
        return job, folder, *run_script(s.code, s.filename, tag_path, folder, args.timeout, shared)

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for future in as_completed([pool.submit(run_one, j) for j in jobs]):
            try:
                (case, region, time_sel, s), folder, proc, secs = future.result()
            except subprocess.TimeoutExpired as exc:
                results.append(('(timeout)', 'FAIL', float(args.timeout), str(exc)))
                continue
            name = case_name(case)
            if proc.returncode != 0:
                problems = ['script exited with an error: ' + (proc.stderr.strip().splitlines() or ['?'])[-1]]
            else:
                expected = reference_columns(te, ds, order, case, time_sel, region)
                try:
                    problems = check_output(te, fits, case, region, expected, tag_path, folder, rootname)
                except Exception as exc:
                    problems = [f'could not check the output: {exc!r}']
                if not problems and 'SEL_REGION' in expected:
                    n = int(expected['SEL_REGION'].sum())
                    if n == 0 or n == len(ds.time):   # a region with nothing, or everything, in it proves little
                        problems.append(f'weak test: the region holds {n:,} events')
                    else:   # the Explorer's binned highlight differs from the exact test only near the region's edge
                        corr = te.Corrections(randcorr=case[3][0], geocorr=case[3][1], flatcorr=case[3][2],
                                              timecorr=case[3][3])
                        binned = to_file_order(ds.region_mask(region, case[0], corr), order)
                        apart = int(np.count_nonzero(binned != expected['SEL_REGION']))
                        if apart > 0.10 * max(n, 1):
                            problems.append(f'the exact region test and the Explorer\'s binned one differ for {apart:,} '
                                            f'of {n:,} selected events: more than an edge effect')
            status = 'PASS' if not problems else 'FAIL'
            results.append((name, status, secs, '; '.join(problems)))
            print(f'  {status}  {name:<42} {secs:6.1f} s  {"; ".join(problems)[:300]}')
            if status == 'PASS' and not args.keep:
                shutil.rmtree(folder, ignore_errors=True)   # each output is as big as the TAG file
            else:
                (folder / 'case.txt').write_text(
                    f'{name}\nregion {region}\ntime {time_sel}\n' + '\n'.join(problems) + '\n')

    # 4. Summary
    counts = {k: sum(1 for r in results if r[1] == k) for k in ('PASS', 'FAIL', 'SKIP')}
    for r in results:
        if r[1] == 'SKIP':
            print(f'  SKIP  {r[0]}: {r[3]}')
    work.mkdir(parents=True, exist_ok=True)
    with open(work / 'results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case', 'status', 'seconds', 'problems'])
        w.writerows(sorted(results))
    print(f'\n{counts["PASS"]} passed, {counts["FAIL"]} failed, {counts["SKIP"]} skipped  (results: {work / "results.csv"})')
    if counts['FAIL']:
        print(f'failing runs are kept in {runs} (script, log.txt, and case.txt with the region and window)')
    return 1 if counts['FAIL'] else 0


if __name__ == '__main__':
    sys.exit(main())
