"""The "Selection code" feature: writes a standalone Python script that reproduces a user's time and region
selections on the original TIME-TAG file and saves the result as flag columns in a copy of it.

Raise SELECTION_SCRIPT_VERSION (constants.py) whenever `build_selection_script` changes what the script computes or
writes."""
from __future__ import annotations

import textwrap
from dataclasses import dataclass
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

from . import config, gdcorr
from .constants import (DEFAULT_TOOL_URL, EDGE_TRIM, EVENT_CHUNK, RANDOM_SEED, SELECTION_SCRIPT_VERSION,
                        SPATIAL_BIN, TIME_DECIMALS, TimeSel)
from .dataset import TagDataset
from .frames import ARCSEC_PER_LOWRES_PIXEL, DETECTOR_CENTRE, Corrections, frame_label, is_map
from .regions import Box, Circle, Polygon, Region

# Selection code: a script that reproduces the user's selections outside the app
# --------------------------------------------------------------------------

PIXEL_FRAMES = ('image', 'rest', 'rest_north')  # frames that need nothing beyond numbers from the header
EDGE_TRIM_TEXT = f'left {EDGE_TRIM.left}, right {EDGE_TRIM.right}, bottom {EDGE_TRIM.bottom}, top {EDGE_TRIM.top}'


@dataclass(frozen=True)
class SelectionScript:
    """A runnable Python script for the current selections, and what its reader needs to know."""
    code: str
    filename: str
    packages: tuple[tuple[str, str], ...]  # (package, what it is for)
    pip_command: str
    notes: tuple[str, ...]      # what the script reproduces and what it leaves out
    warnings: tuple[str, ...]   # things that will not work as the reader might expect


def _region_source(region) -> tuple[list[str], list[str]]:
    """(constants, function) lines that define `region_contains(x, y)` for a region shape,
    with the same tests as Box / Circle / Polygon above."""
    f = lambda v: repr(float(v))  # noqa: E731  (exact literals)
    if isinstance(region, Box):
        return ([f'X0, X1 = {f(region.x0)}, {f(region.x1)}', f'Y0, Y1 = {f(region.y0)}, {f(region.y1)}'],
                ['def region_contains(x, y):',
                 '    """Rectangle: x0 <= x < x1 and y0 <= y < y1."""',
                 '    return (x >= X0) & (x < X1) & (y >= Y0) & (y < Y1)'])
    if isinstance(region, Circle):
        return ([f'CX, CY, R = {f(region.cx)}, {f(region.cy)}, {f(region.r)}'],
                ['def region_contains(x, y):',
                 '    """Circle: distance from the centre <= r."""',
                 '    return (x - CX) ** 2 + (y - CY) ** 2 <= R ** 2'])
    points = [f'({f(px)}, {f(py)})' for px, py in region.points]
    rows = [', '.join(points[i:i + 3]) + ',' for i in range(0, len(points), 3)]
    return (['POLYGON = ['] + ['    ' + row for row in rows] + [']'],
            ['def region_contains(x, y):',
             '    """Free-form region: point-in-polygon test."""',
             '    return MplPath(np.asarray(POLYGON, dtype=float)).contains_points(np.column_stack((x, y)))'])


@lru_cache(maxsize=1)
def gdcorr_source() -> str | None:
    """The text of gdcorr.py, to be included in the scripts that need it (None if it cannot be read)."""
    try:
        return Path(gdcorr.__file__).read_text()
    except (OSError, TypeError):
        return None


def tool_url() -> str:
    """Web address of this tool, for the scripts and FITS headers it writes."""
    return config.PUBLIC_URL or DEFAULT_TOOL_URL


def _region_exact(region) -> str:
    """The region's numbers, unrounded, for the FITS HISTORY."""
    f = lambda v: repr(float(v))  # noqa: E731
    if isinstance(region, Box):
        return f'box x0={f(region.x0)} x1={f(region.x1)} y0={f(region.y0)} y1={f(region.y1)}'
    if isinstance(region, Circle):
        return f'circle cx={f(region.cx)} cy={f(region.cy)} r={f(region.r)}'
    if len(region.points) > 40:
        return f'free-form polygon, {len(region.points)} vertices (listed in the script)'
    return 'free-form polygon, vertices (x,y): ' + ' '.join(f'({f(px)},{f(py)})' for px, py in region.points)


def build_selection_script(ds: TagDataset, time_sel: TimeSel | None, region: Region | None,
                           frame: str, corr: Corrections) -> SelectionScript:
    """A Python script (numpy, astropy, astroquery) that downloads the exposure's TAG file, flags the events
    the Explorer has selected, and writes the file back out with the flags as extra columns.

    The script repeats the Explorer's own arithmetic for the pixel frames: the same seeded random
    offsets, geometric lookup, frame shift and 4-pixel bin test, so it selects the same events.
    Planet frames and maps need ephemerides, so for those it carries the region definition only.
    Change what the script computes or writes only together with SELECTION_SCRIPT_VERSION.
    """
    info = ds.frame_info
    root = (ds.rootname or ds.name.split('_')[0]).lower()
    label = frame_label(info, frame)
    pixel_frame = frame in PIXEL_FRAMES
    placeable = region is not None and pixel_frame   # the script can compute each event's position
    geo = bool(corr.geocorr) and placeable
    jitter = bool(corr.randcorr) and placeable
    rest = placeable and frame != 'image'
    edge = bool(corr.timecorr) and any(EDGE_TRIM)
    grid = ds.grid(frame, corr) if region is not None else None
    url = tool_url()
    version = SELECTION_SCRIPT_VERSION
    geo_source = gdcorr_source() if geo else None
    embed_geo = geo_source is not None  # (else the script imports gdcorr)
    t_text = (f'{time_sel[0]:.{TIME_DECIMALS}f} to {time_sel[1]:.{TIME_DECIMALS}f} s' if time_sel is not None else '')

    # ---- what the reader needs -----------------------------------------------------------------
    packages = [('numpy', 'arrays'),
                ('astropy', 'reading and writing the FITS files'),
                ('astroquery', 'downloading the TAG file from MAST (only used if it is not in the folder already)')]
    if isinstance(region, Polygon):
        packages.append(('matplotlib', 'the point-in-polygon test for a free-form region'))
    if geo and not embed_geo:
        packages.append(('gdcorr.py', 'the STIS FUV geometric distortion polynomial. It comes with the '
                                      'TIME-TAG Explorer source and is not on PyPI: save it next to the script'))
    pip = 'pip install numpy astropy astroquery' + (' matplotlib' if isinstance(region, Polygon) else '')

    ran = []  # corrections the script runs, for the headers
    if jitter:
        ran.append(f'random sub-pixel offset (+/-0.5 hi-res pixel, seed {RANDOM_SEED})')
    if geo:
        ran.append('geometric distortion (gdcorr)')
    if rest:
        ran.append('target rest-frame shift' + (' and north-up rotation' if frame == 'rest_north' else ''))
    if edge:
        ran.append(f'detector edge trim ({EDGE_TRIM_TEXT} pixels)')

    notes, warnings = [], []
    notes.append(f'Run it in an empty folder: it downloads `{root}_tag.fits` from MAST if the file is not there '
                 f'(needs astroquery and a network connection), then writes `{root}_tagsel.fits`, the TAG file’s own '
                 'tables with the flag columns added. That is not an STScI calibrated product, and the '
                 'headers say so. It will not replace an existing file unless `OVERWRITE = True`.')
    if time_sel is not None:
        notes.append('The time cut compares the Explorer’s float32 copy of TIME, so it picks the same events as '
                     'the plots. Cutting the full-precision TIME instead can differ by about one event at an edge.')
    if placeable:
        notes.append(f'The region is tested the way the Explorer tests it: on {SPATIAL_BIN}-pixel bins, an event is '
                     'inside when its bin’s centre is.')
        if jitter or rest:
            notes.append('Which bin a jittered or shifted event lands in depends on the random sub-pixel offsets. '
                         'The script draws them with the same seeded NumPy generator, in the same order, so the '
                         'flags agree with the Explorer’s. That relies on NumPy giving the same stream, which '
                         'it does in practice but does not formally promise between versions.')
    notes.append('TIME is stored in the TAG file as a scaled integer; the output holds the same values as plain '
                 '64-bit floating point, so it is a little larger. Other columns are unchanged.')
    if geo and embed_geo:
        notes.append('The source of `gdcorr.py` (the geometric distortion polynomial) is included in the script.')
    if region is not None and not pixel_frame:
        warnings.append(
            f'The region was drawn in the {label}, which needs JPL Horizons ephemerides'
            + (' and the planet’s rotation model' if is_map(frame) else '')
            + ' to place each event. The script cannot compute that, so it contains the time selection and the '
              'region’s definition, but no region flag: pass your own coordinates in that frame to `in_region(u, v)`.')
    notes.append('Flat-field weights and the exposure-time map (counts to rates) are not applied: the script '
                 'selects events, it does not make images.')
    notes.append(f'Script version {version}. It is written into the output file’s header, so keep the number '
                 'if you report results.')

    # ---- the script ----------------------------------------------------------------------------
    L: list[str] = ['#!/usr/bin/env python']
    add = L.append
    corr_text = ', '.join(f'{name} {"on" if on else "off"}' for name, on in (
        ('random sub-pixel offset', corr.randcorr), ('geometric distortion', corr.geocorr),
        ('exposure-time (edge trim)', corr.timecorr)))
    add(f'"""Reproduce selections made in TIME-TAG Explorer on HST/STIS exposure {root}.')
    add('')
    add(f'TIME-TAG Explorer: {url}')
    add(f'Script version {version}, written {datetime.now(timezone.utc):%Y-%m-%d} (UTC).')
    add('')
    add('Selections in the Explorer')
    if time_sel is not None:
        add(f'  time         {t_text}')
    if region is not None:
        add(f'  region       {region.describe()}  ({label})')
    add(f'  corrections  {corr_text}')
    add('')
    add('What it does')
    add(f'  1. Downloads {root}_tag.fits from MAST into the current folder, unless it is already there.')
    add('  2. Flags the selected events (True = selected).')
    add(f'  3. Writes {root}_tagsel.fits: the TAG file with the flags added as columns, and header')
    add('     keywords and HISTORY that record the script version and the corrections that ran.')
    add('')
    add('Python packages needed')
    for name, why in packages:
        L.extend(textwrap.wrap(f'{name}: {why}', width=92, initial_indent='  ', subsequent_indent='      '))
    add(f'  Install with:  {pip}')
    add('')
    add(f'Run it:  python select_{root}.py')
    add('"""')
    add('import os')
    if rest:
        add('import math')
    add('import re')
    add('import textwrap')
    add('from datetime import datetime, timezone')
    add('')
    add('import numpy as np')
    add('from astropy.io import fits')
    if isinstance(region, Polygon):
        add('from matplotlib.path import Path as MplPath')
    add('')
    add(f"SCRIPT_VERSION = {version!r}   # version of the TIME-TAG Explorer selection script")
    add(f'TOOL_URL = {url!r}')
    add(f'ROOTNAME = {root!r}')
    add("TAG_FILE = ROOTNAME + '_tag.fits'      # the exposure's TIME-TAG file: downloaded from MAST if it is missing")
    add("OUTPUT = ROOTNAME + '_tagsel.fits'      # the TAG file with the flag columns added")
    add('OVERWRITE = False                      # True: replace OUTPUT if it exists already')
    flags = []
    if time_sel is not None:
        flags.append('SEL_TIME')
        add('')
        add('# ---- time selection ------------------------------------------------------------')
        add(f'TIME_WINDOW = ({float(time_sel[0])!r}, {float(time_sel[1])!r})   '
            '# seconds since the exposure start: start <= TIME < end')
    if region is not None:
        constants, function = _region_source(region)
        unit = ('map-plane degrees, as on the Explorer axes' if is_map(frame)
                else '0-indexed hi-res detector pixels (AXIS1 and AXIS2 count from 1)')
        add('')
        add('# ---- region --------------------------------------------------------------------')
        add(f'# Drawn in the {label}; coordinates are {unit}.')
        add(f'# The Explorer tests a region on a grid of {SPATIAL_BIN}-pixel bins: an event is inside when '
            'the centre of its bin is.' if not is_map(frame) else
            '# The Explorer tests a region on a grid of map bins: an event is inside when the centre of its '
            'bin is (and it lies on the planet).')
        L.extend(constants)
        add(f'BIN = {grid.binsize!r}   # bin size of the grid')
        add(f'GRID_X0, GRID_Y0 = {float(grid.xlo)!r}, {float(grid.ylo)!r}   # lower-left edge of the grid, in this frame')
        add(f'GRID_NX, GRID_NY = {int(grid.nx)}, {int(grid.ny)}   # grid size in bins')
        add('')
        add('')
        L.extend(function)
        add('')
        add('')
        add('def in_region(u, v):')
        add(f'    """True for events whose bin centre is inside the region; u, v are positions in the {label}."""')
        add('    centres_x = GRID_X0 + (np.arange(GRID_NX) + 0.5) * BIN')
        add('    centres_y = GRID_Y0 + (np.arange(GRID_NY) + 0.5) * BIN')
        add('    gx, gy = np.meshgrid(centres_x, centres_y)')
        add('    in_bin = region_contains(gx.ravel(), gy.ravel())   # row-major: iy * GRID_NX + ix')
        add('    ix = np.floor((u - GRID_X0) / BIN).astype(np.int32)')
        add('    iy = np.floor((v - GRID_Y0) / BIN).astype(np.int32)')
        add('    ok = (ix >= 0) & (ix < GRID_NX) & (iy >= 0) & (iy < GRID_NY)')
        add('    inside = np.zeros(len(u), dtype=bool)')
        add('    inside[ok] = in_bin[iy[ok] * GRID_NX + ix[ok]]')
        add('    return inside')
    if placeable:
        flags.append('SEL_REGION')
        add('')
        add('')
        add('# ---- where each event is in that frame -------------------------------------------')
        if jitter:
            add(f'SEED, CHUNK = {RANDOM_SEED}, {EVENT_CHUNK}   # the Random sub-pixel offsets of the Explorer: '
                '+/-0.5 pixel, drawn CHUNK events at a time')
        if rest:
            angle = 0.0 if info.slew_angle is None else float(info.slew_angle)
            orient = 0.0 if info.orientat is None else float(info.orientat)
            add(f'ORIENTAT = {orient!r}   # deg, ORIENTAT in the header of the event table')
            add(f'SLEW_ANGLE = {angle!r}   # deg, ANG of the last POS_ANGLE motion-tracking level (MT_LV keywords)')
            add(f'SLEW_RATE = {float(info.slew_rate)!r}   # arcsec/day, R_RAD of that level')
            add(f'ARCSEC_PER_LOWRES_PIXEL, DETECTOR_CENTRE = {ARCSEC_PER_LOWRES_PIXEL!r}, {DETECTOR_CENTRE!r}')
        if geo:
            add('')
            if embed_geo:
                add('# ---- geometric distortion: gdcorr.py, the STIS FUV polynomial, included as it is ------')
                L.extend(line.rstrip() for line in geo_source.strip('\n').split('\n'))
                add('')
                add('')
            else:
                add('try:')
                add('    import gdcorr   # the STIS FUV geometric distortion polynomial (comes with TIME-TAG Explorer)')
                add('except ImportError as exc:')
                add("    raise SystemExit('gdcorr.py not found: save it next to this script') from exc")
            add('_axis = np.arange(2048, dtype=np.float64)')
            add('_gx, _gy = np.meshgrid(_axis, _axis)   # [row = y, column = x]')
            add('_xc, _yc = ' + ('gdcorr(_gx, _gy)' if embed_geo else 'gdcorr.gdcorr(_gx, _gy)'))
            add('GEO_DX, GEO_DY = (_xc - _gx).astype(np.float32), (_yc - _gy).astype(np.float32)   # offset at each pixel')
        add('')
        add('')
        add('def frame_coordinates(x, y, time):')
        add(f'    """(u, v) of each event in the {label}: 0-indexed hi-res pixels, float32 as in the Explorer."""')
        add('    n = len(x)')
        add('    u = np.empty(n, dtype=np.float32)')
        add('    v = np.empty(n, dtype=np.float32)')
        add('    for lo in range(0, n, ' + ('CHUNK' if jitter else str(EVENT_CHUNK)) + '):')
        add('        hi = min(lo + ' + ('CHUNK' if jitter else str(EVENT_CHUNK)) + ', n)')
        add('        uu = x[lo:hi].astype(np.float32) - 1.0   # AXIS1 and AXIS2 count from 1')
        add('        vv = y[lo:hi].astype(np.float32) - 1.0')
        if geo:
            add('        px = np.clip(x[lo:hi].astype(np.int32) - 1, 0, 2047)')
            add('        py = np.clip(y[lo:hi].astype(np.int32) - 1, 0, 2047)')
            add('        uu += GEO_DX[py, px]   # geometric distortion')
            add('        vv += GEO_DY[py, px]')
        if jitter:
            add('        rng = np.random.default_rng(SEED if lo == 0 else [SEED, lo])')
            add('        uu += rng.random(hi - lo, dtype=np.float32) - 0.5   # random sub-pixel offset')
            add('        vv += rng.random(hi - lo, dtype=np.float32) - 0.5')
        if rest:
            add('        slew = SLEW_RATE / 86400.0 / ARCSEC_PER_LOWRES_PIXEL * 2.0 * time[lo:hi]   # drift, hi-res pixels')
            add('        direction = math.radians(SLEW_ANGLE) - math.radians(ORIENTAT)')
            add('        uu = uu + (-slew * math.sin(direction)).astype(np.float32)')
            add('        vv = vv + (slew * math.cos(direction)).astype(np.float32)')
            if frame == 'rest_north':
                add('        orient = math.radians(ORIENTAT)   # rotate so north is up')
                add('        du, dv = uu - DETECTOR_CENTRE, vv - DETECTOR_CENTRE')
                add('        uu = du * math.cos(orient) - dv * math.sin(orient) + DETECTOR_CENTRE')
                add('        vv = du * math.sin(orient) + dv * math.cos(orient) + DETECTOR_CENTRE')
        add('        u[lo:hi], v[lo:hi] = uu, vv')
        add('    return u, v')
    if edge:
        flags.append('EDGE_OK')
        add('')
        add('')
        add('# ---- detector edge -----------------------------------------------------------------')
        add('# With its exposure-time correction on, the Explorer drops events near the detector edges, as far from')
        add(f'# each edge as: {EDGE_TRIM_TEXT} pixels (the image frame\'s left and right are AXIS1, bottom and top AXIS2).')
        add(f'EDGE_X_LOW, EDGE_X_HIGH = {EDGE_TRIM.left + 1}, {2048 - EDGE_TRIM.right}   # AXIS1 values to keep')
        add(f'EDGE_Y_LOW, EDGE_Y_HIGH = {EDGE_TRIM.bottom + 1}, {2048 - EDGE_TRIM.top}   # AXIS2 values to keep')

    # ---- what goes into the headers ------------------------------------------------------------
    history = []
    if time_sel is not None:
        history.append(f'SEL_TIME: True for {time_sel[0]:.{TIME_DECIMALS}f} <= TIME < {time_sel[1]:.{TIME_DECIMALS}f} s '
                       '(seconds since the exposure start)')
    if placeable:
        history.append(f'SEL_REGION: True for events inside this region of the {label}: {region.describe()}')
        history.append('Region, unrounded: ' + _region_exact(region))
    if edge:
        history.append(f'EDGE_OK: True for events outside the detector edge trim ({EDGE_TRIM_TEXT} pixels)')
    history.append('Corrections run: ' + ('; '.join(ran) if ran else 'none'))
    history.append('Not applied: flat-field weights, exposure-time map')
    history.append('Existing columns keep their values; scaled ones (TIME) are stored as plain floating point')
    add('')
    add('')
    add('# ---- notes written into the output file ------------------------------------------------')
    add('HISTORY = [')
    add("    'Columns added by the TIME-TAG Explorer selection script v' + SCRIPT_VERSION + ', run '")
    add("    + datetime.now(timezone.utc).strftime('%Y-%m-%d') + ' UTC',")
    add("    'Tool: ' + TOOL_URL,")
    add("    'Not an STScI calibrated product: the raw _tag table plus the flag columns',")
    for text in history:
        add(f'    {text!r},')
    add(']')
    add('')
    add('')
    add('def download_tag_file(rootname, path):')
    add('    """Fetch the TIME-TAG file from MAST (public HST data) into the current folder."""')
    add('    from astroquery.mast import Observations   # pip install astroquery')
    add('    root = rootname.lower()')
    add("    part = path + '.part'")
    add('    try:')
    add('        # HST products have a MAST URI of this form. Observations.download_file exists in every')
    add('        # astroquery release, unlike the newer MastMissions.download_file.')
    add("        result = Observations.download_file(f'mast:HST/product/{root}_tag.fits', local_path=part,")
    add('                                            cache=False)')
    add("        status = result[0] if isinstance(result, (tuple, list)) else result")
    add("        if status == 'ERROR' or not os.path.exists(part):")
    add("            raise RuntimeError(f'MAST download of {root}_tag.fits failed: {result}')")
    add('        os.replace(part, path)')
    add('    finally:')
    add('        if os.path.exists(part):')
    add('            os.remove(part)')
    add('')
    add('')
    add("STRUCTURE_KEYWORD = re.compile(r'T(TYPE|FORM|UNIT|NULL|SCAL|ZERO|DISP|BCOL|DIM)\\d+')   # astropy writes these afresh")
    add("FITS_FORMATS = {'f8': 'D', 'f4': 'E', 'i2': 'I', 'i4': 'J', 'i8': 'K', 'u1': 'B', 'b1': 'L'}")
    add('# Per-column coordinate keywords (TCTYPn, TCUNIn, TCRPXn, TCRVLn, TCDLTn, TRPOSn). astropy keeps these as')
    add('# attributes of the column and writes the header from them, so a new column has to be given them.')
    add("COORDINATE_ATTRIBUTES = ('coord_type', 'coord_unit', 'coord_ref_point', 'coord_ref_value', 'coord_inc', 'time_ref_pos')")
    add('')
    add('')
    add('def stored_column(hdu, name):')
    add('    """One column of the TAG table as a new, unscaled column.')
    add('')
    add('    Columns such as TIME are stored as scaled integers (TSCAL / TZERO). astropy hands out their')
    add('    real values, so they are written back as plain floating point: re-using the old column')
    add('    definition would make astropy try to scale the values back, and fail.')
    add('    """')
    add('    column = hdu.columns[name]')
    add('    values = np.asarray(hdu.data[name])')
    add("    code = FITS_FORMATS.get(f'{values.dtype.kind}{values.dtype.itemsize}')")
    add('    if code is None or values.ndim != 1:')
    add('        return column   # (not the case for TAG files) anything unusual keeps its own definition')
    add('    coordinates = {attribute: getattr(column, attribute) for attribute in COORDINATE_ATTRIBUTES}')
    add('    return fits.Column(name=name, format=code, unit=column.unit, array=values, **coordinates)')
    add('')
    add('')
    add('def in_file_order(flag, order):')
    add('    """Flags are worked out in time order; put them back in the file\'s row order."""')
    add('    if order is None:')
    add('        return flag')
    add('    restored = np.empty_like(flag)')
    add('    restored[order] = flag')
    add('    return restored')
    add('')
    add('')
    add('# ---- download, read, flag, write -------------------------------------------------------')
    add('if os.path.exists(OUTPUT) and not OVERWRITE:')
    add("    raise SystemExit(f'{OUTPUT} already exists: set OVERWRITE = True to replace it')")
    add('if not os.path.exists(TAG_FILE):')
    add("    print(f'Downloading {TAG_FILE} from MAST ...')")
    add('    download_tag_file(ROOTNAME, TAG_FILE)')
    add('')
    add('with fits.open(TAG_FILE) as hdul:')
    add('    events = hdul[1].data')
    add("    time = np.asarray(events['TIME'], dtype=np.float32)   # the Explorer keeps time as float32")
    add("    x = np.asarray(events['AXIS1'], dtype=np.int16)")
    add("    y = np.asarray(events['AXIS2'], dtype=np.int16)")
    add('    order = None')
    add('    if np.any(time[1:] < time[:-1]):   # the Explorer works on time-ordered events (TAG files usually are)')
    add("        order = np.argsort(time, kind='stable')")
    add('        time, x, y = time[order], x[order], y[order]')
    add('')
    add('    flags = {}')
    if time_sel is not None:
        add('    lo, hi = np.searchsorted(time, TIME_WINDOW)')
        add('    sel_time = np.zeros(len(time), dtype=bool)')
        add('    sel_time[lo:hi] = True')
        add("    flags['SEL_TIME'] = in_file_order(sel_time, order)")
    if placeable:
        add('    u, v = frame_coordinates(x, y, time)')
        add("    flags['SEL_REGION'] = in_file_order(in_region(u, v), order)")
    elif region is not None:
        add('    # No region flag: this script cannot place events in the ' + label + '.')
        add('    # With your own coordinates u, v in that frame:  flags[\'SEL_REGION\'] = in_region(u, v)')
    if edge:
        add("    flags['EDGE_OK'] = in_file_order(")
        add('        (x >= EDGE_X_LOW) & (x <= EDGE_X_HIGH) & (y >= EDGE_Y_LOW) & (y <= EDGE_Y_HIGH), order)')
    add('')
    add("    print(f'{len(time):,} events in {TAG_FILE}')")
    for name in flags:
        add(f"    print(f'  {name:<12}{{flags[\"{name}\"].sum():>12,}}')")
    if len(flags) > 1:
        add('    selected = np.logical_and.reduce(list(flags.values()))')
        add("    print(f'  {\"all of them\":<12}{selected.sum():>12,}')")
    add('')
    add("    # The TAG file's own columns (their real values, scaling applied) plus the flag columns.")
    add('    columns = [stored_column(hdul[1], name) for name in hdul[1].columns.names]')
    add("    columns += [fits.Column(name=name, format='L', array=values) for name, values in flags.items()]")
    add('    header = hdul[1].header.copy()')
    add('    for key in [key for key in header if STRUCTURE_KEYWORD.fullmatch(key)]:')
    add('        del header[key]   # the old column definitions, TSCAL and TZERO included')
    add('    table = fits.BinTableHDU.from_columns(fits.ColDefs(columns), header=header)')
    add('    header = table.header')
    add("    header['SELVER'] = (SCRIPT_VERSION, 'TIME-TAG Explorer script version')")
    add("    header['SELURL'] = (TOOL_URL, 'tool URL')")
    if time_sel is not None:
        add(f"    header['SELTMIN'] = ({float(time_sel[0])!r}, '[s] SEL_TIME window start')")
        add(f"    header['SELTMAX'] = ({float(time_sel[1])!r}, '[s] SEL_TIME window end')")
    if placeable:
        add(f"    header['SELREG'] = ({region.kind!r}, 'SEL_REGION shape')")
        add(f"    header['SELFRAME'] = ({frame!r}, 'SEL_REGION reference frame')")
        add(f"    header['SELRAND'] = ({jitter!r}, 'random sub-pixel offset used for SEL_REGION')")
        add(f"    header['SELGEO'] = ({geo!r}, 'geometric distortion used for SEL_REGION')")
    if edge:
        for key, side, pixels in (('SELEDGEL', 'left', EDGE_TRIM.left), ('SELEDGER', 'right', EDGE_TRIM.right),
                                  ('SELEDGEB', 'bottom', EDGE_TRIM.bottom), ('SELEDGET', 'top', EDGE_TRIM.top)):
            add(f"    header['{key}'] = ({pixels}, '[pix] EDGE_OK trim of the {side} edge')")
    add('    for text in HISTORY:')
    add('        for line in textwrap.wrap(text, 68):')
    add('            header.add_history(line)')
    add("    hdul[0].header.add_history('Flag columns added to extension 1; see its header (SELVER, HISTORY).')")
    add('    hdul[1] = table')
    add("    part = OUTPUT + '.part'   # written aside, so a failure never leaves a broken OUTPUT behind")
    add('    try:')
    add('        hdul.writeto(part, overwrite=True)')
    add('        os.replace(part, OUTPUT)')
    add('    finally:')
    add('        if os.path.exists(part):')
    add('            os.remove(part)')
    add("print(f'Wrote {OUTPUT}')")
    # A commented-out example of using the output: the selected events as an image. It names the flag
    # columns this script wrote, so it works as it stands (a flag that was not written is not listed).
    add('')
    add('')
    add('# ---- Example (commented out): the selected events as an image --------------------------------')
    add('# Needs matplotlib too:  pip install matplotlib')
    add('#')
    add('# from astropy.table import Table')
    add('# from matplotlib import pyplot as plt')
    add('#')
    add('# t = Table.read(OUTPUT, hdu=1)')
    if flags:
        add('# w = np.logical_and.reduce([t[name] for name in ' + repr(tuple(flags)) + '])   # events with every flag set')
    else:
        add('# w = np.ones(len(t), dtype=bool)   # this script wrote no flag columns, so every event')
    add('# BIN = 2    # hi-res pixels per image pixel: 1 gives 2048 x 2048, 2 gives the 1024 x 1024 low-res image')
    add('# n = 2048 // BIN')
    add('# edge = 0.5 + n * BIN   # AXIS1 and AXIS2 are 1-indexed: pixel k covers k - 0.5 to k + 0.5')
    add('# image, _, _ = np.histogram2d(t[\'AXIS2\'][w], t[\'AXIS1\'][w],   # rows are AXIS2 (y), columns AXIS1 (x)')
    add('#                              bins=n, range=[[0.5, edge], [0.5, edge]])   # weights=... takes one value per selected event')
    add("# plt.imshow(image, origin='lower', extent=(0.5, edge, 0.5, edge), cmap='inferno')")
    add("# plt.colorbar(label='counts'); plt.xlabel('AXIS1 (hi-res pixel)'); plt.ylabel('AXIS2 (hi-res pixel)')")
    add('# plt.show()')
    return SelectionScript(
        code='\n'.join(L) + '\n', filename=f'select_{root}.py', packages=tuple(packages), pip_command=pip,
        notes=tuple(notes), warnings=tuple(warnings))
