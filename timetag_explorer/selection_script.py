"""The "Selection code" feature: writes a short Python script that reproduces a user's time and region selections on
the original TIME-TAG file and saves the result as columns in a copy of it.

The script holds only the choices (time window, region, reference frame, corrections) and a call to
`timetag_explorer.tagsel.make_tagsel`, which does the work with this package's own code. The reader installs the
same version of the package that the Explorer ran (`package_install_spec`).

Raise SELECTION_SCRIPT_VERSION (constants.py) whenever `build_selection_script` or `tagsel` changes what the script
computes or writes."""
from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import metadata

from . import __version__, config
from .constants import DEFAULT_TOOL_URL, EDGE_TRIM, RANDOM_SEED, SELECTION_SCRIPT_VERSION, TIME_DECIMALS, TimeSel
from .dataset import TagDataset
from .frames import Corrections, frame_label, is_map
from .regions import Region
from .tagsel import EDGE_TRIM_TEXT

PACKAGE = 'timetag-explorer'


@dataclass(frozen=True)
class SelectionScript:
    """A runnable Python script for the current selections, and what its reader needs to know."""
    code: str
    filename: str
    packages: tuple[tuple[str, str], ...]  # (package, what it is for)
    pip_command: str                       # installs the exact version of the science code the Explorer ran
    notes: tuple[str, ...]      # what the script reproduces and what it leaves out
    warnings: tuple[str, ...]   # things that will not work as the reader might expect


def tool_url() -> str:
    """Web address of this tool, for the scripts and FITS headers it writes."""
    return config.PUBLIC_URL or DEFAULT_TOOL_URL


def package_install_spec() -> tuple[str, str | None]:
    """(what to give pip, a warning or None) to install the copy of this package that is running.

    A release from PyPI is `timetag-explorer==X.Y.Z`. A copy installed from a git branch, a commit or an archive
    address (a feature branch the web app is testing) is installed from that same address, which pip records.
    A local development copy cannot be reproduced from outside: the version is given, with a warning.
    """
    version_spec = f'{PACKAGE}=={__version__}'
    try:
        raw = metadata.distribution(PACKAGE).read_text('direct_url.json')
    except metadata.PackageNotFoundError:
        raw = None
    if not raw:
        return version_spec, None
    try:
        info = json.loads(raw)
    except ValueError:
        return version_spec, None
    url = str(info.get('url', ''))
    vcs = info.get('vcs_info') or {}
    if url.startswith(('http://', 'https://')):
        if vcs.get('vcs') == 'git' and vcs.get('commit_id'):
            return f'{PACKAGE} @ git+{url}@{vcs["commit_id"]}', None
        if 'archive_info' in info:
            return f'{PACKAGE} @ {url}', None
    return version_spec, (f'The Explorer is running a development copy of the science code (version {__version__}) that '
                          'is not published, so the exact code cannot be installed from the instruction below; the '
                          'version number may match a release that differs from it.')


def pip_command(spec: str) -> str:
    return f'pip install "{spec}"' if (' ' in spec or '@' in spec) else f'pip install {spec}'


def _corrections_source(ds: TagDataset, corr: Corrections) -> str:
    """`Corrections(...)` with the switches that matter for a script (the planet-map offsets do not)."""
    args = [f'randcorr={bool(corr.randcorr)}', f'geocorr={bool(corr.geocorr)}', f'flatcorr={bool(corr.flatcorr)}',
            f'timecorr={bool(corr.timecorr)}']
    if corr.vigcorr and ds.vignette_applies():
        args.append('vigcorr=True')
        for name in ('vigbot_start', 'vigbot_end', 'vigbot_loss', 'vigtop_start', 'vigtop_end', 'vigtop_loss'):
            args.append(f'{name}={float(getattr(corr, name))!r}')
    return 'Corrections(' + ', '.join(args) + ')'


def build_selection_script(ds: TagDataset, time_sel: TimeSel | None, region: Region | None,
                           frame: str, corr: Corrections) -> SelectionScript:
    """A short Python script that downloads the exposure's TAG file, flags the events the Explorer has selected,
    adds the corrected positions and weights, and writes the file back out (`<rootname>_tagsel.fits`).

    The script calls this package, so it needs the same version of it as the Explorer is running (the script says
    which, and how to install it, and refuses to run with another unless told to). The region is tested at each
    event's own position, not on the Explorer's 4-pixel bins, so events within about two pixels of its edge can differ
    from the highlighted bins. Regions drawn on planet maps (experimental) are not carried over.
    Change what the script computes or writes only together with SELECTION_SCRIPT_VERSION.
    """
    root = (ds.rootname or ds.name.split('_')[0]).lower()
    label = frame_label(ds.frame_info, frame)
    url = tool_url()
    version = SELECTION_SCRIPT_VERSION
    spec, spec_warning = package_install_spec()
    install = pip_command(spec)
    t_text = (f'{time_sel[0]:.{TIME_DECIMALS}f} to {time_sel[1]:.{TIME_DECIMALS}f} s' if time_sel is not None else '')
    output = f'{root}_tagsel.fits'

    warnings = [spec_warning] if spec_warning else []
    if region is not None and is_map(frame):
        warnings.append(f'The region was drawn on a planet map ({label}), which is experimental and not carried over: '
                        'the script contains the time selection and the other columns, but no region flag.')
        region = None

    # which columns the script will write (the frames are found at run time; they depend on the files)
    flags = (['SEL_TIME'] if time_sel is not None else []) + (['SEL_REGION'] if region is not None else [])
    flags += ['EDGE_OK'] if corr.timecorr and any(EDGE_TRIM) else []
    corrected = bool(corr.randcorr or corr.geocorr)
    weights = []
    if corr.flatcorr and ds.weight is not None:
        weights.append('FLATW')
    if corr.vigcorr and ds.vignette_applies():
        weights.append('VIGW')
    if weights:
        weights.append('WEIGHT')
    if corr.flatcorr and ds.weight is None:
        warnings.append(f'No flat-field weights: {ds.flat_problem or "the flat field is not available here"}. The script '
                        'downloads the flat field itself, so FLATW and WEIGHT are written if it can.')

    notes = [
        f'Install the science code first: `{install}`. That is the version the Explorer ran ({__version__}); the script '
        'stops if another is installed, because other versions can give different numbers.',
        f'Run it in an empty folder: it downloads `{root}_tag.fits` and `{root}_spt.fits` (pointing, for the rest-frame '
        f'columns) from MAST and the flat field from the Explorer’s site if they are not there (astroquery and a '
        f'network connection), then writes `{output}`: the TAG file’s own tables with columns added. That is not an '
        'STScI calibrated product, and the headers say so. It will not replace an existing file unless '
        '`OVERWRITE = True`.',
        'AXIS1 and AXIS2 are unchanged (the recorded integer pixels).']
    if corrected:
        notes.append('`XCORR`, `YCORR`: the image-frame position with ' + ' and '.join(
            name for name, on in (('the random sub-pixel offset', corr.randcorr), ('the geometric distortion correction',
                                                                                  corr.geocorr)) if on)
            + ' applied (hi-res pixels counted from 1, like AXIS1 and AXIS2).')
    notes.append('One pair of columns per reference frame the exposure supports, from those positions: `XREST`/`YREST` '
                 '(target rest frame), `XRESTN`/`YRESTN` (rest frame, north up) and, for moon observations, '
                 '`XPLAN`/`YPLAN`, `XPLANN`/`YPLANN` (planet rest frame, and north up). The planet frames need JPL '
                 'Horizons ephemerides (a network connection). Planet maps are not included.')
    if weights:
        notes.append('Weights: ' + ', '.join(f'`{w}`' for w in weights) + '. Sum `WEIGHT` over events for counts and '
                     'divide by the exposure for rates; `FLATW` is the flat field at the recorded pixel, `VIGW` the '
                     'experimental vignetting model (F25SRF2 only).')
    if time_sel is not None:
        notes.append('The time cut compares the Explorer’s float32 copy of TIME, so it picks the same events as the plots.')
    if region is not None:
        notes.append(f'The region is tested at each event’s own position in the {label}. The Explorer tests it on '
                     f'4-pixel bins (an event is inside when its bin’s centre is), so events within about two pixels '
                     'of the edge can differ from the highlighted bins.')
        if corr.randcorr or frame != 'image':
            notes.append('The random sub-pixel offsets are drawn with the same seeded NumPy generator as the Explorer '
                         f'(seed {RANDOM_SEED}), so positions agree with it. That relies on NumPy giving the same '
                         'stream, which it does in practice but does not formally promise between versions.')
    if flags and 'EDGE_OK' in flags:
        notes.append(f'`EDGE_OK` is the Explorer’s detector-edge trim ({EDGE_TRIM_TEXT} pixels); it is not applied to '
                     'the other columns, which cover every event.')
    notes.append('TIME is stored in the TAG file as a scaled integer; the output holds the same values as plain 64-bit '
                 'floating point, so it is a little larger. Other original columns are unchanged.')
    notes.append(f'Script version {version}. It is written into the output file’s header (with the package version), so '
                 'keep the numbers if you report results.')

    # ---- the script ----------------------------------------------------------------------------
    L: list[str] = ['#!/usr/bin/env python']
    add = L.append
    corr_text = ', '.join(f'{name} {"on" if on else "off"}' for name, on in (
        ('flat field', corr.flatcorr), ('random sub-pixel offset', corr.randcorr),
        ('geometric distortion', corr.geocorr), ('exposure-time (edge trim)', corr.timecorr))
        + ((('vignetting model', corr.vigcorr),) if ds.vignette_applies() else ()))
    add(f'"""Reproduce selections made in TIME-TAG Explorer on HST/STIS exposure {root}.')
    add('')
    add(f'TIME-TAG Explorer: {url}')
    add(f'Script version {version}, written {datetime.now(timezone.utc):%Y-%m-%d} (UTC), for timetag-explorer {__version__}.')
    add('')
    add('Selections in the Explorer')
    if time_sel is not None:
        add(f'  time         {t_text}')
    if region is not None:
        add(f'  region       {region.describe()}  ({label})')
    add(f'  corrections  {corr_text}')
    add('')
    add('Install the science code (the version the Explorer ran), then run this script')
    add(f'  {install}')
    add(f'  python select_{root}.py')
    add('')
    add('What it does')
    add(f'  1. Downloads {root}_tag.fits and {root}_spt.fits (MAST) and the flat field, unless they are here.')
    add('  2. Flags the selected events (True = selected) and adds the corrected positions and weights.')
    add(f'  3. Writes {output}: the TAG file with those columns added, and header keywords and')
    add('     HISTORY that record the script and package versions and the corrections that ran.')
    add('"""')
    add('import timetag_explorer')
    add('from timetag_explorer.frames import Corrections')
    if region is not None:
        add(f'from timetag_explorer.regions import {type(region).__name__}')
    add('from timetag_explorer.tagsel import make_tagsel')
    add('')
    add(f"SCRIPT_VERSION = {version!r}   # version of the TIME-TAG Explorer selection script")
    add(f"PACKAGE_VERSION = {__version__!r}   # the timetag-explorer version the Explorer ran")
    add(f'INSTALL = {install!r}')
    add(f'TOOL_URL = {url!r}')
    add(f'ROOTNAME = {root!r}')
    add("FOLDER = '.'                  # where the files are (downloaded here if missing) and the output is written")
    add('OVERWRITE = False             # True: replace the output if it exists already')
    add('ALLOW_OTHER_VERSION = False   # True: run with another timetag-explorer version (the numbers may differ)')
    add('')
    add('# ---- the selections ------------------------------------------------------------------')
    if time_sel is not None:
        add(f'TIME_WINDOW = ({float(time_sel[0])!r}, {float(time_sel[1])!r})   # seconds since the exposure start: start <= TIME < end')
    else:
        add('TIME_WINDOW = None')
    if region is not None:
        add(f'# Drawn in the {label}; coordinates are 0-indexed hi-res detector pixels (AXIS1 and AXIS2 count from 1).')
        add(f'REGION = {region!r}')
    else:
        add('REGION = None')
    add(f'REGION_FRAME = {frame if region is not None else "image"!r}')
    add('')
    add('# ---- the corrections, as the Explorer had them ---------------------------------------')
    add(f'CORRECTIONS = {_corrections_source(ds, corr)}')
    add('')
    add('')
    add('if timetag_explorer.__version__ != PACKAGE_VERSION and not ALLOW_OTHER_VERSION:')
    add("    raise SystemExit(f'This script needs timetag-explorer {PACKAGE_VERSION} (found '")
    add("                     f'{timetag_explorer.__version__}). Install it with:  {INSTALL}')")
    add('')
    add('OUTPUT = make_tagsel(')
    add('    ROOTNAME, time_sel=TIME_WINDOW, region=REGION, frame=REGION_FRAME, corr=CORRECTIONS,')
    add('    script_version=SCRIPT_VERSION, tool_url=TOOL_URL, folder=FOLDER, overwrite=OVERWRITE)')
    # A commented-out example of using the output: an image from the added columns. It names columns this script
    # writes, so it works as it stands.
    add('')
    add('')
    add('# ---- Example (commented out): the selected events as an image --------------------------------')
    add('# Needs matplotlib:  pip install matplotlib')
    add('#')
    add('# import numpy as np')
    add('# from astropy.table import Table')
    add('# from matplotlib import pyplot as plt')
    add('#')
    add('# t = Table.read(OUTPUT, hdu=1)')
    if flags:
        add('# keep = np.logical_and.reduce([t[name] for name in ' + repr(tuple(flags)) + '])   # every flag set')
    else:
        add('# keep = np.ones(len(t), dtype=bool)   # this script wrote no flag columns, so every event')
    x, y = ('XCORR', 'YCORR') if corrected else ('AXIS1', 'AXIS2')
    add(f"# x, y = t['{x}'][keep], t['{y}'][keep]   # hi-res pixels, counted from 1")
    add("# w = t['WEIGHT'][keep]" if weights else '# w = None   # no weights were written')
    add('# BIN = 2    # hi-res pixels per image pixel: 1 gives 2048 x 2048, 2 gives the 1024 x 1024 low-res image')
    add('# n = 2048 // BIN')
    add('# edge = 0.5 + n * BIN   # pixel k covers k - 0.5 to k + 0.5')
    add('# image, _, _ = np.histogram2d(y, x, bins=n, range=[[0.5, edge], [0.5, edge]], weights=w)   # rows are y')
    add("# plt.imshow(image, origin='lower', extent=(0.5, edge, 0.5, edge), cmap='inferno')")
    add("# plt.colorbar(label='counts'); plt.xlabel('x (hi-res pixel)'); plt.ylabel('y (hi-res pixel)')")
    add('# plt.show()')
    packages = (('timetag-explorer', 'the Explorer’s own code: corrections, reference frames, planet frames; '
                                     'it brings numpy, astropy, matplotlib and astroquery'),)
    return SelectionScript(code='\n'.join(L) + '\n', filename=f'select_{root}.py', packages=packages,
                           pip_command=install, notes=tuple(notes), warnings=tuple(warnings))
