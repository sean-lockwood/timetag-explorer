"""The "Selection code" feature: writes a short Python script that reproduces a user's named time ranges and regions on
the original TIME-TAG file and saves the result as columns in a copy of it.

The script holds only the choices (time ranges, regions and their frames, corrections) and a call to
`timetag_explorer.tagsel.make_tagsel`, which does the work with this package's own code. The reader installs the
same version of the package that the Explorer ran (`package_install_spec`).

Raise SELECTION_SCRIPT_VERSION (constants.py) whenever `build_selection_script` or `tagsel` changes what the script
computes or writes."""
from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from . import __version__, config
from .constants import DEFAULT_TOOL_URL, EDGE_TRIM, RANDOM_SEED, SELECTION_SCRIPT_VERSION, TIME_DECIMALS, TimeSel
from .dataset import TagDataset
from .frames import Corrections, frame_label, is_map
from .regions import Region, region_source
from .selections import RegionTable, TimeTable, single_selections
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


def package_version_label() -> str:
    """Which copy of this package is running, for a page footer: `1.5.0` for a release, `1.5.0 (branch @ abc1234)` for a
    copy installed from git (the branch or tag pip was given, and the commit), `1.5.0 (archive)` for an archive address,
    and for a development copy `1.5.0 (dev, branch @ abc1234)` when its folder is a git checkout, else `1.5.0 (dev)`."""
    try:
        raw = metadata.distribution(PACKAGE).read_text('direct_url.json')
        info = json.loads(raw) if raw else None
    except (metadata.PackageNotFoundError, ValueError):
        info = None
    if not info:
        return __version__
    vcs = info.get('vcs_info') or {}
    if str(info.get('url', '')).startswith(('http://', 'https://')):
        if vcs.get('vcs') == 'git' and vcs.get('commit_id'):
            ref = vcs.get('requested_revision')
            commit = vcs['commit_id'][:7]
            return f'{__version__} ({ref} @ {commit})' if ref and not str(ref).startswith(commit) else f'{__version__} (@ {commit})'
        if 'archive_info' in info:
            return f'{__version__} (archive)'
    return f'{__version__} (dev{_checkout_label()})'


def _checkout_label() -> str:
    """`, branch @ abc1234` when this package's folder is a git checkout (a development copy), else nothing."""
    folder = Path(__file__).resolve().parent
    try:
        out = [subprocess.run(['git', '-C', str(folder), *args], capture_output=True, text=True, timeout=5,
                              check=True).stdout.strip() for args in (('rev-parse', '--abbrev-ref', 'HEAD'),
                                                                       ('rev-parse', '--short=7', 'HEAD'))]
    except (OSError, subprocess.SubprocessError):
        return ''
    return f', {out[0]} @ {out[1]}' if all(out) else ''


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


def _selection_source(row, kind: str) -> str:
    """One `TimeRange(...)` or `RegionSelection(...)` line of the script; empty note and default flags are left out."""
    parts = [f'name={row.name!r}']
    if row.derived:
        parts.append(f'expr={row.expr!r}')
    elif kind == 'time':
        parts += [f'start={float(row.start)!r}', f'stop={float(row.stop)!r}']
    else:
        parts += [f'frame={row.frame!r}', f'region={region_source(row.region)}']
    if not row.derived:
        parts.append(f'bit={row.bit}')
    parts.append(f'color={row.color!r}')
    if not row.visible:
        parts.append('visible=False')
    if row.note:
        parts.append(f'note={row.note!r}')
    cls = 'TimeRange' if kind == 'time' else 'RegionSelection'
    return f'    {cls}({", ".join(parts)}),'


def build_selection_script(ds: TagDataset, time_sel: TimeSel | None, region: Region | None,
                           frame: str, corr: Corrections) -> SelectionScript:
    """The script for the Explorer's current selection: one time window and one region, the region as drawn on the
    Explorer's 0-indexed arrays (it is moved to the file's coordinates, hi-res pixels counted from 1). The time range is
    named `T1` and the region `R1`. See `build_selections_script` for any number of named selections.
    Regions drawn on planet maps (not supported yet) are not carried over, and a warning says so."""
    warnings = []
    if region is not None and is_map(frame):
        warnings.append(f'The region was drawn on a planet map ({frame_label(ds.frame_info, frame)}), which is not '
                        'supported yet: the script contains the time selection and the other columns, but no region flag.')
        region = None
    times, regions = single_selections(time_sel, region, frame)
    script = build_selections_script(ds, corr, times, regions)
    return replace(script, warnings=tuple(warnings) + script.warnings)


def build_selections_script(ds: TagDataset, corr: Corrections, times: TimeTable | None = None,
                            regions: RegionTable | None = None) -> SelectionScript:
    """A short Python script that downloads the exposure's TAG file, flags the events in each named time range and region,
    adds the corrected positions and weights, and writes the file back out (`<rootname>_tagsel.fits`).

    `times` and `regions` are `TimeTable` / `RegionTable` (regions in file coordinates, hi-res pixels counted from 1).
    The script calls this package, so it needs the same version of it as the Explorer is running (the script says which,
    and how to install it, and refuses to run with another unless told to). Each region is tested at each event's own
    position, not on the Explorer's 4-pixel bins, so events within about two pixels of its edge can differ from the
    highlighted bins. Change what the script computes or writes only together with SELECTION_SCRIPT_VERSION.
    """
    times = times if times is not None else TimeTable()
    regions = regions if regions is not None else RegionTable()
    root = (ds.rootname or ds.name.split('_')[0]).lower()
    url = tool_url()
    version = SELECTION_SCRIPT_VERSION
    spec, spec_warning = package_install_spec()
    install = pip_command(spec)
    output = f'{root}_tagsel.fits'

    warnings = [spec_warning] if spec_warning else []

    # which columns the script will write (the frames are found at run time; they depend on the files)
    flags = (['SEL_TIME'] if times else []) + (['SEL_REGION'] if regions else [])
    flags += ['EDGE_OK'] if corr.timecorr and any(EDGE_TRIM) else []
    corrected = bool(corr.randcorr or corr.geocorr)
    weights = []
    if corr.flatcorr and ds.weight is not None:
        weights.append('FLATW')
    if corr.vigcorr and ds.vignette_applies():
        weights.append('VIGW')
    if corr.flatcorr and ds.weight is None:
        warnings.append(f'No flat-field weights: {ds.flat_problem or "the flat field is not available here"}. The script '
                        'downloads the flat field itself, so FLATW is written if it can.')
    many = len(times.primitives) > 1 or len(regions.primitives) > 1
    labels = {row.name: frame_label(ds.frame_info, row.frame) for row in regions.primitives}

    notes = [
        f'Install the science code first: `{install}`. That is the version the Explorer ran ({__version__}); the script '
        'stops if another is installed, because other versions can give different numbers.',
        f'Run it in an empty folder: it downloads `{root}_tag.fits` and `{root}_spt.fits` (pointing, for the rest-frame '
        f'columns) from MAST and the flat field from the Explorer’s site if they are not there (astroquery and a '
        f'network connection), then writes `{output}`: the TAG file’s own tables with columns added. That is not an '
        'STScI calibrated product, and the headers say so. It will not replace an existing file unless '
        '`OVERWRITE = True`.',
        'AXIS1 and AXIS2 are unchanged (the recorded integer pixels).']
    if flags[:2] and flags[0] in ('SEL_TIME', 'SEL_REGION'):
        notes.append('`SEL_TIME` and `SEL_REGION` hold one bit per time range or region (16-bit words): bit 0 is the first '
                     'one listed in the script, bit 1 the next, and so on, and an event inside several has several bits '
                     'set. The `TIMESEL` and `REGSEL` extensions list them with their names, shapes, frames, notes and '
                     'colors, and `read_selections(OUTPUT)` picks events by name (any case) or bit position, for example '
                     '`sel.in_region(t, "R1")`. A selection written as an expression of others (such as `oval - io`) '
                     'has no bit of its own.')
    if corrected:
        notes.append('`XCORR`, `YCORR`: the image-frame position with ' + ' and '.join(
            name for name, on in (('the random sub-pixel offset', corr.randcorr), ('the geometric distortion correction',
                                                                                  corr.geocorr)) if on)
            + ' applied (hi-res pixels counted from 1, like AXIS1 and AXIS2).')
    notes.append('One pair of columns per reference frame the exposure supports, from those positions: `XREST`/`YREST` '
                 '(target rest frame) and, for moon observations, `XPLAN`/`YPLAN` (planet rest frame). The north-up '
                 'versions (`XRESTN`/`YRESTN`, `XPLANN`/`YPLANN`) make the file much larger, so they are left out unless '
                 'you set `INCLUDE_NORTH_UP = True` in the script or a region was drawn in one. The planet frames need '
                 'JPL Horizons ephemerides (a network connection). Planet maps are not supported.')
    if weights:
        notes.append('Weights: ' + ', '.join(f'`{w}`' for w in weights) + '. There is no combined weight column: an '
                     'event counts as the product of the weights that exist (`FLATW` is the flat field at the '
                     'recorded pixel, `VIGW` the experimental vignetting model, F25SRF2 only). The commented examples '
                     'at the end of the script multiply them; sum the product over events for counts and divide by the '
                     'exposure for rates.')
    if times:
        notes.append('The time cut compares the Explorer’s float32 copy of TIME, so it picks the same events as the plots. '
                     'The exposure of a time range is not stored; `sel.time_exposure(name, gti)` works it out from the '
                     'good-time intervals when asked.')
    if regions:
        notes.append('Each region is tested at each event’s own position, on the values written to the position columns of '
                     'its frame, so applying the region from `REGSEL` to those columns by hand gives the same flags. The '
                     'Explorer tests regions on 4-pixel bins (an event is inside when its bin’s centre is), so events '
                     'within about two pixels of the edge can differ from the highlighted bins. A region’s area and its '
                     'exposure (which varies across it) are not stored; the package works them out when asked.')
        notes.append('Every frame a region was drawn in gets its position columns. If one of those frames cannot be built '
                     '(for example without the SPT file), the script stops before it writes anything.')
        if corr.randcorr or any(row.frame != 'image' for row in regions.primitives):
            notes.append('The random sub-pixel offsets are drawn with the same seeded NumPy generator as the Explorer '
                         f'(seed {RANDOM_SEED}), so positions agree with it. That relies on NumPy giving the same '
                         'stream, which it does in practice but does not formally promise between versions.')
    if 'EDGE_OK' in flags:
        notes.append(f'`EDGE_OK` is the Explorer’s detector-edge trim ({EDGE_TRIM_TEXT} pixels); it is not applied to '
                     'the other columns, which cover every event.')
    notes.append('The original columns are unchanged.')
    if not times and not regions:
        notes.append('No time range or region was selected, so the file has no `SEL_TIME` / `SEL_REGION` columns or '
                     '`TIMESEL` / `REGSEL` extensions (`TTETSEL` and `TTERSEL` are False): the corrections, positions and '
                     'weights are still written.')
    notes.append('The primary header ends with one block of this tool’s keywords, under `/ TIMETAG-EXPLORER KEYWORDS`. '
                 '`TTERAND`, `TTEGEO`, `TTEFLAT`, `TTEVIGN`, `TTEEDGE` say which steps ran (`DONE` or `OMIT`; a step is '
                 '`OMIT` when it was off or could not take effect, such as the vignetting model for an aperture other '
                 'than F25SRF2), `TTEFFILE` names the flat-field file used, and `TTETSEL` and `TTERSEL` are True when time '
                 'ranges or regions were selected. They are separate from the calibration pipeline’s own keywords '
                 '(`RANDCORR`, `FLATCORR` ...), which are left as they were, so `stistools` is unaffected.')
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
    if not times and not regions:
        add('  none (the corrections are still applied; the file gets no selection columns)')
    for row in times:
        what = f'= {row.expr}' if row.derived else \
            f'{row.start:.{TIME_DECIMALS}f} to {row.stop:.{TIME_DECIMALS}f} s   (bit {row.bit})'
        add(f'  time    {row.name:<10}{what}')
    for row in regions:
        what = f'= {row.expr}' if row.derived else f'{row.region.describe()}   ({labels[row.name]}, bit {row.bit})'
        add(f'  region  {row.name:<10}{what}')
    add(f'  corrections  {corr_text}')
    add('')
    add('Install the science code (the version the Explorer ran), then run this script')
    add(f'  {install}')
    add(f'  python select_{root}.py')
    add('')
    add('What it does')
    add(f'  1. Downloads {root}_tag.fits and {root}_spt.fits (MAST) and the flat field, unless they are here.')
    add('  2. Flags the events in each time range and region (one bit each in SEL_TIME and SEL_REGION) and adds')
    add('     the corrected positions and weights.')
    add(f'  3. Writes {output}: the TAG file with those columns added, and the TIMESEL and REGSEL tables that')
    add('     name the bits. The primary header records which steps were done (TTERAND, TTEGEO, TTEFLAT, TTEVIGN,')
    add('     TTEEDGE), the flat file (TTEFFILE), whether time or space selections were made (TTETSEL, TTERSEL),')
    add('     and the script and package versions.')
    add('"""')
    add('import timetag_explorer')
    add('from timetag_explorer.frames import Corrections')
    shapes = sorted({type(row.region).__name__ for row in regions.primitives})
    if shapes:
        add(f'from timetag_explorer.regions import {", ".join(shapes)}')
    selection_classes = (['RegionSelection'] if regions else []) + (['TimeRange'] if times else [])
    if selection_classes:
        add(f'from timetag_explorer.selections import {", ".join(sorted(selection_classes))}')
    add('from timetag_explorer.tagsel import make_tagsel, read_selections')
    add('')
    add(f"SCRIPT_VERSION = {version!r}   # version of the TIME-TAG Explorer selection script")
    add(f"PACKAGE_VERSION = {__version__!r}   # the timetag-explorer version the Explorer ran")
    add(f'INSTALL = {install!r}')
    add(f'TOOL_URL = {url!r}')
    add(f'ROOTNAME = {root!r}')
    add("FOLDER = '.'                  # where the files are (downloaded here if missing) and the output is written")
    add('OVERWRITE = False             # True: replace the output if it exists already')
    add('INCLUDE_NORTH_UP = False      # True: also write the north-up frame columns (XRESTN, YRESTN, ...); a larger file')
    add('ALLOW_OTHER_VERSION = False   # True: run with another timetag-explorer version (the numbers may differ)')
    add('')
    add('# ---- the selections ------------------------------------------------------------------')
    add('# Names are case-insensitive. Each shape gets a bit (0, 1, 2, ...) in SEL_TIME or SEL_REGION; a selection')
    add('# written as an expression of earlier names (| union, & intersection, - difference, ~ not) has no bit.')
    if times:
        add('# Time ranges: seconds since the exposure start, start <= TIME < stop.')
        add('TIMES = [')
        L.extend(_selection_source(row, 'time') for row in times)
        add(']')
    else:
        add('TIMES = []')
    if regions:
        add('# Regions: hi-res detector pixels counted from 1, like AXIS1 and the XCORR, XREST ... columns, tested in the')
        add('# reference frame each was drawn in.')
        add('REGIONS = [')
        for row in regions:
            L.append(_selection_source(row, 'region') + (f'   # {labels[row.name]}' if not row.derived else ''))
        add(']')
    else:
        add('REGIONS = []')
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
    add('    ROOTNAME, times=TIMES, regions=REGIONS, corr=CORRECTIONS,')
    add('    script_version=SCRIPT_VERSION, tool_url=TOOL_URL, north_up=INCLUDE_NORTH_UP, folder=FOLDER,')
    add('    overwrite=OVERWRITE)')
    # Commented-out examples of using the output. They name only columns and selections this script writes, so they work
    # as they stand. Both build `keep` one line per selection (`keep &= ...`, so a line can be negated or commented out)
    # and weight events by FLATW x VIGW, whichever of those exist. Rates need the exposure, which comes from the GTI table
    # (and the first time range, if there is one).
    x, y = ('XCORR', 'YCORR') if corrected or any(row.frame == 'image' for row in regions.primitives) else ('AXIS1', 'AXIS2')
    first_time = times.rows[0].name if times else None

    def keep_lines():
        out = ['# keep = np.ones(len(t), dtype=bool)']
        out += [f"# keep &= sel.in_time(t, {row.name!r})" for row in times]
        out += [f"# keep &= sel.in_region(t, {row.name!r})" for row in regions]
        out += ["# keep &= t['EDGE_OK']"] if 'EDGE_OK' in flags else []
        if not flags:
            out[0] += '   # this script wrote no flag columns, so every event'
        elif many:
            out.append('# # (to combine selections, use | and ~ on the masks, or write an expression in the lists above)')
        return out

    def weight_lines():
        out = ['# w = np.ones(len(t))   # an event counts as the product of its weights']
        out += [f"# w *= t['{name}']" for name in weights]
        if not weights:
            out[0] += ' (this script wrote no weights, so each event counts 1)'
        return out

    def window_lines():
        if times:
            return [f"# # the exposure is only the time range {first_time!r} (switch this with its keep line)",
                    f"# t0, t1 = sel.times.window({first_time!r})   # seconds since the exposure start",
                    '# # without the time line, use the whole exposure instead:  t0, t1 = 0.0, gti_rows[:, 1].max()']
        return ['# t0, t1 = 0.0, gti_rows[:, 1].max()   # the whole exposure']

    imports = ['# import warnings', '# import numpy as np', '# from astropy.io import fits',
               '# from astropy.table import Table', '# from astropy.units import UnitsWarning',
               '# from matplotlib import pyplot as plt']
    read = ['# # STIS writes the unit "seconds", which the FITS standard does not know; astropy warns, harmlessly',
            '# with warnings.catch_warnings():',
            "#     warnings.filterwarnings('ignore', category=UnitsWarning, message='.*did not parse as fits unit.*')",
            "#     t = Table.read(OUTPUT, hdu=1)",
            "#     gti = Table.read(OUTPUT, hdu='GTI')",
            "# gti_rows = np.column_stack((gti['START'], gti['STOP']))   # good-time intervals, seconds",
            "# sel = read_selections(OUTPUT)   # the time ranges and regions: by name (any case) or bit position"]
    L.extend(['', '',
              '# ---- Example 1 (commented out): a light curve, weighted counts per second --------------------',
              '# Needs matplotlib:  pip install matplotlib', '#'])
    L.extend(imports)
    L.append('#')
    L.extend(read)
    L.append('#')
    L.extend(keep_lines())
    L.extend(weight_lines())
    L.append('#')
    L.extend(window_lines())
    L.append('# BIN = 3.0   # seconds per light-curve bin')
    L.append('# edges = np.append(np.arange(t0, t1, BIN), t1)   # the last bin is shorter if BIN does not divide t1 - t0')
    L.append('# lo, hi = edges[:-1], edges[1:]')
    L.append("# counts, _ = np.histogram(t['TIME'][keep], bins=edges, weights=w[keep])   # weighted counts per bin")
    L.append("# # exposure per bin: the overlap of each good time interval with the bin")
    L.append("# exposure = np.clip(np.minimum(gti['STOP'][:, None], hi) - np.maximum(gti['START'][:, None], lo),"
             " 0, None).sum(axis=0)")
    L.append("# rate = np.divide(counts, exposure, out=np.full(len(counts), np.nan), where=exposure > 0)   # counts / s")
    L.append("# plt.step(lo, rate, where='post'); plt.xlabel('time since exposure start (s)'); plt.ylabel('count rate (counts / s)')")
    L.append('# plt.show()')
    L.extend(['', '',
              '# ---- Example 2 (commented out): the selected events as an image, in counts per second ----------',
              '# Needs matplotlib:  pip install matplotlib', '#'])
    L.extend(imports)
    L.append('#')
    L.extend(read)
    L.append('#')
    L.extend(keep_lines())
    L.extend(weight_lines())
    L.append('#')
    if times:
        L.append(f"# # the exposure (s) in the time range {first_time!r}: its overlap with the good time intervals")
        L.append(f"# exposure = sel.time_exposure({first_time!r}, gti_rows)")
    else:
        L.append('# # the exposure (s): all the good time')
        L.append("# exposure = float(np.sum(gti_rows[:, 1] - gti_rows[:, 0]))")
    L.append('#')
    L.append(f"# xs, ys = t['{x}'][keep], t['{y}'][keep]   # hi-res pixels, counted from 1")
    L.append('# BIN = 2    # hi-res pixels per image pixel: 1 gives 2048 x 2048, 2 gives the 1024 x 1024 low-res image')
    L.append('# n = 2048 // BIN')
    L.append('# edge = 0.5 + n * BIN   # pixel k covers k - 0.5 to k + 0.5')
    L.append('# counts, _, _ = np.histogram2d(ys, xs, bins=n, range=[[0.5, edge], [0.5, edge]], weights=w[keep])   # rows are y')
    L.append('# image = counts / exposure   # counts / s per image pixel')
    L.append("# plt.imshow(image, origin='lower', extent=(0.5, edge, 0.5, edge), cmap='inferno')")
    L.append("# plt.colorbar(label='counts / s'); plt.xlabel('x (hi-res pixel)'); plt.ylabel('y (hi-res pixel)')")
    L.append('# plt.show()')
    packages = (('timetag-explorer', 'the Explorer’s own code: corrections, reference frames, planet frames; '
                                     'it brings numpy, astropy, matplotlib and astroquery'),)
    return SelectionScript(code='\n'.join(L) + '\n', filename=f'select_{root}.py', packages=packages,
                           pip_command=install, notes=tuple(notes), warnings=tuple(warnings))
