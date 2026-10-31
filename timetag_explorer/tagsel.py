"""The work a selection script does: flag the selected events, add the corrected positions and weights, and write
`<rootname>_tagsel.fits`.

A script written by `build_selection_script` is only the choices (time window, region, frame, corrections) and a call
to `make_tagsel`, so the numbers come from the same code as the Explorer's own images and light curves.

The output is the exposure's TIME-TAG table with columns added; the columns that were there (AXIS1 and AXIS2 included)
are not changed. Added columns, in time-tag row order:

* flags (`L`): `SEL_TIME`, `SEL_REGION` and `EDGE_OK`, for the selections made and the corrections that are on;
* positions (`E`, float32, hi-res pixels counted from 1 like AXIS1 and AXIS2): `XCORR`, `YCORR` when the random
  sub-pixel offset or the geometric distortion correction is on (the image-frame position with both applied), and one
  pair per reference frame the exposure supports, from those corrected positions: `XREST`/`YREST` (target rest
  frame) and `XPLAN`/`YPLAN` (planet rest frame, for moon observations). The north-up frames, `XRESTN`/`YRESTN` and
  `XPLANN`/`YPLANN`, are written only with `north_up=True` (they make the file larger);
* weights (`E`): `FLATW` (flat field) and `VIGW` (vignetting model, F25SRF2 only), for the corrections that are on.
  An event's weight is their product; they are separate columns so either can be left out of a calculation.
"""
from __future__ import annotations

import functools
import os
import re
import textwrap
import urllib.request
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.units import UnitsWarning

from . import config
from .constants import EDGE_TRIM, FLAT_URL, RANDOM_SEED, TIME_DECIMALS, TimeSel
from .dataset import TagDataset, _is_sorted, load_dataset
from .frames import FRAMES, Corrections, frame_label, is_map, is_north, is_planet
from .regions import Region, region_source
from .transforms import frame_coordinates, load_geo_lut

# Columns of the positions in each reference frame (the image frame is XCORR / YCORR).
POSITION_COLUMNS = {'rest': ('XREST', 'YREST'), 'rest_north': ('XRESTN', 'YRESTN'),
                    'planet': ('XPLAN', 'YPLAN'), 'planet_north': ('XPLANN', 'YPLANN')}
CORRECTED_COLUMNS = ('XCORR', 'YCORR')
PIXEL_FRAMES = ('image',) + tuple(POSITION_COLUMNS)    # the frames with position columns; maps are not written
EDGE_TRIM_TEXT = f'left {EDGE_TRIM.left}, right {EDGE_TRIM.right}, bottom {EDGE_TRIM.bottom}, top {EDGE_TRIM.top}'


@dataclass
class TagselColumns:
    """The columns to add, in time order (the dataset's order), and what was left out and why."""
    arrays: dict = field(default_factory=dict)     # name -> array, in the order they are written
    flags: tuple = ()                              # names of the flag columns
    skipped: list = field(default_factory=list)    # (what, reason)
    history: list = field(default_factory=list)    # lines for the FITS HISTORY


def _positions(ds: TagDataset, frame: str, corr: Corrections):
    """(x, y) of every event in a frame as the Explorer computes them: 0-indexed hi-res pixels, float32."""
    track = ds.planet_track() if is_planet(frame) else None
    geo_lut = load_geo_lut() if corr.geocorr else None
    return frame_coordinates(ds.x, ds.y, ds.time, ds.frame_info, frame, corr,
                             seed=RANDOM_SEED, geo_lut=geo_lut, track=track)


def _one_based(a: np.ndarray) -> np.ndarray:
    """0-indexed position to the 1-indexed convention of AXIS1 and AXIS2."""
    return a + np.float32(1.0)


def selection_columns(ds: TagDataset, time_sel: TimeSel | None, region: Region | None, frame: str,
                      corr: Corrections, north_up: bool = False) -> TagselColumns:
    """Work out the columns a selection script adds (time order). Needs ephemerides (the network) only for planet
    frames. A frame that cannot be built is skipped with its reason, unless the region was drawn in it. The north-up
    frames' position columns are written only when `north_up` (a region drawn in one is still tested in it)."""
    out = TagselColumns()
    flags, positions = {}, {}
    n = len(ds.time)
    if time_sel is not None:
        sel = np.zeros(n, dtype=bool)
        sel[ds.time_slice(time_sel)] = True
        flags['SEL_TIME'] = sel
        out.history.append(f'SEL_TIME: True for {time_sel[0]:.{TIME_DECIMALS}f} <= TIME < {time_sel[1]:.{TIME_DECIMALS}f} s '
                           '(seconds since the exposure start)')

    # Positions in every frame that can be built, one frame at a time (only the region's frame is kept for the test).
    # The region's frame is built first, so a failure there is an error, not a gap.
    corrected = bool(corr.randcorr or corr.geocorr)
    todo = list(dict.fromkeys((['image'] if frame == 'image' else [frame] if region is not None else [])
                              + ['image'] + list(POSITION_COLUMNS)))
    for f in todo:
        try:
            if f != 'image':
                ds.prepare_frame(f)
            u, v = _positions(ds, f, corr)
        except Exception as exc:   # an unavailable frame (ValueError), or Horizons not reachable
            if region is not None and f == frame:
                raise RuntimeError(f'the region was drawn in the {frame_label(ds.frame_info, f)}, which cannot be '
                                   f'built: {exc}') from exc
            if north_up or not is_north(f):
                out.skipped.append((frame_label(ds.frame_info, f), str(exc)))
            continue
        if region is not None and f == frame:
            flags['SEL_REGION'] = np.asarray(region.contains(u, v), dtype=bool)
            out.history.append(f'SEL_REGION: True for events inside this region of the {frame_label(ds.frame_info, f)}, '
                               f'tested at each event\'s own position: {region.describe()}')
            out.history.append('Region, unrounded: ' + region_source(region))
        if f == 'image':
            if corrected:
                positions['XCORR'], positions['YCORR'] = _one_based(u), _one_based(v)
                what = ' and '.join(name for name, on in (('random sub-pixel offset', corr.randcorr),
                                                         ('geometric distortion', corr.geocorr)) if on)
                out.history.append(f'XCORR, YCORR: the image-frame position with {what} applied')
        elif north_up or not is_north(f):
            xname, yname = POSITION_COLUMNS[f]
            positions[xname], positions[yname] = _one_based(u), _one_based(v)
            out.history.append(f'{xname}, {yname}: the position in the {frame_label(ds.frame_info, f)}, from the '
                               'corrected positions')
        del u, v
    valid = ds.valid_events(corr)
    if valid is not None:
        flags['EDGE_OK'] = np.array(valid, dtype=bool)
        out.history.append(f'EDGE_OK: True for events outside the detector edge trim ({EDGE_TRIM_TEXT} pixels)')
    arrays = out.arrays
    arrays.update(flags)
    arrays.update(positions)
    out.flags = tuple(flags)

    # weights
    flat = ds.weight if corr.flatcorr and ds.weight is not None else None
    gain = ds.vignette_gain(corr)
    if corr.flatcorr and flat is None:
        out.skipped.append(('flat field', ds.flat_problem or 'not available'))
    if flat is not None:
        arrays['FLATW'] = np.asarray(flat, dtype=np.float32)
        out.history.append('FLATW: flat-field weight at the recorded pixel; an event counts as FLATW x VIGW (where both exist)')
    if gain is not None:
        arrays['VIGW'] = gain
        out.history.append('VIGW: vignetting-model weight (experimental; F25SRF2 only) at the distortion-corrected row')
    return out


# ---- getting the files ---------------------------------------------------------------------------------------------

def download_mast(rootname: str, suffix: str, dest: Path) -> None:
    """Fetch one product of an HST exposure (public data) from MAST."""
    from astroquery.mast import Observations   # download_file exists in every astroquery release
    root = rootname.lower()
    part = dest.with_name(dest.name + '.part')
    try:
        result = Observations.download_file(f'mast:HST/product/{root}_{suffix}.fits', local_path=str(part), cache=False)
        status = result[0] if isinstance(result, (tuple, list)) else result
        if status == 'ERROR' or not part.exists():
            raise RuntimeError(f'MAST download of {root}_{suffix}.fits failed: {result}')
        os.replace(part, dest)
    finally:
        part.unlink(missing_ok=True)


def download_flat(dest: Path) -> None:
    """Fetch the flat field from the address the Explorer publishes it at."""
    part = dest.with_name(dest.name + '.part')
    try:
        urllib.request.urlretrieve(FLAT_URL, part)
        os.replace(part, dest)
    finally:
        part.unlink(missing_ok=True)


def fetch_inputs(rootname: str, folder: Path, corr: Corrections, log=print) -> Path:
    """The TAG file (required), the SPT file (for the rest frames) and the flat field (when that switch is on) in
    `folder`, downloading what is missing. Returns the TAG file's path."""
    root = rootname.lower()
    tag = folder / f'{root}_tag.fits'
    if not tag.exists():
        log(f'Downloading {tag.name} from MAST ...')
        download_mast(root, 'tag', tag)
    spt = folder / f'{root}_spt.fits'
    if not spt.exists() and not (folder / f'{root}_spt.fits.gz').exists():
        try:
            log(f'Downloading {spt.name} from MAST (pointing information, for the rest-frame columns) ...')
            download_mast(root, 'spt', spt)
        except Exception as exc:
            log(f'  no SPT file ({exc}); the rest- and planet-frame columns will be left out')
    if corr.flatcorr and not config.FLAT_FILE.is_file():
        try:
            log(f'Downloading the flat field ({FLAT_URL}) ...')
            config.FLAT_FILE.parent.mkdir(parents=True, exist_ok=True)
            download_flat(config.FLAT_FILE)
        except Exception as exc:
            log(f'  no flat field ({exc}); FLATW will be left out')
    return tag


# ---- writing the file -----------------------------------------------------------------------------------------------

STRUCTURE_KEYWORD = re.compile(r'T(TYPE|FORM|UNIT|NULL|SCAL|ZERO|DISP|BCOL|DIM)\d+')   # astropy writes these afresh
FITS_FORMATS = {'f8': 'D', 'f4': 'E', 'i2': 'I', 'i4': 'J', 'i8': 'K', 'u1': 'B', 'b1': 'L'}
# Per-column coordinate keywords (TCTYPn, TCUNIn, TCRPXn, TCRVLn, TCDLTn, TRPOSn). astropy keeps these as attributes
# of the column and writes the header from them, so a rebuilt column has to be given them.
COORDINATE_ATTRIBUTES = ('coord_type', 'coord_unit', 'coord_ref_point', 'coord_ref_value', 'coord_inc', 'time_ref_pos')


def stored_column(hdu, name):
    """One column of the TAG table as a new, unscaled column.

    Columns such as TIME are stored as scaled integers (TSCAL / TZERO). astropy hands out their real values, so they
    are written back as plain floating point: re-using the old column definition would make astropy try to scale the
    values back, and fail.
    """
    column = hdu.columns[name]
    values = np.asarray(hdu.data[name])
    code = FITS_FORMATS.get(f'{values.dtype.kind}{values.dtype.itemsize}')
    if code is None or values.ndim != 1:
        return column   # (not the case for TAG files) anything unusual keeps its own definition
    coordinates = {attribute: getattr(column, attribute) for attribute in COORDINATE_ATTRIBUTES}
    return fits.Column(name=name, format=code, unit=column.unit, array=values, **coordinates)


def in_file_order(values: np.ndarray, order) -> np.ndarray:
    """Columns are worked out in time order; put them back in the file's row order."""
    if order is None:
        return values
    restored = np.empty_like(values)
    restored[order] = values
    return restored


def column_definition(name: str, values: np.ndarray):
    """A new column: a flag (L), or a float32 position (pixel) or weight (E)."""
    if values.dtype == bool:
        return fits.Column(name=name, format='L', array=values)
    return fits.Column(name=name, format='E', array=values.astype(np.float32),
                       unit='pixel' if name[0] in 'XY' else None)


STEP_KEYWORDS = ('RANDCORR', 'GEOCORR', 'FLATCORR', 'VIGNCORR', 'EXPTCORR')


def primary_keywords(ds: TagDataset, time_sel: TimeSel | None, region: Region | None, frame: str, corr: Corrections,
                     columns: TagselColumns, script_version: str, tool_url: str) -> list:
    """(keyword, value, comment) for the primary header: which steps were performed, which selections were made, and
    the versions and numbers behind them. A step is PERFORM only when it took effect (the flat field must have been
    found, the vignetting model must apply to this aperture); otherwise OMIT. Needs no FITS library."""
    from . import __version__
    done = {'RANDCORR': bool(corr.randcorr), 'GEOCORR': bool(corr.geocorr), 'FLATCORR': 'FLATW' in columns.arrays,
            'VIGNCORR': 'VIGW' in columns.arrays, 'EXPTCORR': 'EDGE_OK' in columns.arrays}
    what = {'RANDCORR': 'random sub-pixel offset (XCORR, YCORR and the frame columns)',
            'GEOCORR': 'geometric distortion correction (XCORR, YCORR and the frame columns)',
            'FLATCORR': 'flat-field weights (FLATW)', 'VIGNCORR': 'vignetting model, F25SRF2 only (VIGW)',
            'EXPTCORR': 'detector-edge trim flag (EDGE_OK)'}
    keys = [(key, 'PERFORM' if done[key] else 'OMIT', f'{what[key]}: performed or omitted') for key in STEP_KEYWORDS]
    keys += [('TIMESEL', time_sel is not None, 'a time selection was made (SEL_TIME)'),
             ('SPCSEL', region is not None, 'a region was selected (SEL_REGION)')]
    if time_sel is not None:
        keys += [('SELTMIN', float(time_sel[0]), '[s] SEL_TIME window start'),
                 ('SELTMAX', float(time_sel[1]), '[s] SEL_TIME window end')]
    if region is not None:
        keys += [('SELREG', region.kind, 'SEL_REGION shape'), ('SELFRAME', frame, 'SEL_REGION reference frame')]
    if done['RANDCORR']:
        keys.append(('RANDSEED', RANDOM_SEED, 'seed of the random sub-pixel offsets'))
    if done['VIGNCORR']:
        keys += [('VIGBSTRT', float(corr.vigbot_start), '[pix] bottom ramp starts at this row (0-indexed)'),
                 ('VIGBEND', float(corr.vigbot_end), '[pix] bottom ramp: row where the loss is reached'),
                 ('VIGBLOSS', float(corr.vigbot_loss), 'bottom ramp: fractional loss at that row'),
                 ('VIGTSTRT', float(corr.vigtop_start), '[pix] top ramp starts at this row (0-indexed)'),
                 ('VIGTEND', float(corr.vigtop_end), '[pix] top ramp: row where the loss is reached'),
                 ('VIGTLOSS', float(corr.vigtop_loss), 'top ramp: fractional loss at that row')]
    if done['EXPTCORR']:
        keys += [(key, pixels, f'[pix] EDGE_OK trim of the {side} edge') for key, side, pixels in (
            ('SELEDGEL', 'left', EDGE_TRIM.left), ('SELEDGER', 'right', EDGE_TRIM.right),
            ('SELEDGEB', 'bottom', EDGE_TRIM.bottom), ('SELEDGET', 'top', EDGE_TRIM.top))]
    keys += [('SELVER', script_version, 'TIME-TAG Explorer selection script version'),
             ('TTEVER', __version__, 'timetag-explorer package version'),
             ('SELURL', tool_url, 'tool URL')]
    return keys


def _quiet_unit_warning(function):
    """Run `function` without astropy's warning that 'seconds' (the unit STIS writes for TIME) is not a FITS-standard
    unit name. The FITS unit parser accepts only the standard names, so a custom unit cannot be enabled for it; the
    warning is harmless (the values are read as written) and appears for every file."""
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UnitsWarning, message=r'.*did not parse as fits unit.*')
            return function(*args, **kwargs)
    return wrapper


@_quiet_unit_warning
def make_tagsel(rootname: str, *, time_sel: TimeSel | None = None, region: Region | None = None,
                frame: str = 'image', corr: Corrections, script_version: str, tool_url: str,
                north_up: bool = False, folder='.', overwrite: bool = False, log=print) -> Path:
    """Download what is missing, work out the columns and write `<rootname>_tagsel.fits` in `folder`. Returns its path.

    Refuses to replace an existing output unless `overwrite`. Needs a network connection for files that are not in
    `folder` yet (MAST, the flat field) and, for planet frames, for JPL Horizons.
    """
    from . import __version__
    if is_map(frame):
        raise ValueError('regions drawn on a planet map are not supported (maps are experimental)')
    if frame not in FRAMES:
        raise ValueError(f'unknown reference frame {frame!r}')
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    if not os.environ.get('TAG_EXPLORER_FLAT'):
        config.FLAT_FILE = folder / 'combined_flat.fits'     # (unless the reader keeps the flat elsewhere)
    if not os.environ.get('TAG_EXPLORER_HORIZONS'):
        config.HORIZONS_DIR = folder / 'horizons'
    root = rootname.lower()
    output = folder / f'{root}_tagsel.fits'
    if output.exists() and not overwrite:
        raise SystemExit(f'{output.name} already exists: set OVERWRITE = True to replace it')
    tag_path = fetch_inputs(root, folder, corr, log)

    log(f'Reading {tag_path.name} ...')
    ds = load_dataset(tag_path)
    columns = selection_columns(ds, time_sel, region, frame, corr, north_up)
    for what, why in columns.skipped:
        log(f'  left out: {what} ({why})')
    log(f'{len(ds.time):,} events')
    for name in columns.flags:
        log(f'  {name:<12}{int(columns.arrays[name].sum()):>12,}')
    if len(columns.flags) > 1:
        both = np.logical_and.reduce([columns.arrays[name] for name in columns.flags])
        log(f'  {"all of them":<12}{int(both.sum()):>12,}')

    history = [f'Columns added by the TIME-TAG Explorer selection script v{script_version} (timetag-explorer {__version__}), '
               f'run {datetime.now(timezone.utc):%Y-%m-%d} UTC',
               f'Tool: {tool_url}',
               'Not an STScI calibrated product: the raw _tag table plus the added columns',
               'Existing columns keep their values (AXIS1 and AXIS2 are the recorded integer pixels); '
               'scaled ones (TIME) are stored as plain floating point']
    history += columns.history
    history.append('Positions are hi-res pixels counted from 1, like AXIS1 and AXIS2; float32')
    history.append('Which steps were performed, the selections made and the versions are keywords of the primary header')

    with fits.open(tag_path) as hdul:
        raw_time = np.asarray(hdul[1].data['TIME'], dtype=np.float32)
        order = None if _is_sorted(raw_time) else np.argsort(raw_time, kind='stable')
        del raw_time
        table_columns = [stored_column(hdul[1], name) for name in hdul[1].columns.names]
        table_columns += [column_definition(name, in_file_order(values, order)) for name, values in columns.arrays.items()]
        header = hdul[1].header.copy()
        for key in [key for key in header if STRUCTURE_KEYWORD.fullmatch(key)]:
            del header[key]   # the old column definitions, TSCAL and TZERO included
        table = fits.BinTableHDU.from_columns(fits.ColDefs(table_columns), header=header)
        header = table.header
        for text in history:
            for line in textwrap.wrap(text, 68):
                header.add_history(line)
        primary = hdul[0].header
        replaced = []
        for key, value, comment in primary_keywords(ds, time_sel, region, frame, corr, columns, script_version, tool_url):
            if key in STEP_KEYWORDS and key in primary and primary[key] != value:
                replaced.append(f'{key} was {primary[key]!r}')
            primary[key] = (value, comment)
        primary.add_history('Columns added to extension 1; see its HISTORY. The processing steps and selections are '
                            'the keywords here (RANDCORR ... EXPTCORR, TIMESEL, SPCSEL, SELVER, TTEVER).')
        if replaced:   # the calibration pipeline's own flags of these names now say what this script did instead
            for line in textwrap.wrap('Values these step keywords had in the original file (the calibration '
                                      'pipeline\'s): ' + '; '.join(replaced), 68):
                primary.add_history(line)
        hdul[1] = table
        part = output.with_name(output.name + '.part')   # written aside, so a failure never leaves a broken output
        try:
            hdul.writeto(part, overwrite=True)
            os.replace(part, output)
        finally:
            part.unlink(missing_ok=True)
    log(f'Wrote {output}')
    return output
