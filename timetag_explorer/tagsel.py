"""The work a selection script does: flag the selected events, add the corrected positions and weights, and write
`<rootname>_tagsel.fits`.

A script written by `build_selection_script` is only the choices (named time ranges and regions, corrections) and a call
to `make_tagsel`, so the numbers come from the same code as the Explorer's own images and light curves.

The output is the exposure's TIME-TAG table with columns added; the columns that were there (AXIS1 and AXIS2 included)
are not changed. Added columns, in time-tag row order:

* flags: `SEL_TIME` and `SEL_REGION` (`I` with TZERO 32768, unsigned 16-bit): one bit per time range or region, set
  where the event is inside it, as the `TIMESEL` and `REGSEL` tables say (SELECTIONS.md); `EDGE_OK` (`L`);
* positions (`E`, float32, hi-res pixels counted from 1 like AXIS1 and AXIS2): `XCORR`, `YCORR` when the random
  sub-pixel offset or the geometric distortion correction is on or a region was drawn in the image frame (the
  image-frame position with both applied), and one pair per reference frame the exposure supports, from those corrected
  positions: `XREST`/`YREST` (target rest frame) and `XPLAN`/`YPLAN` (planet rest frame, for moon observations). The
  north-up frames, `XRESTN`/`YRESTN` and `XPLANN`/`YPLANN`, are written only with `north_up=True` or when a region was
  drawn in one (they make the file larger);
* planet-map positions (`E`, float32, degrees; experimental, Jupiter and Saturn observations): `MAPLON`/`MAPLAT`, the
  west longitude (unwrapped around the map's centre longitude, so it can lie outside 0-360) and planetographic latitude
  of the point of the planet's emission surface each event is seen at, NaN for events that miss the planet. The polar
  maps, `MAPNX`/`MAPNY` (north) and `MAPSX`/`MAPSY` (south), the azimuthal-equidistant planes the Explorer's polar maps
  are drawn in (colatitude in degrees from the pole, longitude 180 at the top), only when a region was drawn in one.
  They depend on the planet registration (the centre offsets and the emission altitude), which the primary header
  records (`MAPDX`, `MAPDY`, `MAPALT`, `MAPLCEN`);
* weights (`E`): `FLATW` (flat field) and `VIGW` (vignetting model, F25SRF2 only), for the corrections that are on.
  An event's weight is their product; they are separate columns so either can be left out of a calculation.

Regions are tested on exactly the values written to the position columns, so applying a region from `REGSEL` to those
columns by hand gives the same flags.
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

from . import config, gdcorr
from .constants import EDGE_TRIM, FLAT_URL, RANDOM_SEED, SELECTION_SCRIPT_VERSION, TIME_DECIMALS

PACKAGE_URL = 'https://github.com/sean-lockwood/timetag-explorer'    # the SELURL of a file made by calling the package directly
from .dataset import TagDataset, _is_sorted, load_dataset
from .frames import FRAMES, Corrections, frame_label, is_map, is_north, is_planet
from .regions import region_source
from .selections import RegionTable, Selections, TimeTable
from .transforms import frame_coordinates, load_geo_lut, map_coordinates, map_longitude_centre

# Columns of the positions in each reference frame (the image frame is XCORR / YCORR).
POSITION_COLUMNS = {'rest': ('XREST', 'YREST'), 'rest_north': ('XRESTN', 'YRESTN'),
                    'planet': ('XPLAN', 'YPLAN'), 'planet_north': ('XPLANN', 'YPLANN')}
CORRECTED_COLUMNS = ('XCORR', 'YCORR')
MAP_COLUMNS = {'map_rect': ('MAPLON', 'MAPLAT'), 'map_north': ('MAPNX', 'MAPNY'), 'map_south': ('MAPSX', 'MAPSY')}
MAP_KEYWORDS = ('MAPDX', 'MAPDY', 'MAPALT', 'MAPLCEN')   # the planet registration the map columns were made with
PIXEL_FRAMES = ('image',) + tuple(POSITION_COLUMNS)    # the frames with pixel position columns
EDGE_TRIM_TEXT = f'left {EDGE_TRIM.left}, right {EDGE_TRIM.right}, bottom {EDGE_TRIM.bottom}, top {EDGE_TRIM.top}'


@dataclass
class TagselColumns:
    """The columns to add, in time order (the dataset's order), and what was left out and why."""
    arrays: dict = field(default_factory=dict)     # name -> array, in the order they are written
    flags: tuple = ()                              # names of the flag columns
    skipped: list = field(default_factory=list)    # (what, reason)
    history: list = field(default_factory=list)    # lines for the FITS HISTORY
    times: TimeTable = field(default_factory=TimeTable)         # the time ranges, as flagged in SEL_TIME
    regions: RegionTable = field(default_factory=RegionTable)   # the regions, as flagged in SEL_REGION
    map_registration: dict = field(default_factory=dict)        # MAPDX, MAPDY, MAPALT, MAPLCEN, when map columns are written


def _positions(ds: TagDataset, frame: str, corr: Corrections):
    """(x, y) of every event in a frame as the Explorer computes them: 0-indexed hi-res pixels, float32."""
    track = ds.planet_track() if is_planet(frame) else None
    geo_lut = load_geo_lut(ds.opt_elem) if corr.geocorr else None
    return frame_coordinates(ds.x, ds.y, ds.time, ds.frame_info, frame, corr,
                             seed=RANDOM_SEED, geo_lut=geo_lut, track=track)


def _map_positions(ds: TagDataset, frame: str, corr: Corrections):
    """(u, v) of every event in a planet-map plane (degrees, float32; NaN for events that miss the planet), as the
    Explorer's map frames compute them."""
    reason = ds.frame_info.unavailable_reason(frame)
    if reason:
        raise ValueError(f'{frame_label(ds.frame_info, frame)} unavailable: {reason}')
    u, v, keep = map_coordinates(ds, frame, corr)
    u, v = np.where(keep, u, np.float32('nan')).astype(np.float32), np.where(keep, v, np.float32('nan')).astype(np.float32)
    return u, v


def map_registration(ds: TagDataset, corr: Corrections) -> dict:
    """What the planet-map columns depend on besides the pixels: the planet-centre offsets (hi-res pixels) and the emission
    altitude (km) the exposure was registered with, and the centre longitude of the rectilinear map."""
    geom = ds.planet_geometry()
    return {'MAPDX': float(corr.centre_dx), 'MAPDY': float(corr.centre_dy), 'MAPALT': float(geom.altitude_km(corr)),
            'MAPLCEN': float(map_longitude_centre(geom.orient))}


def map_history(ds: TagDataset, corr: Corrections) -> list:
    """HISTORY lines saying what the planet-map columns are made with."""
    r = map_registration(ds, corr)
    return [f'Planet-map columns (experimental): planet-centre offset ({r["MAPDX"]:g}, {r["MAPDY"]:g}) hi-res pixels and '
            f'emission altitude {r["MAPALT"]:g} km, as MAPDX, MAPDY and MAPALT; the longitude MAPLON is unwrapped around '
            f'{r["MAPLCEN"]:.2f} deg (MAPLCEN). A different registration gives different columns; the default is only a '
            'starting point and must be checked against the data.']


def _one_based(a: np.ndarray) -> np.ndarray:
    """0-indexed position to the 1-indexed convention of AXIS1 and AXIS2."""
    return a + np.float32(1.0)


def selection_columns(ds: TagDataset, times: TimeTable | None, regions: RegionTable | None, corr: Corrections,
                      north_up: bool = False) -> TagselColumns:
    """Work out the columns a selection script adds (time order). Needs ephemerides (the network) only for planet
    frames. A frame that cannot be built is skipped with its reason, unless a region was drawn in it: then this raises
    (RuntimeError), so nothing is written. Every frame a region was drawn in gets its position columns, north-up or
    not, so the region can be applied by hand; the other north-up frames only when `north_up`."""
    times = times if times is not None else TimeTable()
    regions = regions if regions is not None else RegionTable()
    out = TagselColumns(times=times, regions=regions)
    flags, positions = {}, {}
    n = len(ds.time)
    if times:
        bits = np.zeros(n, dtype=np.uint16)
        for row in times.primitives:
            bits[ds.time_slice((row.start, row.stop))] |= np.uint16(row.mask)
            out.history.append(f'SEL_TIME bit {row.bit} ({row.name}): set for {row.start:.{TIME_DECIMALS}f} <= TIME < '
                               f'{row.stop:.{TIME_DECIMALS}f} s (seconds since the exposure start)')
        flags['SEL_TIME'] = bits

    # Positions in every frame that can be built, one frame at a time. The frames regions were drawn in come first:
    # a failure there is an error, not a gap.
    needed = regions.used_frames
    region_bits = np.zeros(n, dtype=np.uint16) if regions else None
    corrected = bool(corr.randcorr or corr.geocorr)
    todo = list(dict.fromkeys(list(needed) + ['image'] + list(POSITION_COLUMNS) + ['map_rect']))
    for f in todo:
        try:
            if f != 'image':
                ds.prepare_frame(f)
            u, v = _map_positions(ds, f, corr) if is_map(f) else _positions(ds, f, corr)
        except Exception as exc:   # an unavailable frame (ValueError), or Horizons not reachable
            if f in needed:
                raise RuntimeError(f'a region was drawn in the {frame_label(ds.frame_info, f)}, which cannot be '
                                   f'built: {exc}') from exc
            if is_map(f) and ds.frame_info.unavailable_reason(f):
                continue      # not a Jupiter or Saturn exposure: nothing was expected, so nothing is reported as left out
            if north_up or not is_north(f):
                out.skipped.append((frame_label(ds.frame_info, f), str(exc)))
            continue
        if is_map(f):    # degrees as they are (a region on a map is not shifted by one), NaN: off the planet
            x1, y1 = u, v
        else:
            x1, y1 = _one_based(u), _one_based(v)    # what the position columns hold; the regions are tested on these
        del u, v
        for row in regions.primitives:
            if row.frame == f:
                inside = np.asarray(row.region.contains(x1, y1), dtype=bool)
                if is_map(f):
                    inside &= np.isfinite(x1) & np.isfinite(y1)      # (events that miss the planet are in no region)
                region_bits[inside] |= np.uint16(row.mask)
                out.history.append(f'SEL_REGION bit {row.bit} ({row.name}): set for events inside this region of the '
                                   f'{frame_label(ds.frame_info, f)}, tested at each event\'s own position as written '
                                   f'in the position columns: {row.region.describe()}')
                out.history.append(f'Region {row.name}, unrounded: ' + region_source(row.region))
        if f == 'image':
            if corrected or f in needed:
                positions['XCORR'], positions['YCORR'] = x1, y1
                what = ' and '.join(name for name, on in (('random sub-pixel offset', corr.randcorr),
                                                         ('geometric distortion', corr.geocorr)) if on)
                out.history.append('XCORR, YCORR: the image-frame position' +
                                   (f' with {what} applied' if what else ' (no correction applied)'))
        elif is_map(f):
            if f == 'map_rect' or f in needed:
                positions[MAP_COLUMNS[f][0]], positions[MAP_COLUMNS[f][1]] = x1, y1
                where = {'map_rect': 'west longitude (unwrapped around the map centre) and planetographic latitude',
                         'map_north': 'x and y in the north polar map plane (degrees of colatitude)',
                         'map_south': 'x and y in the south polar map plane (degrees of colatitude)'}[f]
                out.history.append(f'{MAP_COLUMNS[f][0]}, {MAP_COLUMNS[f][1]}: {where}, in degrees, on the emission surface '
                                   '(experimental); NaN where the event misses the planet')
        elif north_up or not is_north(f) or f in needed:
            xname, yname = POSITION_COLUMNS[f]
            positions[xname], positions[yname] = x1, y1
            out.history.append(f'{xname}, {yname}: the position in the {frame_label(ds.frame_info, f)}, from the '
                               'corrected positions')
        del x1, y1
    if regions:
        flags['SEL_REGION'] = region_bits
    # Always the same order (XCORR, YCORR, then the frames as listed in POSITION_COLUMNS), whichever frame was built first.
    order = (['XCORR', 'YCORR'] + [name for pair in POSITION_COLUMNS.values() for name in pair]
             + [name for pair in MAP_COLUMNS.values() for name in pair])
    positions = {name: positions[name] for name in order if name in positions}
    if any(name in positions for pair in MAP_COLUMNS.values() for name in pair):
        out.history += map_history(ds, corr)
        out.map_registration = map_registration(ds, corr)
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

# STIS tables write the units 'seconds' and 'pixels', which the FITS standard does not list; astropy warns once per column
UNIT_WARNING_PATTERN = r'.*did not parse as fits unit.*'


def _quiet_unit_warning(function):
    """Run `function` without astropy's warning that 'seconds' (the unit STIS writes for TIME) is not a FITS-standard
    unit name. The FITS unit parser accepts only the standard names, so a custom unit cannot be enabled for it; the
    warning is harmless (the values are read as written) and appears for every file."""
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UnitsWarning, message=UNIT_WARNING_PATTERN)
            return function(*args, **kwargs)
    return wrapper


STRUCTURE_KEYWORD = re.compile(r'T(TYPE|FORM|UNIT|NULL|SCAL|ZERO|DISP|BCOL|DIM)\d+')   # astropy writes these afresh
FITS_FORMATS = {'f8': 'D', 'f4': 'E', 'i2': 'I', 'i4': 'J', 'i8': 'K', 'u1': 'B', 'b1': 'L'}
# Per-column coordinate keywords (TCTYPn, TCUNIn, TCRPXn, TCRVLn, TCDLTn, TRPOSn). astropy keeps these as attributes
# of the column and writes the header from them, so a rebuilt column has to be given them.
COORDINATE_ATTRIBUTES = ('coord_type', 'coord_unit', 'coord_ref_point', 'coord_ref_value', 'coord_inc', 'time_ref_pos')


def stored_column(hdu, name):
    """One column of the TAG table as a new, unscaled column.

    A column stored with a scale (TSCAL / TZERO) is handed out by astropy with its real values, so it is written
    back as plain floating point: re-using the old column definition would make astropy try to scale the
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
    """A new column: a bit-flag word (unsigned 16-bit, stored as `I` with TZERO 32768), a flag (L), or a float32
    position (pixel) or weight (E)."""
    if values.dtype == np.uint16:
        return fits.Column(name=name, format='I', bzero=32768, array=values)
    if values.dtype == bool:
        return fits.Column(name=name, format='L', array=values)
    unit = 'pixel' if name[0] in 'XY' else 'deg' if name.startswith('MAP') else None
    return fits.Column(name=name, format='E', array=values.astype(np.float32), unit=unit)


TABLE_UNITS = {'TSTART': 's', 'TSTOP': 's', 'PARAMS': 'pixel'}
TABLE_FORMATS = {'BIT': 'I', 'MASK': 'J', 'TSTART': 'D', 'TSTOP': 'D', 'VISIBLE': 'L', 'CURVE': 'L'}


def selection_table_hdu(table, extname: str):
    """The `TIMESEL` or `REGSEL` extension for a table of selections (see SELECTIONS.md)."""
    columns = []
    for name, values in table.to_columns().items():
        if name == 'PARAMS':
            columns.append(fits.Column(name=name, format=f'{values.shape[1]}D', unit=TABLE_UNITS[name], array=values))
        elif name in TABLE_FORMATS:
            columns.append(fits.Column(name=name, format=TABLE_FORMATS[name], unit=TABLE_UNITS.get(name),
                                       array=np.asarray(values)))
        else:   # text
            width = max([1] + [len(text) for text in values])
            columns.append(fits.Column(name=name, format=f'{width}A', array=np.array(values, dtype=f'U{width}')))
    hdu = fits.BinTableHDU.from_columns(fits.ColDefs(columns), name=extname)
    hdu.header['EXTVER'] = (1, 'extension version')
    hdu.header.add_comment('The selections that SEL_TIME / SEL_REGION flag, one bit each; see SELECTIONS.md in')
    hdu.header.add_comment('https://github.com/sean-lockwood/timetag-explorer.')
    hdu.header.add_comment('Times are seconds since the exposure start.' if extname == 'TIMESEL'
                           else 'Positions are hi-res pixels counted from 1, like AXIS1.')
    return hdu


@_quiet_unit_warning
def read_selections(path) -> Selections:
    """The time ranges and regions recorded in a selection file (`TIMESEL` and `REGSEL`); an extension that is not there
    gives an empty table. Use `.in_time(t, name_or_bit)` and `.in_region(t, name_or_bit)` on the events table."""
    def columns(hdu):
        return {name: np.array(hdu.data[name]) for name in hdu.columns.names}
    with fits.open(path) as hdul:
        names = [hdu.name for hdu in hdul]
        times = TimeTable.from_columns(columns(hdul['TIMESEL'])) if 'TIMESEL' in names else TimeTable()
        regions = RegionTable.from_columns(columns(hdul['REGSEL'])) if 'REGSEL' in names else RegionTable()
    return Selections(times, regions)


# The primary-header keywords are this tool's own (prefix TTE), never the calibration pipeline's step keywords
# (RANDCORR, FLATCORR ... mean PERFORM/OMIT/COMPLETE to calstis and stistools), so later processing is unaffected.
STEP_KEYWORDS = ('TTERAND', 'TTEGEO', 'TTEFLAT', 'TTEVIGN', 'TTEEDGE', 'TTEMAP')   # DONE or OMIT
KEYWORD_GROUP = 'TIMETAG-EXPLORER KEYWORDS'                                # the label above the block, `/ TIMETAG-EXPLORER KEYWORDS`


def primary_keywords(times: TimeTable, regions: RegionTable, corr: Corrections, columns: TagselColumns,
                     script_version: str, tool_url: str) -> list:
    """(keyword, value, comment) for the primary header, in the order they are written as one block: which steps were
    done (DONE or OMIT), the flat file used, whether time ranges and regions were selected (the TIMESEL and REGSEL
    extensions hold them), and the numbers and versions behind them. A step is DONE only when it took effect (the flat
    field must have been found, the vignetting model must apply to this aperture); otherwise OMIT. Needs no FITS
    library."""
    from . import __version__
    done = {'TTERAND': bool(corr.randcorr), 'TTEGEO': bool(corr.geocorr), 'TTEFLAT': 'FLATW' in columns.arrays,
            'TTEVIGN': 'VIGW' in columns.arrays, 'TTEEDGE': 'EDGE_OK' in columns.arrays,
            'TTEMAP': bool(columns.map_registration)}
    what = {'TTERAND': 'random sub-pixel offset', 'TTEGEO': 'geometric distortion',
            'TTEFLAT': 'flat-field weights FLATW', 'TTEVIGN': 'vignetting model VIGW',
            'TTEEDGE': 'edge-trim flag EDGE_OK', 'TTEMAP': 'planet-map columns MAP... (experimental)'}   # (short: a FITS card holds 80 characters, value included)
    keys = [(key, 'DONE' if done[key] else 'OMIT', what[key]) for key in STEP_KEYWORDS]
    if done['TTEFLAT']:
        keys.append(('TTEFFILE', config.FLAT_FILE.name, 'flat-field file used for FLATW'))
    if done['TTEGEO']:
        keys.append(('TTEGFILE', gdcorr.idc_file().name, 'IDC table used for the distortion'))
    keys += [('TTETSEL', bool(times), 'time ranges selected: SEL_TIME, TIMESEL'),
             ('TTERSEL', bool(regions), 'regions selected: SEL_REGION, REGSEL')]
    if done['TTERAND']:
        keys.append(('RANDSEED', RANDOM_SEED, 'seed of the random sub-pixel offsets'))
    if done['TTEVIGN']:
        keys += [('VIGBSTRT', float(corr.vigbot_start), '[pix] bottom ramp starts at this row (0-indexed)'),
                 ('VIGBEND', float(corr.vigbot_end), '[pix] bottom ramp: row where the loss is reached'),
                 ('VIGBLOSS', float(corr.vigbot_loss), 'bottom ramp: fractional loss at that row'),
                 ('VIGTSTRT', float(corr.vigtop_start), '[pix] top ramp starts at this row (0-indexed)'),
                 ('VIGTEND', float(corr.vigtop_end), '[pix] top ramp: row where the loss is reached'),
                 ('VIGTLOSS', float(corr.vigtop_loss), 'top ramp: fractional loss at that row')]
    if done['TTEMAP']:
        r = columns.map_registration
        keys += [('MAPDX', r['MAPDX'], '[pix] planet-centre offset in x, hi-res pixels'),
                 ('MAPDY', r['MAPDY'], '[pix] planet-centre offset in y, hi-res pixels'),
                 ('MAPALT', r['MAPALT'], '[km] emission altitude of the map projection'),
                 ('MAPLCEN', r['MAPLCEN'], '[deg] centre longitude of MAPLON (west)')]
    if done['TTEEDGE']:
        keys += [(key, pixels, f'[pix] EDGE_OK trim of the {side} edge') for key, side, pixels in (
            ('SELEDGEL', 'left', EDGE_TRIM.left), ('SELEDGER', 'right', EDGE_TRIM.right),
            ('SELEDGEB', 'bottom', EDGE_TRIM.bottom), ('SELEDGET', 'top', EDGE_TRIM.top))]
    keys += [('SELVER', script_version, 'TIME-TAG Explorer selection script version'),
             ('TTEVER', __version__, 'timetag-explorer package version'),
             ('SELURL', tool_url, 'tool URL')]
    return keys


def _as_table(table, table_class):
    return table if isinstance(table, table_class) else table_class(table or ())


@_quiet_unit_warning
def make_tagsel(rootname: str, *, times=None, regions=None, corr: Corrections, script_version: str = SELECTION_SCRIPT_VERSION,
                tool_url: str = PACKAGE_URL, north_up: bool = False, folder='.', overwrite: bool = False, log=print) -> Path:
    """Download what is missing, work out the columns and write `<rootname>_tagsel.fits` in `folder`. Returns its path.

    Called directly (not from a generated script), only `rootname`, `corr` and the selections are needed: `script_version`
    and `tool_url` default to this package's own, and every selection may leave out its name, bit, color and visibility.

    `times` and `regions` are the named selections (a `TimeTable` / `RegionTable`, or a list of `TimeRange` /
    `RegionSelection`), the regions in file coordinates (hi-res pixels counted from 1). Anything wrong with them, such as
    a region on a planet map with no planet in the exposure, stops the script here, before anything is downloaded or
    written; so does a region whose frame cannot be built. Regions on planet maps (experimental) are written in degrees
    and tested on the `MAP...` columns, which depend on `corr`'s planet registration (`centre_dx`, `centre_dy`,
    `altitude_km`). Refuses to replace an existing output unless `overwrite`. Needs a network
    connection for files that are not in `folder` yet (MAST, the flat field) and, for planet frames, for JPL Horizons.
    """
    from . import __version__
    times, regions = _as_table(times, TimeTable), _as_table(regions, RegionTable)   # (checks them)
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
    columns = selection_columns(ds, times, regions, corr, north_up)
    for what, why in columns.skipped:
        log(f'  left out: {what} ({why})')
    log(f'{len(ds.time):,} events')
    for row in times.primitives:
        log(f'  time    {row.name:<12}{int(np.count_nonzero(columns.arrays["SEL_TIME"] & np.uint16(row.mask))):>12,}')
    for row in regions.primitives:
        log(f'  region  {row.name:<12}{int(np.count_nonzero(columns.arrays["SEL_REGION"] & np.uint16(row.mask))):>12,}')
    if 'EDGE_OK' in columns.arrays:
        log(f'  EDGE_OK {"":<12}{int(columns.arrays["EDGE_OK"].sum()):>12,}')

    history = [f'Columns added by the TIME-TAG Explorer selection script v{script_version} (timetag-explorer {__version__}), '
               f'run {datetime.now(timezone.utc):%Y-%m-%d} UTC',
               f'Tool: {tool_url}',
               'Not an STScI calibrated product: the raw _tag table plus the added columns',
               'Existing columns keep their values (AXIS1 and AXIS2 are the recorded integer pixels)']
    history += columns.history
    history.append('Positions are hi-res pixels counted from 1, like AXIS1 and AXIS2; float32')
    history.append('Which steps were performed and the versions are keywords of the primary header; the time ranges and '
                   'regions are the TIMESEL and REGSEL extensions')

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
        # One block at the end of the header, under a label card. Every card is appended with end=True, at the very end:
        # astropy puts an ordinary keyword before any trailing blank cards, and the label is a blank-keyword card, so
        # without it each keyword would land above the label. The calibration pipeline's own keywords are left exactly
        # as they were.
        keywords = primary_keywords(times, regions, corr, columns, script_version, tool_url)
        for key, _, _ in keywords:
            if key in primary:
                del primary[key]
        primary.append(('', f'/ {KEYWORD_GROUP}'), end=True)
        for card in keywords:
            primary.append(card, end=True)
        primary.add_history('Columns added to extension 1; see its HISTORY. The processing steps, the flat file and '
                            f'the selections are the {KEYWORD_GROUP.lower()} (TTE...) of this header.')
        hdul[1] = table
        if times:
            hdul.append(selection_table_hdu(times, 'TIMESEL'))
        if regions:
            hdul.append(selection_table_hdu(regions, 'REGSEL'))
        part = output.with_name(output.name + '.part')   # written aside, so a failure never leaves a broken output
        try:
            hdul.writeto(part, overwrite=True)
            os.replace(part, output)
        finally:
            part.unlink(missing_ok=True)
    log(f'Wrote {output}')
    return output
