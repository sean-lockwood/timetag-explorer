"""The processing switches, the reference frames, and the pointing information read from a TIME-TAG (and SPT) header."""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

from astropy.io import fits

from . import gdcorr
from .constants import DEFAULT_TIME_BIN, MAX_LIGHT_CURVE_BINS, MIN_TIME_BIN
from .log import logger
from .planet import BODIES



# --------------------------------------------------------------------------
# Corrections and reference frames
# --------------------------------------------------------------------------

def geocorr_available() -> bool:
    """False while gdcorr.py holds placeholder coefficients: the correction would change nothing, so it is left off."""
    return not gdcorr.PLACEHOLDER


@dataclass(frozen=True)
class Corrections:
    """Processing steps the user can switch on and off (hashable: used in cache keys)."""

    randcorr: bool = True   # +/-0.5 hi-res-pixel random jitter of each event position
    geocorr: bool = True    # geometric distortion correction (gdcorr.py, via a lookup table)
    flatcorr: bool = True   # flat-field weights
    timecorr: bool = True   # divide by the time each frame pixel was actually on the detector
    # Planet maps and the limb overlay: offset added to the planet centre guessed from the
    # header (hi-res pixels), and the emission altitude of the projection (None = the
    # planet's default). They live here so they are part of every cache key.
    centre_dx: float = 0.0
    centre_dy: float = 0.0
    altitude_km: float | None = None

    def map_key(self) -> tuple:
        return (self.centre_dx, self.centre_dy, self.altitude_km)


# geocorr is off (and its switch disabled) while gdcorr.py holds placeholder coefficients.
DEFAULT_CORRECTIONS = Corrections(geocorr=geocorr_available())

FRAMES = {
    'image': 'Image frame',
    'rest': 'Target rest frame',
    'rest_north': 'Target rest frame, north up',
    'planet': 'Planet rest frame',
    'planet_north': 'Planet rest frame, north up',
    'map_rect': 'Planet map (longitude / latitude)',
    'map_north': 'Planet map, north polar',
    'map_south': 'Planet map, south polar',
}
ARCSEC_PER_LOWRES_PIXEL = 0.024468  # 1 low-res pixel = 2 hi-res pixels
DETECTOR_CENTRE = 1023.5            # centre of the 2048-pixel detector, 0-indexed hi-res pixels

# JPL Horizons major-body ids for the names used in the SPT motion-tracking levels.
HORIZONS_IDS = {
    'JUPITER': '599', 'IO': '501', 'EUROPA': '502', 'GANYMEDE': '503', 'CALLISTO': '504',
    'SATURN': '699', 'MIMAS': '601', 'ENCELADUS': '602', 'TETHYS': '603', 'DIONE': '604',
    'RHEA': '605', 'TITAN': '606', 'IAPETUS': '608',
}


def is_north(frame: str) -> bool:
    return frame in ('rest_north', 'planet_north')


def is_planet(frame: str) -> bool:
    return frame in ('planet', 'planet_north')


def is_map(frame: str) -> bool:
    return frame.startswith('map_')


MAX_ALTITUDE_KM = 100_000.0     # planet-registration inputs: keep them physically sensible
MAX_CENTRE_SHIFT_PX = 2048.0


def _clip_finite(value, low: float, high: float) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return max(low, min(0.0, high))
    return min(max(value, low), high) if math.isfinite(value) else max(low, min(0.0, high))


def clamp_binsize(value, tmax: float | None = None) -> float:
    """A usable light-curve bin: finite, at least MIN_TIME_BIN, and (once the exposure length is
    known) not so small that the light curve would need more than MAX_LIGHT_CURVE_BINS bins.
    A bin of 1e-9 s typed into the box would otherwise ask for billions of bins."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        value = DEFAULT_TIME_BIN
    if not math.isfinite(value):
        value = DEFAULT_TIME_BIN
    floor = MIN_TIME_BIN if not tmax else max(MIN_TIME_BIN, tmax / MAX_LIGHT_CURVE_BINS)
    return max(value, floor)


@dataclass(frozen=True)
class FrameInfo:
    """What is needed to move events from the image frame to the rest frames."""

    orientat: float | None    # deg, ORIENTAT from the tag file's ext=1 header
    slew_angle: float | None  # deg, ANG of the last POS_ANGLE motion-tracking level
    slew_rate: float          # arcsec/day, R_RAD of that level (0 = no time-dependent slew)
    levels: tuple = ()        # parsed MT_LV levels (list of dicts), for reference
    problem: str | None = None  # why the rest frames can't be built (e.g. no SPT file)
    planet: str | None = None   # first STD level, when a second (moon) STD level follows
    target: str | None = None   # last STD level (what HST tracked)
    cd: tuple | None = None     # WCS CD matrix (1_1, 1_2, 2_1, 2_2), deg per hi-res pixel
    expstart: float | None = None  # MJD (UTC)
    expend: float | None = None
    wcs_ref: tuple | None = None   # (CRPIX1, CRPIX2, CRVAL1, CRVAL2) of the event-table WCS

    @property
    def body(self) -> str | None:
        """The planet the map frames project onto (Jupiter or Saturn), if this is such an observation."""
        for name in (self.planet, self.target):
            if name in BODIES:
                return name
        return None

    def map_base_frame(self) -> str:
        """The pixel frame in which the planet's centre stays fixed."""
        return 'planet' if self.planet else 'rest'

    def unavailable_reason(self, frame: str) -> str | None:
        if frame == 'image':
            return None
        if self.problem:
            return self.problem
        if is_map(frame):
            if self.body is None:
                return 'maps need a Jupiter or Saturn (or one of their moons) observation'
            reason = self.unavailable_reason(self.map_base_frame())
            if reason:
                return reason
            if self.cd is None:
                return 'tag file header has no usable WCS (TC2_2, TC2_3, TC3_2, TC3_3)'
            if self.expstart is None or self.expend is None:
                return 'EXPSTART/EXPEND missing from the tag file header'
            return None
        moving = self.slew_rate != 0
        if moving and self.slew_angle is None:
            return 'slew level has R_RAD but no ANG'
        if (moving or is_north(frame)) and self.orientat is None:
            return 'ORIENTAT is missing from the tag file header'
        if is_planet(frame):
            if not self.planet or not self.target:
                return 'not a moon observation (needs STD=<planet> then STD=<moon> levels)'
            if self.planet not in HORIZONS_IDS or self.target not in HORIZONS_IDS:
                return f'no ephemeris id for {self.planet}/{self.target}'
            if self.cd is None:
                return 'tag file header has no usable WCS (TC2_2, TC2_3, TC3_2, TC3_3)'
            if self.expstart is None or self.expend is None:
                return 'EXPSTART/EXPEND missing from the tag file header'
        return None


def frame_label(info: FrameInfo, frame: str) -> str:
    target = (info.target or 'target').title()
    planet = (info.planet or 'planet').title()
    body = (info.body or 'planet').title()
    return {
        'image': 'Image frame',
        'rest': f'{target} rest frame',
        'rest_north': f'{target} rest frame, north up',
        'planet': f'{planet} rest frame',
        'planet_north': f'{planet} rest frame, north up',
        'map_rect': f'{body} map (longitude / latitude)',
        'map_north': f'{body} map, north polar',
        'map_south': f'{body} map, south polar',
    }[frame]


def parse_mt_levels(header) -> list[dict[str, str]]:
    """Motion-tracking levels from MT_LV<level>_<part> keywords.

    A trailing part number > 1 continues the previous part of the same level
    (FITS values are length-limited), so the parts are concatenated in order.
    Each level becomes a dict such as {'TYPE': 'POS_ANGLE', 'ANG': '44.0', ...}.
    """
    parts: dict[int, dict[int, str]] = {}
    for key in header:
        m = re.fullmatch(r'MT_LV(\d+)_(\d+)', str(key))
        if m:
            parts.setdefault(int(m[1]), {})[int(m[2])] = str(header[key])
    levels = []
    for level in sorted(parts):
        text = ''.join(parts[level][k].strip() for k in sorted(parts[level]))
        fields = {}
        for item in text.split(','):
            if '=' in item:
                name, value = item.split('=', 1)
                fields[name.strip().upper()] = value.strip()
        levels.append(fields)
    return levels


def slew_from_levels(levels: list[dict[str, str]]) -> tuple[float | None, float]:
    """(angle [deg], rate [arcsec/day]) of the last level if it is a POS_ANGLE slew.

    Backing out that last level takes the data from the image frame to the
    target's rest frame. Any other last level (e.g. plain STD=<target>
    tracking) means there is no extra slew to remove: (None, 0.0).
    """
    if not levels or levels[-1].get('TYPE', '').upper() != 'POS_ANGLE':
        return None, 0.0
    last = levels[-1]
    angle = float(last['ANG']) if 'ANG' in last else None
    rate = float(last['R_RAD']) if 'R_RAD' in last else 0.0
    if float(last.get('R_ANG', 0.0)) != 0.0:
        logger.warning('R_ANG=%s is not 0 and is ignored (rate assumed along ANG)', last['R_ANG'])
    return angle, rate


def names_from_levels(levels: list[dict[str, str]]) -> tuple[str | None, str | None]:
    """(planet, target) from the STD levels: STD=JUPITER, STD=EUROPA -> ('JUPITER', 'EUROPA')."""
    stds = [lv['STD'].upper() for lv in levels if 'STD' in lv and 'TYPE' not in lv]
    target = stds[-1] if stds else None
    planet = stds[0] if len(stds) >= 2 and stds[0] != stds[-1] else None
    return planet, target


def cd_from_header(header) -> tuple | None:
    """WCS CD matrix of the event table (the TC2_2, TC2_3, TC3_2, TC3_3 keywords),
    treated as degrees per hi-res pixel as parse_ttag_wcs does. None if absent or implausible."""
    try:
        cd = tuple(float(header[k]) for k in ('TC2_2', 'TC2_3', 'TC3_2', 'TC3_3'))
    except (KeyError, TypeError, ValueError):
        return None
    scale = math.sqrt(abs(cd[0] * cd[3] - cd[1] * cd[2])) * 3600.0  # arcsec per hi-res pixel
    if not 0.002 < scale < 0.1:
        logger.warning('WCS scale %.4g arcsec/pixel is implausible; ignoring the header WCS', scale)
        return None
    return cd


def wcs_ref_from_header(header) -> tuple | None:
    """(CRPIX1, CRPIX2, CRVAL1, CRVAL2) of the event table's WCS (TCRPX2/3, TCRVL2/3), or None."""
    try:
        return tuple(float(header[k]) for k in ('TCRPX2', 'TCRPX3', 'TCRVL2', 'TCRVL3'))
    except (KeyError, TypeError, ValueError):
        return None


def load_frame_info(tag_path: Path, tag_header) -> FrameInfo:
    """Slew + WCS description: <rootname>_spt.fits(.gz) next to the tag file, plus the tag's ext=1 header."""
    orientat = tag_header.get('ORIENTAT')
    base = dict(
        orientat=None if orientat is None else float(orientat),
        cd=cd_from_header(tag_header), wcs_ref=wcs_ref_from_header(tag_header),
        expstart=tag_header.get('EXPSTART'), expend=tag_header.get('EXPEND'),
    )
    rootname = tag_path.name.split('_', 1)[0]
    candidates = [tag_path.with_name(f'{rootname}_spt.fits'),
                  tag_path.with_name(f'{rootname}_spt.fits.gz')]
    spt = next((c for c in candidates if c.is_file()), None)
    if spt is None:
        return FrameInfo(**base, slew_angle=None, slew_rate=0.0,
                         problem=f'{rootname}_spt.fits not found next to the tag file')
    try:
        levels = parse_mt_levels(fits.getheader(spt, ext=0))
        if not levels:
            return FrameInfo(**base, slew_angle=None, slew_rate=0.0,
                             problem=f'no MT_LV keywords in {spt.name}')
        angle, rate = slew_from_levels(levels)
    except (OSError, ValueError) as exc:
        return FrameInfo(**base, slew_angle=None, slew_rate=0.0,
                         problem=f'could not read slew from {spt.name}: {exc}')
    planet, target = names_from_levels(levels)
    return FrameInfo(**base, slew_angle=angle, slew_rate=rate, levels=tuple(levels),
                     planet=planet, target=target)
