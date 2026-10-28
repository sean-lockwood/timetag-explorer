"""Coordinate transformations: detector pixels to reference frames and planet maps, the exposure-time maps that go
with them, and the geometric-distortion lookup table."""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING

import numpy as np

from . import gdcorr
from .constants import (DETECTOR_LIMITS, EXPOSURE_SUBSAMPLE, MAP_BIN_DEGREES, MAP_SEGMENTS, MAP_SLICE_SECONDS,
                        PLANET_SEGMENTS, RANDOM_SEED, SPATIAL_BIN, _chunks)
from .frames import (ARCSEC_PER_LOWRES_PIXEL, DETECTOR_CENTRE, FRAMES, Corrections, FrameInfo, frame_label, is_map,
                     is_north, is_planet)
from .horizons import PlanetTrack
from .log import logger
from .planet import (BODIES, MAP_MIN_MU, Body, PlanetOrientation, limb_outline, map_bin_area, map_plane_coordinates,
                     map_plane_to_surface, pixels_to_screen, planet_orientation, screen_to_pixels, screen_to_surface,
                     spheroid_radii, surface_to_screen)

if TYPE_CHECKING:   # only in annotations: dataset.py imports this module
    from .dataset import TagDataset



@dataclass(frozen=True)
class PlanetGeometry:
    """Everything needed to project between the detector and a planet's surface."""

    body: Body
    orient: PlanetOrientation
    centre0: tuple            # planet centre guessed from the header, in planet-frame pixels
    cd: tuple
    cd_inv: np.ndarray
    pixel_area_arcsec2: float

    def centre(self, corr: Corrections) -> tuple[float, float]:
        return self.centre0[0] + corr.centre_dx, self.centre0[1] + corr.centre_dy

    def altitude_km(self, corr: Corrections) -> float:
        return self.body.altitude_km if corr.altitude_km is None else corr.altitude_km


def _interp_angle(t, degrees, tq):
    return float(np.degrees(np.interp(tq, t, np.unwrap(np.radians(degrees)))) % 360.0)


def build_planet_geometry(info: FrameInfo, positions: np.ndarray) -> PlanetGeometry:
    """Orientation from the IAU model, and a first guess of the planet's centre.

    The guess is the planet's ephemeris position at the exposure start relative to the
    event table's WCS reference point (TCRVL / TCRPX in the header). It is only as good
    as the pointing, so the user can add an offset.
    """
    body = BODIES[info.body]
    orient = planet_orientation(positions, body, info.expstart)
    cd = np.array(info.cd, dtype=float).reshape(2, 2)
    cd_inv = np.linalg.inv(cd)
    if info.wcs_ref is not None:
        crpix1, crpix2, crval1, crval2 = info.wcs_ref
        d_xi = ((_interp_angle(orient.t, orient.ra, 0.0) - crval1 + 180.0) % 360.0 - 180.0) \
            * math.cos(math.radians(crval2))
        d_eta = float(np.interp(0.0, orient.t, orient.dec)) - crval2
        off = cd_inv @ np.array([d_xi, d_eta])
        centre0 = (crpix1 - 1.0 + off[0], crpix2 - 1.0 + off[1])
    else:
        logger.warning('no WCS reference in the header; starting the planet centre at the detector centre')
        centre0 = (DETECTOR_CENTRE, DETECTOR_CENTRE)
    return PlanetGeometry(body, orient, centre0, tuple(cd.ravel()), cd_inv,
                          abs(float(np.linalg.det(cd))) * 3600.0 ** 2)


# ---- geometric distortion correction (gdcorr.py) as a lookup table --------------------------

@lru_cache(maxsize=1)
def load_geo_lut():
    """Distortion offsets (dx, dy) in hi-res pixels at every integer detector pixel.

    gdcorr is a smooth polynomial, so a per-pixel table gives the same answer as
    evaluating it per event (the offset changes by < 0.01 px across a pixel) at the
    cost of one gather per event. Built once per process (about half a second).
    """
    axis = np.arange(2048, dtype=np.float64)
    gx, gy = np.meshgrid(axis, axis)  # [row = y, col = x]
    xc, yc = gdcorr.gdcorr(gx, gy)
    lut = ((xc - gx).astype(np.float32), (yc - gy).astype(np.float32))
    for arr in lut:
        arr.flags.writeable = False
    return lut


# ---- frame transform ------------------------------------------------------------------------------

def frame_shift(info: FrameInfo, frame: str, times, track: PlanetTrack | None = None):
    """Shift (sx, sy) in hi-res pixels, in the image orientation, that takes a
    detector position to its position in the frame at the given times.

    The slew part is the user's rate [arcsec/day] * t, along (ANG - ORIENTAT) in the
    image frame. Planet frames also subtract the planet's motion relative to the moon.
    """
    times = np.asarray(times)
    if frame == 'image':
        return np.zeros(times.shape), np.zeros(times.shape)
    orient = math.radians(info.orientat) if info.orientat is not None else 0.0
    angle = math.radians(info.slew_angle) if info.slew_angle is not None else 0.0
    slew = info.slew_rate / 86400.0 / ARCSEC_PER_LOWRES_PIXEL * 2.0 * times  # hi-res pixels
    direction = angle - orient
    sx, sy = -slew * math.sin(direction), slew * math.cos(direction)
    if is_planet(frame):
        px, py = track.offset(times)
        sx, sy = sx - px, sy - py
    return sx, sy


def frame_coordinates(x, y, time, info: FrameInfo, frame: str, corr: Corrections, *,
                      seed: int = 0, geo_lut=None, track: PlanetTrack | None = None):
    """Per-event position (u, v) in a reference frame, in 0-indexed hi-res pixels (float32).

    Done EVENT_CHUNK events at a time, so the temporary arrays stay small. The random
    sub-pixel offsets are drawn per chunk from a generator seeded with (seed, first event), so
    repeated calls agree; a file of at most one chunk draws exactly as an unchunked run would.
    """
    n = len(x)
    u = np.empty(n, dtype=np.float32)
    v = np.empty(n, dtype=np.float32)
    for lo, hi in _chunks(n):
        rng = np.random.default_rng(seed if lo == 0 else [seed, lo]) if corr.randcorr else None
        u[lo:hi], v[lo:hi] = _frame_coordinates_block(
            x[lo:hi], y[lo:hi], time[lo:hi], info, frame, corr, rng, geo_lut, track)
    return u, v


def _frame_coordinates_block(x, y, time, info: FrameInfo, frame: str, corr: Corrections, rng,
                             geo_lut, track: PlanetTrack | None):
    """frame_coordinates for one chunk of events.

    Steps: (the flat-field weight is taken from the raw pixel when the file is loaded, so it
    comes first), random sub-pixel jitter, geometric correction, then the frame shift, then
    (north-up frames) a rotation by ORIENTAT about the detector centre. The distortion offset is
    looked up at each event's raw integer pixel, so adding it before or after the jitter gives
    the same position; the exposure-time division happens when the images are binned.
    """
    u = x.astype(np.float32) - 1.0  # AXIS1/AXIS2 are 1-indexed
    v = y.astype(np.float32) - 1.0
    if corr.geocorr:
        dxl, dyl = geo_lut
        ix = np.clip(x.astype(np.int32) - 1, 0, 2047)
        iy = np.clip(y.astype(np.int32) - 1, 0, 2047)
        u += dxl[iy, ix]
        v += dyl[iy, ix]
    if rng is not None:
        u += rng.random(len(u), dtype=np.float32) - 0.5
        v += rng.random(len(v), dtype=np.float32) - 0.5
    if frame == 'image':
        return u, v

    reason = info.unavailable_reason(frame)
    if reason:
        raise ValueError(f'{frame_label(info, frame)} unavailable: {reason}')
    sx, sy = frame_shift(info, frame, time, track)
    u = u + sx.astype(np.float32)
    v = v + sy.astype(np.float32)
    if is_north(frame):
        orient = math.radians(info.orientat)
        du, dv = u - DETECTOR_CENTRE, v - DETECTOR_CENTRE
        u = du * math.cos(orient) - dv * math.sin(orient) + DETECTOR_CENTRE
        v = du * math.sin(orient) + dv * math.cos(orient) + DETECTOR_CENTRE
    return u, v


@dataclass(eq=False)
class FrameGrid:
    """Image-bin grid for one frame; events are stored as flat bin indices."""

    bin_index: np.ndarray  # int32 per event: iy * nx + ix
    xlo: float             # left edge of the grid in frame coordinates
    ylo: float
    nx: int                # grid size in bins
    ny: int
    binsize: float = SPATIAL_BIN         # hi-res pixels for pixel frames, degrees for maps
    keep: np.ndarray | None = None       # maps: events that land on the planet (None = all)

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """(left, right, bottom, top) in frame coordinates (hi-res pixels, or map degrees)."""
        return (self.xlo, self.xlo + self.nx * self.binsize,
                self.ylo, self.ylo + self.ny * self.binsize)


def map_longitude_centre(orient: PlanetOrientation) -> float:
    """Circular mean of the sub-observer longitude: the longitude the rectilinear map is centred on."""
    lon = np.radians(orient.lon_obs)
    return float(np.degrees(np.arctan2(np.sin(lon).mean(), np.cos(lon).mean())) % 360.0)


def map_coordinates(ds: TagDataset, frame: str, corr: Corrections):
    """Per-event position in a map plane, and which events land on the planet.

    Steps: the pixel frame in which the planet is fixed (so the shifts, jitter and
    geometric correction are already applied), then the pixel -> sky -> screen conversion
    about the planet's centre, then the inverse projection onto the (raised) spheroid, then
    the map plane. Events are handled in time slices so the planet's orientation can be
    taken at that time; the sub-observer longitude advances linearly within each slice.
    """
    info, geom = ds.frame_info, ds.planet_geometry()
    base = info.map_base_frame()
    track = ds.planet_track() if is_planet(base) else None
    geo_lut = load_geo_lut() if corr.geocorr else None
    u, v = frame_coordinates(ds.x, ds.y, ds.time, info, base, corr,
                             seed=RANDOM_SEED, geo_lut=geo_lut, track=track)
    cx, cy = geom.centre(corr)
    a, b = spheroid_radii(geom.body, geom.altitude_km(corr))
    lon_centre = map_longitude_centre(geom.orient)

    out_u = np.zeros(len(u), dtype=np.float32)
    out_v = np.zeros(len(u), dtype=np.float32)
    keep = np.zeros(len(u), dtype=bool)
    edges = np.arange(0.0, ds.tmax + MAP_SLICE_SECONDS, MAP_SLICE_SECONDS)
    bounds = np.searchsorted(ds.time, edges)
    for lo, hi, t0, t1 in zip(bounds[:-1], bounds[1:], edges[:-1], edges[1:]):
        if hi <= lo:
            continue
        og = geom.orient.at(np.array([t0, t1]))
        frac = (ds.time[lo:hi].astype(np.float64) - t0) / (t1 - t0)
        lon_obs = og['lon_obs'][0] + frac * (og['lon_obs'][1] - og['lon_obs'][0])
        x, y = pixels_to_screen(u[lo:hi].astype(np.float64), v[lo:hi].astype(np.float64), cx, cy,
                                geom.cd, og['pole_pa'].mean(), og['radius_arcsec'].mean())
        lon, lat, mu, on = screen_to_surface(x, y, og['lat_obs_c'].mean(), lon_obs, a, b)
        mu_u, mu_v, in_map = map_plane_coordinates(lon, lat, frame, lon_centre)
        out_u[lo:hi], out_v[lo:hi] = mu_u, mu_v
        keep[lo:hi] = on & (mu >= MAP_MIN_MU) & in_map
    return out_u, out_v, keep


def _bin_events(u, v, xlo: float, ylo: float, binsize: float, keep=None):
    """(flat bin index per event, nx, ny) on a grid whose lower-left edge is (xlo, ylo). Events that are
    not `keep` go to bin 0 (and are ignored later through the grid's `keep` mask). Chunked: the
    per-event bin numbers are held as int16 until the grid width is known."""
    n = len(u)
    ix = np.zeros(n, dtype=np.int16)
    iy = np.zeros(n, dtype=np.int16)
    max_x = max_y = 0
    for lo, hi in _chunks(n):
        cx = np.floor((u[lo:hi] - xlo) / binsize)
        cy = np.floor((v[lo:hi] - ylo) / binsize)
        if keep is not None:
            cx, cy = np.where(keep[lo:hi], cx, 0), np.where(keep[lo:hi], cy, 0)
        cx, cy = cx.astype(np.int32), cy.astype(np.int32)
        max_x, max_y = max(max_x, int(cx.max())), max(max_y, int(cy.max()))
        if max(max_x, max_y) > np.iinfo(np.int16).max:
            raise ValueError('The image grid would be too large: check the reference frame and the planet-centre offsets.')
        ix[lo:hi], iy[lo:hi] = cx, cy
    nx, ny = max_x + 1, max_y + 1
    bin_index = np.empty(n, dtype=np.int32)
    for lo, hi in _chunks(n):
        bin_index[lo:hi] = iy[lo:hi].astype(np.int32) * nx + ix[lo:hi]
    return bin_index, nx, ny


def make_frame_grid(ds: TagDataset, frame: str, corr: Corrections) -> FrameGrid:
    """Bin all events in a frame. Bin edges sit at k*SPATIAL_BIN - 0.5 so a
    detector pixel (centred on its integer coordinate) never straddles a bin.
    Map frames are binned in degrees on a grid anchored at multiples of MAP_BIN_DEGREES.
    """
    if is_map(frame):
        u, v, keep = map_coordinates(ds, frame, corr)
        if not keep.any():
            raise ValueError('no events land on the planet in this projection '
                             '(check the planet-centre offsets and the emission altitude)')
        B = MAP_BIN_DEGREES
        xlo = math.floor(float(np.min(u, where=keep, initial=np.inf)) / B) * B
        ylo = math.floor(float(np.min(v, where=keep, initial=np.inf)) / B) * B
        bin_index, nx, ny = _bin_events(u, v, xlo, ylo, B, keep)
        bin_index.flags.writeable = False
        keep.flags.writeable = False
        return FrameGrid(bin_index, xlo, ylo, nx, ny, binsize=B, keep=keep)

    track = ds.planet_track() if is_planet(frame) else None
    geo_lut = load_geo_lut() if corr.geocorr else None
    u, v = frame_coordinates(ds.x, ds.y, ds.time, ds.frame_info, frame, corr,
                             seed=RANDOM_SEED, geo_lut=geo_lut, track=track)
    xlo = math.floor((float(u.min()) + 0.5) / SPATIAL_BIN) * SPATIAL_BIN - 0.5
    ylo = math.floor((float(v.min()) + 0.5) / SPATIAL_BIN) * SPATIAL_BIN - 0.5
    bin_index, nx, ny = _bin_events(u, v, xlo, ylo, SPATIAL_BIN)
    bin_index.flags.writeable = False
    return FrameGrid(bin_index, xlo, ylo, nx, ny)


# ---- exposure-time map (analytic) ---------------------------------------------------------------

def cumulative_exposure(gti: np.ndarray):
    """Function C(t): good exposure accumulated up to time t (piecewise linear)."""
    a, b = gti[:, 0], gti[:, 1]
    if len(a) == 1:
        a0, b0 = float(a[0]), float(b[0])
        return lambda t: np.clip(t, a0, b0) - a0
    knots = np.empty(2 * len(a))
    knots[0::2], knots[1::2] = a, b
    total = np.cumsum(b - a)
    values = np.empty_like(knots)
    values[0::2], values[1::2] = np.concatenate(([0.0], total[:-1])), total
    return lambda t: np.interp(t, knots, values)


def _axis_interval(a, w, d0, d1):
    """Times tau where d0 <= a - w*tau <= d1: the interval [lo, hi] (empty if lo > hi)."""
    if abs(w) < 1e-12:  # not moving along this axis: always or never inside
        inside = (a >= d0) & (a <= d1)
        return np.where(inside, -np.inf, np.inf), np.where(inside, np.inf, -np.inf)
    t1, t2 = (a - d1) / w, (a - d0) / w
    return np.minimum(t1, t2), np.maximum(t1, t2)


def compute_exposure_map(ds: TagDataset, grid: FrameGrid, frame: str, window,
                         track: PlanetTrack | None = None, geo_lut=None) -> np.ndarray:
    """Good exposure time (s) during which each frame pixel was on the detector.

    The field of view moves across the frame as the shift grows with time. For a
    straight-line shift the times a pixel is on the detector are an interval, found
    in closed form per axis and intersected with the good-time intervals through
    the cumulative exposure function, so there is no event simulation. The planet's
    curved path is split into PLANET_SEGMENTS straight segments. Each bin is sampled
    on an EXPOSURE_SUBSAMPLE x EXPOSURE_SUBSAMPLE grid, which turns the hard
    detector edge into fractional coverage.

    With the geometric correction on (geo_lut given) the detector's edges are
    curved in the corrected frame by up to ~13 pixels. The distortion varies slowly
    (about 1% per pixel), so it is taken as constant over each straight segment, at
    the segment's mid-time, which shifts each bin's footprint test by the local offset.
    """
    info, B = ds.frame_info, SPATIAL_BIN
    lo = max(window[0], float(ds.gti[0, 0]))
    hi = min(window[1], float(ds.gti[-1, 1]))
    if hi <= lo:
        return np.zeros((grid.ny, grid.nx))

    d0, d1 = DETECTOR_LIMITS
    nseg = PLANET_SEGMENTS if is_planet(frame) else 1
    tk = np.linspace(lo, hi, nseg + 1)
    sx, sy = frame_shift(info, frame, tk, track)
    # A footprint that never moves has hard edges cutting through bins, so sample
    # it more finely; a moving one is smoothed by the time average anyway.
    moving = np.ptp(sx) + np.ptp(sy) > B
    nsub = EXPOSURE_SUBSAMPLE if moving else 2 * EXPOSURE_SUBSAMPLE
    offsets = ((np.arange(nsub) + 0.5) / nsub - 0.5) * B
    qx = (grid.xlo + (np.arange(grid.nx) + 0.5)[:, None] * B + offsets[None, :]).ravel()
    qy = (grid.ylo + (np.arange(grid.ny) + 0.5)[:, None] * B + offsets[None, :]).ravel()
    cum = cumulative_exposure(ds.gti)
    if is_north(frame):
        th = math.radians(info.orientat)
        cos_t, sin_t = math.cos(th), math.sin(th)

    out = np.empty((qy.size, qx.size))
    chunk = 256  # rows at a time, to bound memory
    for r0 in range(0, qy.size, chunk):
        qyc = qy[r0:r0 + chunk]
        if is_north(frame):  # undo the rotation to get back to image orientation
            p0x = (cos_t * (qx - DETECTOR_CENTRE))[None, :] + (sin_t * (qyc - DETECTOR_CENTRE))[:, None] + DETECTOR_CENTRE
            p0y = (-sin_t * (qx - DETECTOR_CENTRE))[None, :] + (cos_t * (qyc - DETECTOR_CENTRE))[:, None] + DETECTOR_CENTRE
        else:
            p0x, p0y = qx[None, :], qyc[:, None]
        acc = np.zeros((qyc.size, qx.size))
        for j in range(nseg):
            dt = tk[j + 1] - tk[j]
            wx, wy = (sx[j + 1] - sx[j]) / dt, (sy[j + 1] - sy[j]) / dt
            start_x, start_y = p0x - sx[j], p0y - sy[j]
            if geo_lut is not None:
                smx, smy = frame_shift(info, frame, 0.5 * (tk[j] + tk[j + 1]), track)
                ix = np.clip(np.rint(p0x - smx).astype(np.intp), 0, 2047)
                iy = np.clip(np.rint(p0y - smy).astype(np.intp), 0, 2047)
                start_x = start_x - geo_lut[0][iy, ix]  # raw position = corrected - distortion
                start_y = start_y - geo_lut[1][iy, ix]
            lo_x, hi_x = _axis_interval(start_x, wx, d0, d1)
            lo_y, hi_y = _axis_interval(start_y, wy, d0, d1)
            t_lo = np.maximum(np.maximum(lo_x, lo_y), 0.0)
            t_hi = np.minimum(np.minimum(hi_x, hi_y), dt)
            acc += np.where(t_hi > t_lo, cum(tk[j] + t_hi) - cum(tk[j] + t_lo), 0.0)
        out[r0:r0 + chunk] = acc
    return out.reshape(grid.ny, nsub, grid.nx, nsub).mean(axis=(1, 3))


def _interval_vec(a, w, d0, d1):
    """Vectorised _axis_interval: times tau with d0 <= a - w*tau <= d1, for arrays a and w."""
    small = np.abs(w) < 1e-12
    w_safe = np.where(small, 1.0, w)
    t1, t2 = (a - d1) / w_safe, (a - d0) / w_safe
    inside = (a >= d0) & (a <= d1)
    lo = np.where(small, np.where(inside, -np.inf, np.inf), np.minimum(t1, t2))
    hi = np.where(small, np.where(inside, np.inf, -np.inf), np.maximum(t1, t2))
    return lo, hi


def map_knots(ds: TagDataset, grid: FrameGrid, frame: str, corr: Corrections,
              track: PlanetTrack | None = None, geo_lut=None) -> dict:
    """Where every map bin falls on the detector, and its foreshortening, at MAP_SEGMENTS + 1
    times spanning the good exposure. Independent of any time window, so it is cached and
    each new time selection only needs map_exposure_from_knots.

    Per bin: detector position (px, py; undoing the frame shift and the geometric
    correction), mu = cos(emission angle), and the surface area in km^2.
    """
    info, geom = ds.frame_info, ds.planet_geometry()
    base = info.map_base_frame()
    body = geom.body
    altitude = geom.altitude_km(corr)
    a, b = spheroid_radii(body, altitude)

    centres_u = grid.xlo + (np.arange(grid.nx) + 0.5) * grid.binsize
    centres_v = grid.ylo + (np.arange(grid.ny) + 0.5) * grid.binsize
    gu, gv = np.meshgrid(centres_u, centres_v)
    gu, gv = gu.ravel(), gv.ravel()
    lon, lat = map_plane_to_surface(gu, gv, frame)
    _, _, in_map = map_plane_coordinates(lon, lat, frame, map_longitude_centre(geom.orient))
    in_map &= np.abs(lat) <= 90.0
    area_km2 = np.where(in_map, map_bin_area(gu, gv, frame, grid.binsize,
                                             body.radius_eq_km + altitude, body.radius_pol_km + altitude), 0.0)

    tk = np.linspace(float(ds.gti[0, 0]), float(ds.gti[-1, 1]), MAP_SEGMENTS + 1)
    og = geom.orient.at(tk)
    sx, sy = frame_shift(info, base, tk, track)
    cx, cy = geom.centre(corr)
    px = np.empty((len(tk), gu.size), dtype=np.float32)
    py = np.empty_like(px)
    mu = np.empty_like(px)
    for k in range(len(tk)):
        x, y, m = surface_to_screen(lon, lat, og['lat_obs_c'][k], og['lon_obs'][k], a, b)
        ru, rv = screen_to_pixels(x, y, cx, cy, geom.cd_inv, og['pole_pa'][k], og['radius_arcsec'][k])
        pu, pv = ru - sx[k], rv - sy[k]  # undo the frame shift: detector position
        if geo_lut is not None:
            ix = np.clip(np.rint(pu).astype(np.intp), 0, 2047)
            iy = np.clip(np.rint(pv).astype(np.intp), 0, 2047)
            pu, pv = pu - geo_lut[0][iy, ix], pv - geo_lut[1][iy, ix]
        px[k], py[k], mu[k] = pu, pv, m
    return dict(tk=tk, px=px, py=py, mu=mu, area_km2=area_km2, shape=(grid.ny, grid.nx),
                km_per_arcsec=body.radius_eq_km / og['radius_arcsec'], pixel_area_arcsec2=geom.pixel_area_arcsec2)


def map_exposure_from_knots(knots: dict, gti: np.ndarray, window):
    """(pixel_seconds, full) per map bin for the good time inside `window`.

    A map bin is seen by a number of detector pixels that depends on how foreshortened it
    is (mu) and on whether it is on the detector at all. Mean surface brightness is the
    counts divided by the pixel-seconds of coverage:

        pixel_seconds = integral over good time of  mu * (bin area / pixel area) * [on the detector]
        full          = the same without the on-the-detector condition (still mu >= MAP_MIN_MU)

    Between the knot times each bin's detector position and mu are straight lines (the
    same trick as the planet's path in pixel frames), so the time on the detector inside
    each segment follows in closed form and is intersected with the good-time intervals.
    """
    tk, px, py, mu = knots['tk'], knots['px'], knots['py'], knots['mu']
    nb = px.shape[1]
    cum = cumulative_exposure(gti)
    d0, d1 = DETECTOR_LIMITS
    pixel_seconds = np.zeros(nb)
    full = np.zeros(nb)
    for j in range(len(tk) - 1):
        dt = tk[j + 1] - tk[j]
        w_lo, w_hi = max(window[0] - tk[j], 0.0), min(window[1] - tk[j], dt)  # the window inside this segment
        if w_hi <= w_lo:
            continue
        kmas = 0.5 * (knots['km_per_arcsec'][j] + knots['km_per_arcsec'][j + 1])
        pixels_face_on = knots['area_km2'] / (kmas ** 2 * knots['pixel_area_arcsec2'])  # pixels covered if mu = 1
        mu_a, mu_b = mu[j].astype(np.float64), mu[j + 1].astype(np.float64)
        f_lo, f_hi = _interval_vec(mu_a, -(mu_b - mu_a) / dt, MAP_MIN_MU, 10.0)
        f_lo, f_hi = np.maximum(f_lo, w_lo), np.minimum(f_hi, w_hi)
        mid = 0.5 * (f_lo + f_hi)
        full += np.where(f_hi > f_lo, (cum(tk[j] + f_hi) - cum(tk[j] + f_lo))
                         * (mu_a + (mu_b - mu_a) * mid / dt) * pixels_face_on, 0.0)
        xa, xb = px[j].astype(np.float64), px[j + 1].astype(np.float64)
        ya, yb = py[j].astype(np.float64), py[j + 1].astype(np.float64)
        x_lo, x_hi = _interval_vec(xa, -(xb - xa) / dt, d0, d1)
        y_lo, y_hi = _interval_vec(ya, -(yb - ya) / dt, d0, d1)
        t_lo = np.maximum(np.maximum(f_lo, x_lo), y_lo)
        t_hi = np.minimum(np.minimum(f_hi, x_hi), y_hi)
        mid = 0.5 * (t_lo + t_hi)
        pixel_seconds += np.where(t_hi > t_lo, (cum(tk[j] + t_hi) - cum(tk[j] + t_lo))
                                  * (mu_a + (mu_b - mu_a) * mid / dt) * pixels_face_on, 0.0)
    ny, nx = knots['shape']
    return pixel_seconds.reshape(ny, nx), full.reshape(ny, nx)


def _compute_map_exposure(ds: TagDataset, grid: FrameGrid, frame: str, window, corr: Corrections,
                          track: PlanetTrack | None = None, geo_lut=None):
    """Both steps in one call (the dataset caches the first and reuses it across windows)."""
    return map_exposure_from_knots(map_knots(ds, grid, frame, corr, track, geo_lut), ds.gti, window)


def compute_map_exposure(*args, **kwargs):
    """See map_exposure_from_knots (inf - inf in the interval arithmetic is expected and harmless)."""
    with np.errstate(invalid='ignore', over='ignore'):
        return _compute_map_exposure(*args, **kwargs)


GRATICULE_STEP = 10  # degrees between graticule lines drawn on the pixel frames


def _rotate_north(u, v, info: FrameInfo):
    """Rotate image-orientation frame coordinates by ORIENTAT about the detector centre."""
    orient = math.radians(info.orientat)
    du, dv = u - DETECTOR_CENTRE, v - DETECTOR_CENTRE
    return (du * math.cos(orient) - dv * math.sin(orient) + DETECTOR_CENTRE,
            du * math.sin(orient) + dv * math.cos(orient) + DETECTOR_CENTRE)


def _visible_runs(x, y, visible):
    """Split a polyline into the runs where `visible` is True."""
    runs, start = [], None
    for i, ok in enumerate(visible):
        if ok and start is None:
            start = i
        elif not ok and start is not None:
            if i - start > 1:
                runs.append((x[start:i], y[start:i]))
            start = None
    if start is not None and len(x) - start > 1:
        runs.append((x[start:], y[start:]))
    return runs


def planet_overlay(ds: TagDataset, frame: str, corr: Corrections):
    """The planet's limb and a lon/lat graticule in the coordinates of a pixel frame in
    which the planet is fixed (for registering the planet centre), or None if not applicable.

    Drawn at mid-exposure. The limb is shown at the 1-bar level and at the emission
    altitude; the graticule is on the emission-altitude surface the maps use.
    """
    info = ds.frame_info
    if info.body is None or is_map(frame) or frame == 'image':
        return None
    base = info.map_base_frame()
    if frame not in (base, base + '_north'):
        return None
    geom = ds.planet_geometry()
    og = {k: float(v) for k, v in geom.orient.at(0.5 * ds.tmax).items()}
    cx, cy = geom.centre(corr)
    altitude = geom.altitude_km(corr)

    def to_frame(x, y):
        u, v = screen_to_pixels(x, y, cx, cy, geom.cd_inv, og['pole_pa'], og['radius_arcsec'])
        return _rotate_north(u, v, info) if is_north(frame) else (u, v)

    def to_display(u, v):
        return _rotate_north(u, v, info) if is_north(frame) else (u, v)

    lines = []
    for kind, h in (('limb', 0.0), ('limb_alt', altitude)):
        x, y = limb_outline(og['lat_obs_c'], *spheroid_radii(geom.body, h))
        u, v = to_frame(x, y)
        lines.append({'kind': kind, 'x': u, 'y': v})
    a, b = spheroid_radii(geom.body, altitude)
    sweeps = [(np.full(361, float(lat)), np.linspace(0.0, 360.0, 361)) for lat in range(-80, 81, GRATICULE_STEP)]
    sweeps += [(np.linspace(-90.0, 90.0, 181), np.full(181, float(lon))) for lon in range(0, 360, GRATICULE_STEP)]
    for lat, lon in sweeps:
        x, y, mu = surface_to_screen(lon, lat, og['lat_obs_c'], og['lon_obs'], a, b)
        u, v = to_frame(x, y)
        for ru, rv in _visible_runs(u, v, mu > 0.0):
            lines.append({'kind': 'grid', 'x': ru, 'y': rv})
    # Where the registration offsets point: +x and +y of the planet-frame pixels (the Centre
    # delta-x / delta-y fields) as displayed, which is rotated in the north-up frame.
    c0, cx1, cy1 = to_display(cx, cy), to_display(cx + 1.0, cy), to_display(cx, cy + 1.0)
    lines.append({'kind': 'compass', 'centre': c0,
                  'x_dir': (cx1[0] - c0[0], cx1[1] - c0[1]),
                  'y_dir': (cy1[0] - c0[0], cy1[1] - c0[1])})
    return lines


def available_frames(ds: TagDataset) -> dict[str, str]:
    """{frame key: menu label} for the frames this dataset supports."""
    return {k: frame_label(ds.frame_info, k) for k in FRAMES
            if ds.frame_info.unavailable_reason(k) is None}


def counts_label(ds: TagDataset, corr: Corrections) -> str:
    return 'flat-fielded counts' if (corr.flatcorr and ds.weight is not None) else 'counts'
