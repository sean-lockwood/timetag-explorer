"""The event list of one exposure (pure numpy): reading a TIME-TAG file, slicing in time, binning into frames,
flat-field weights and exposure maps. No web server and no plotting."""
from __future__ import annotations

import threading
from collections import OrderedDict
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
from astropy.io import fits

from . import config
from .constants import (EDGE_TRIM_PIXELS, EVENT_CHUNK, FLAT_SHAPE, FLAT_URL, MAX_CACHED_FRAMES, SPATIAL_BIN,
                        TIME_CORR_MIN_FRACTION, TIME_COL, X_COL, Y_COL, TimeSel, _chunks)
from .frames import (DEFAULT_CORRECTIONS, Corrections, FrameInfo, clamp_binsize, frame_label, is_map, is_planet,
                     load_frame_info)
from .horizons import (PlanetTrack, load_horizons_positions, planet_track_from_offsets,
                       sky_offsets_from_positions)
from .log import logger
from .regions import Region
from .transforms import (FrameGrid, PlanetGeometry, build_planet_geometry, compute_exposure_map, load_geo_lut,
                         make_frame_grid, map_exposure_from_knots, map_knots)



# --------------------------------------------------------------------------
# Data layer (no UI, no matplotlib)
# --------------------------------------------------------------------------

@dataclass(eq=False)
class TagDataset:
    """Immutable, time-sorted photon event list.

    Positions are kept as the recorded 1-indexed pixel coordinates; each
    reference frame bins them on demand into a FrameGrid (cached, a few
    frames per dataset). Each event also carries its flat-field weight (the
    flat file's value at its *image-frame* pixel) when the flat is available,
    so images and light curves are weighted sums whatever frame is displayed.
    """

    name: str
    time: np.ndarray       # float32, seconds, sorted ascending
    x: np.ndarray          # int16 AXIS1 (1-indexed hi-res pixel)
    y: np.ndarray          # int16 AXIS2
    gti: np.ndarray        # (n, 2) good-time intervals [start, stop] in seconds
    frame_info: FrameInfo  # slew / WCS description for the reference frames
    weight: np.ndarray | None = None     # float32 per-event flat-field weight
    flat_problem: str | None = None      # why there is no flat, if there is none
    rootname: str = ''     # header metadata used for website links
    proposid: int | None = None
    obset_id: str | None = None
    targname: str = ''
    linenum: str = ''
    tdateobs: str = ''     # start date / time (UTC) and exposure time, for display
    ttimeobs: str = ''
    texptime: float | None = None
    propaper: str = ''     # proposed aperture / filter
    _grids: OrderedDict = field(default_factory=OrderedDict, init=False, repr=False)
    _expo: OrderedDict = field(default_factory=OrderedDict, init=False, repr=False)
    _track: PlanetTrack | None = field(default=None, init=False, repr=False)
    _positions: np.ndarray | None = field(default=None, init=False, repr=False)
    _geom: PlanetGeometry | None = field(default=None, init=False, repr=False)
    _valid: np.ndarray | None = field(default=None, init=False, repr=False)
    _cum: tuple | None = field(default=None, init=False, repr=False)  # (key, running sum)
    _knots: tuple | None = field(default=None, init=False, repr=False)  # (key, map_knots result)
    _lock: threading.RLock = field(default_factory=threading.RLock, init=False, repr=False)

    @property
    def tmax(self) -> float:
        return float(self.time[-1])

    def positions(self) -> np.ndarray:
        """Planet (and moon) positions from the saved Horizons file or a fresh query (once)."""
        with self._lock:
            if self._positions is None:
                self._positions = load_horizons_positions(self.frame_info, self.rootname)
            return self._positions

    def planet_track(self) -> PlanetTrack:
        """The planet's motion relative to the moon (fetched once, then kept)."""
        with self._lock:
            if self._track is None:
                t, d_xi, d_eta = sky_offsets_from_positions(self.positions(), self.frame_info.expstart)
                self._track = planet_track_from_offsets(t, d_xi, d_eta, self.frame_info.cd)
            return self._track

    def planet_geometry(self) -> PlanetGeometry:
        """Orientation, scale and centre guess for the planet maps (built once)."""
        with self._lock:
            if self._geom is None:
                self._geom = build_planet_geometry(self.frame_info, self.positions())
            return self._geom

    def prepare_overlay(self) -> None:
        """Load what the limb/graticule overlay needs (may use the network); ValueError if it can't be shown."""
        reason = self.frame_info.unavailable_reason('map_rect')
        if reason:
            raise ValueError(f'planet limb unavailable: {reason}')
        self.planet_geometry()

    def prepare_frame(self, frame: str) -> None:
        """Check a frame is available and load what it needs (may use the network).
        Raises ValueError with a reason if it cannot be built."""
        reason = self.frame_info.unavailable_reason(frame)
        if reason:
            raise ValueError(f'{frame_label(self.frame_info, frame)} unavailable: {reason}')
        if is_planet(frame) or (is_map(frame) and is_planet(self.frame_info.map_base_frame())):
            self.planet_track()
        if is_map(frame):
            self.planet_geometry()

    def grid(self, frame: str = 'image', corr: Corrections = DEFAULT_CORRECTIONS) -> FrameGrid:
        """Bin grid for a frame (computed on first use, then cached)."""
        key = (frame, corr.randcorr, corr.geocorr, corr.map_key() if is_map(frame) else None)
        with self._lock:
            grid = self._grids.get(key)
            if grid is None:
                grid = self._grids[key] = make_frame_grid(self, frame, corr)
                while len(self._grids) > MAX_CACHED_FRAMES:
                    self._grids.popitem(last=False)  # drop the least recently used
            else:
                self._grids.move_to_end(key)
            return grid

    def exposure(self, t0: float = -np.inf, t1: float = np.inf) -> float:
        """Good exposure time (s) inside [t0, t1]: only time within a GTI row counts."""
        a, b = self.gti[:, 0], self.gti[:, 1]
        return float(np.sum(np.clip(b, t0, t1) - np.clip(a, t0, t1)))

    def bin_exposure(self, edges: np.ndarray) -> np.ndarray:
        """Good exposure time (s) in each bin [edges[i], edges[i+1])."""
        a, b = self.gti[:, 0], self.gti[:, 1]
        accumulated = (np.clip(edges[:, None], a, b) - a).sum(axis=1)  # exposure up to each edge
        return np.diff(accumulated)

    def exposure_map(self, time_sel: TimeSel | None, frame: str = 'image',
                     corr: Corrections = DEFAULT_CORRECTIONS) -> np.ndarray:
        """(ny, nx) good exposure time each frame pixel was on the detector in the window."""
        window = (-np.inf, np.inf) if time_sel is None else tuple(time_sel)
        key = (frame, corr.randcorr, corr.geocorr, window)
        grid = self.grid(frame, corr)
        with self._lock:
            emap = self._expo.get(key)
            if emap is not None:
                self._expo.move_to_end(key)
                return emap
        track = self.planet_track() if is_planet(frame) else None
        geo_lut = load_geo_lut() if corr.geocorr else None
        emap = compute_exposure_map(self, grid, frame, window, track, geo_lut)
        emap.flags.writeable = False
        with self._lock:
            self._expo[key] = emap
            while len(self._expo) > 6:
                self._expo.popitem(last=False)
        return emap

    def valid_events(self, corr: Corrections) -> np.ndarray | None:
        """Boolean per-event mask of the events kept, or None to keep all.

        With the exposure-time correction on, events within EDGE_TRIM_PIXELS of a
        detector edge (vignetted, so not a perfect rectangle) are discarded. This
        uses the recorded pixel, i.e. the detector position, whatever the frame.
        """
        if not corr.timecorr or EDGE_TRIM_PIXELS <= 0:
            return None
        with self._lock:
            if self._valid is None:
                lo, hi = EDGE_TRIM_PIXELS + 1, 2048 - EDGE_TRIM_PIXELS  # AXIS values are 1-indexed
                valid = (self.x >= lo) & (self.x <= hi) & (self.y >= lo) & (self.y <= hi)
                valid.flags.writeable = False
                self._valid = valid
            return self._valid

    def _running_sum(self, weighted: bool, valid: np.ndarray | None) -> np.ndarray:
        """[0, cumsum] over events of the flat weight (or 1) of the kept events, so any
        time range sums in O(1). One is kept at a time (8 bytes per event); switching
        corrections rebuilds it in tens of milliseconds per million events."""
        key = (weighted, valid is not None)
        with self._lock:
            if self._cum is None or self._cum[0] != key:
                n = len(self.time)
                running = np.zeros(n + 1)  # float64, filled a chunk at a time
                for lo, hi in _chunks(n):
                    if weighted:
                        per_event = self.weight[lo:hi] * valid[lo:hi] if valid is not None else self.weight[lo:hi]
                    else:
                        per_event = valid[lo:hi]
                    running[lo + 1:hi + 1] = np.cumsum(per_event, dtype=np.float64) + running[lo]
                running.flags.writeable = False
                self._cum = (key, running)
            return self._cum[1]

    def map_exposure(self, time_sel: TimeSel | None, frame: str, corr: Corrections):
        """(pixel_seconds, full) per bin of a map frame for the time window (cached)."""
        window = (-np.inf, np.inf) if time_sel is None else tuple(time_sel)
        key = ('map', frame, corr.randcorr, corr.geocorr, corr.map_key(), window)
        grid = self.grid(frame, corr)
        with self._lock:
            hit = self._expo.get(key)
            if hit is not None:
                self._expo.move_to_end(key)
                return hit
        knots_key = (frame, corr.randcorr, corr.geocorr, corr.map_key())
        with self._lock:
            if self._knots is None or self._knots[0] != knots_key:
                track = self.planet_track() if is_planet(self.frame_info.map_base_frame()) else None
                geo_lut = load_geo_lut() if corr.geocorr else None
                self._knots = (knots_key, map_knots(self, grid, frame, corr, track, geo_lut))
            knots = self._knots[1]
        with np.errstate(invalid='ignore', over='ignore'):  # inf - inf in the interval arithmetic is harmless
            result = map_exposure_from_knots(knots, self.gti, window)
        with self._lock:
            self._expo[key] = result
            while len(self._expo) > 6:
                self._expo.popitem(last=False)
        return result

    def time_slice(self, time_sel: TimeSel | None) -> slice:
        """Events are time-sorted, so a time cut is a cheap slice (no copy)."""
        if time_sel is None:
            return slice(None)
        lo, hi = np.searchsorted(self.time, time_sel)
        return slice(int(lo), int(hi))

    def _weighted(self, corr: Corrections) -> bool:
        return corr.flatcorr and self.weight is not None

    def region_grid(self, region: Region, frame: str = 'image',
                    corr: Corrections = DEFAULT_CORRECTIONS) -> np.ndarray:
        """Boolean (ny, nx) image grid: True where the bin centre is in the region."""
        grid = self.grid(frame, corr)
        xc = grid.xlo + (np.arange(grid.nx) + 0.5) * grid.binsize
        yc = grid.ylo + (np.arange(grid.ny) + 0.5) * grid.binsize
        gx, gy = np.meshgrid(xc, yc)  # shape (ny, nx), row-major like bin_index
        return region.contains(gx.ravel(), gy.ravel()).reshape(grid.ny, grid.nx)

    def region_mask(self, region: Region, frame: str = 'image',
                    corr: Corrections = DEFAULT_CORRECTIONS) -> np.ndarray:
        """Per-event boolean mask for a region drawn in a frame.

        The shape is tested on the image-bin grid (a bin is in if its centre
        is inside), so the cut matches exactly what the user sees highlighted,
        and the cost depends on the image size, not the number of events.
        """
        grid = self.grid(frame, corr)
        in_bin = self.region_grid(region, frame, corr).ravel()
        mask = np.empty(len(grid.bin_index), dtype=bool)
        for lo, hi in _chunks(len(mask)):  # (an int32 index array is copied to int64 when indexing)
            mask[lo:hi] = in_bin[grid.bin_index[lo:hi]]
        if grid.keep is not None:  # maps: events that miss the planet are in no region
            mask &= grid.keep
        return mask

    def light_curve(self, binsize: float, region: Region | None = None,
                    frame: str = 'image', corr: Corrections = DEFAULT_CORRECTIONS):
        """Return (bin edges, count rate). Optionally restricted to a region.

        Rates divide by the *good* exposure in each bin (GTI), so a bin that
        straddles a gap is normalised correctly and a bin entirely inside a
        gap is NaN. The region is interpreted in the given frame. With the
        exposure-time correction on, the vignetted detector border is excluded
        (as in the images); its per-pixel exposure map is not applied here.
        """
        binsize = clamp_binsize(binsize, self.tmax)  # whatever the caller passed, the arrays stay small
        nbins = int(np.floor(self.tmax / binsize)) + 1
        edges = np.arange(nbins + 1) * binsize
        weighted = self._weighted(corr)
        valid = self.valid_events(corr)
        if region is None:
            # Sorted times: the histogram is the difference of the running
            # sum (of events, or of weights) at the bin edges.
            idx = np.searchsorted(self.time, edges)
            if weighted or valid is not None:
                counts = np.diff(self._running_sum(weighted, valid)[idx])
            else:
                counts = np.diff(idx)
        else:
            mask = self.region_mask(region, frame, corr)
            if valid is not None:
                mask &= valid
            # float64 weights: np.histogram accumulates in the weights' dtype,
            # and float32 sums over ~1e6 events lose ~1e-6 of precision.
            counts = np.zeros(nbins)
            for lo, hi in _chunks(len(mask)):
                chosen = mask[lo:hi]
                if not chosen.any():
                    continue
                w = self.weight[lo:hi][chosen].astype(np.float64) if weighted else None
                counts += np.histogram(self.time[lo:hi][chosen], bins=edges, weights=w)[0]
        exposure = self.bin_exposure(edges)
        rate = np.divide(counts, exposure, out=np.full(len(counts), np.nan), where=exposure > 1e-9)
        return edges, rate

    def image(self, time_sel: TimeSel | None = None, frame: str = 'image',
              corr: Corrections = DEFAULT_CORRECTIONS) -> np.ndarray:
        """Count-rate image (ny, nx) in a frame.

        With the exposure-time correction, events in the vignetted border of the
        detector are discarded and each bin is divided by the time it was on the
        rest of the detector (so the edges of a drifting field of view are not dim);
        bins seen for less than TIME_CORR_MIN_FRACTION of the good exposure are
        left at 0. Without it, every bin is divided by the total good exposure.
        """
        grid = self.grid(frame, corr)
        start, stop, _ = self.time_slice(time_sel).indices(len(self.time))
        weighted = self._weighted(corr)
        keep = self.valid_events(corr)       # the trimmed detector border (exposure-time correction)
        if grid.keep is not None:            # maps: only events that land on the planet
            keep = grid.keep if keep is None else keep & grid.keep
        counts = np.zeros(grid.nx * grid.ny)
        for lo in range(start, stop, EVENT_CHUNK):  # a chunk at a time: bincount copies its inputs to int64 / float64
            hi = min(lo + EVENT_CHUNK, stop)
            w = self.weight[lo:hi] if weighted else None
            if keep is not None:  # a zero weight discards the event
                w = keep[lo:hi].astype(np.float64) if w is None else w * keep[lo:hi]
            counts += np.bincount(grid.bin_index[lo:hi], weights=w, minlength=grid.nx * grid.ny)
        t0, t1 = (-np.inf, np.inf) if time_sel is None else time_sel
        total = self.exposure(t0, t1)
        if total <= 0:
            return np.zeros((grid.ny, grid.nx))
        if is_map(frame):
            # Brightness per 4x4 hi-res pixels: counts over the pixel-seconds that saw each bin.
            pixel_seconds, full = self.map_exposure(time_sel, frame, corr)
            denom = (pixel_seconds if corr.timecorr else full).ravel()
            ok = (full.ravel() > 0) & (denom >= TIME_CORR_MIN_FRACTION * full.ravel())
            rate = np.divide(counts * SPATIAL_BIN ** 2, denom, out=np.zeros(counts.shape), where=ok)
            return rate.reshape(grid.ny, grid.nx)
        if corr.timecorr:
            emap = self.exposure_map(time_sel, frame, corr).ravel()
            rate = np.divide(counts, emap, out=np.zeros(counts.shape),
                             where=emap >= TIME_CORR_MIN_FRACTION * total)
        else:
            rate = counts / total
        return rate.reshape(grid.ny, grid.nx)


def _as_float(value) -> float | None:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def _read_gti(path: Path, tmax: float) -> np.ndarray:
    """Good-time intervals (ext 'GTI') as an (n, 2) array; whole exposure if absent."""
    gti = np.empty((0, 2))
    try:
        table = fits.getdata(path, extname='GTI')
        gti = np.column_stack((np.asarray(table['START'], dtype=np.float64),
                               np.asarray(table['STOP'], dtype=np.float64)))
        gti = gti[gti[:, 1] > gti[:, 0]]
    except KeyError:
        logger.warning('%s has no GTI extension; assuming one interval [0, tmax]', path.name)
    if gti.size == 0:
        gti = np.array([[0.0, tmax]])
    return gti[np.argsort(gti[:, 0])]


@lru_cache(maxsize=1)
def load_flat_weights() -> np.ndarray:
    """Per-pixel weight map flat as float32, shared by all sessions.

    Pixels with a non-finite or non-positive flat value get weight 0, so
    those events drop out instead of becoming inf/NaN.
    """
    if not config.FLAT_FILE.is_file():
        raise FileNotFoundError(f'Flat field not found: {config.FLAT_FILE} (it can be downloaded from {FLAT_URL})')
    flat = np.asarray(fits.getdata(config.FLAT_FILE, ext=0), dtype=np.float64)
    if flat.shape != FLAT_SHAPE:
        raise ValueError(f'Flat field has shape {flat.shape}, expected {FLAT_SHAPE}')
    weights = np.zeros(flat.shape, dtype=np.float32)
    ok = np.isfinite(flat) & (flat > 0)
    weights[ok] = flat[ok]
    weights.flags.writeable = False
    return weights


def _is_sorted(a: np.ndarray) -> bool:
    """True if `a` never decreases (checked a chunk at a time, without a temporary of the whole array)."""
    for lo, hi in _chunks(len(a)):
        seg = a[max(lo - 1, 0):hi]
        if np.any(seg[1:] < seg[:-1]):
            return False
    return True


def _load_dataset_uncached(path: str) -> TagDataset:
    """Read a TAG file. Blocking: call via asyncio.to_thread."""
    tag_path = Path(path)
    data = fits.getdata(tag_path, ext=1)

    missing = [c for c in (TIME_COL, X_COL, Y_COL) if c not in data.names]
    if missing:
        raise ValueError(f'Missing column(s) in event table: {", ".join(missing)}')
    if len(data) == 0:
        raise ValueError('Event table is empty')

    # FITS arrays are big-endian and the table is memory-mapped: copy each column straight into a
    # compact native dtype (no wider temporaries), then let go of the mapping.
    time = np.asarray(data[TIME_COL], dtype=np.float32)
    x = np.asarray(data[X_COL], dtype=np.int16)
    y = np.asarray(data[Y_COL], dtype=np.int16)
    del data

    if not _is_sorted(time):
        order = np.argsort(time, kind='stable')
        time, x, y = time[order], x[order], y[order]
        del order

    weight = flat_problem = None
    try:
        wmap = load_flat_weights()
    except (OSError, ValueError) as exc:  # no flat: carry on with unweighted counts
        flat_problem = str(exc)
        logger.warning('flat field unavailable (%s); counts will be unweighted', exc)
    else:
        weight = np.zeros(len(x), dtype=np.float32)
        outside = 0
        for lo, hi in _chunks(len(x)):
            xc, yc = x[lo:hi], y[lo:hi]
            inside = (xc >= 1) & (xc < wmap.shape[1] + 1) & (yc >= 1) & (yc < wmap.shape[0] + 1)
            weight[lo:hi][inside] = wmap[yc[inside] - 1, xc[inside] - 1]  # 1-indexed AXIS -> 0-indexed [row=AXIS2, col=AXIS1]
            outside += int((~inside).sum())
        if outside:
            logger.warning('%d events fall outside the flat field and get weight 0', outside)

    gti = _read_gti(tag_path, float(time[-1]))
    header = fits.getheader(tag_path, ext=0)
    frame_info = load_frame_info(tag_path, fits.getheader(tag_path, ext=1))

    for arr in (time, x, y, weight, gti):
        if arr is not None:
            arr.flags.writeable = False

    return TagDataset(
        name=tag_path.name, time=time, x=x, y=y, gti=gti, frame_info=frame_info,
        weight=weight, flat_problem=flat_problem,
        rootname=tag_path.name.split('_', 1)[0].lower(),
        proposid=header.get('PROPOSID'),
        obset_id=None if header.get('OBSET_ID') is None else str(header['OBSET_ID']).strip(),
        targname=str(header.get('TARGNAME', '')).strip(),
        linenum=str(header.get('LINENUM', '')).strip(),
        tdateobs=str(header.get('TDATEOBS', '')).strip(),
        ttimeobs=str(header.get('TTIMEOBS', '')).strip(),
        texptime=_as_float(header.get('TEXPTIME')),
        propaper=str(header.get('PROPAPER', '')).strip(),
    )


# --------------------------------------------------------------------------


def load_dataset(path) -> TagDataset:
    """Read a TIME-TAG file into a TagDataset. If `<rootname>_spt.fits` is in the same folder, the pointing
    information for the rest frames and planet maps is read from it too. The flat field (`config.FLAT_FILE`) is
    applied when it exists."""
    return _load_dataset_uncached(str(path))
