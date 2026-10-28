"""Planet and moon positions from JPL Horizons, saved in a local folder (`config.HORIZONS_DIR`) and reused."""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from . import config
from .frames import HORIZONS_IDS, FrameInfo
from .log import audit, logger



# ---- planet rest frame: the planet's motion relative to the moon ----------------------------

@dataclass(frozen=True)
class PlanetTrack:
    """Planet's displacement relative to the moon during the exposure, in hi-res
    pixels (image orientation), sampled at t seconds after the exposure start."""

    t: np.ndarray
    dx: np.ndarray
    dy: np.ndarray

    def offset(self, times):
        """Displacement since t = 0 at the given times (linear interpolation)."""
        x0, y0 = np.interp(0.0, self.t, self.dx), np.interp(0.0, self.t, self.dy)
        return np.interp(times, self.t, self.dx) - x0, np.interp(times, self.t, self.dy) - y0


def planet_track_from_offsets(t, d_xi, d_eta, cd) -> PlanetTrack:
    """Turn sky offsets of the planet minus the moon (arcsec east, arcsec north) into
    detector pixels with the inverse WCS CD matrix. Offsets are tiny next to the
    field, so this linearisation replaces a per-event pixel -> RA/Dec transform."""
    cd_inv = np.linalg.inv(np.array(cd, dtype=float).reshape(2, 2))
    dx, dy = cd_inv @ (np.vstack((d_xi, d_eta)) / 3600.0)  # degrees -> pixels
    return PlanetTrack(np.asarray(t, dtype=float), dx, dy)


def query_horizons(planet: str, moon: str | None, expstart: float, expend: float) -> np.ndarray:
    """Planet (and moon) positions as seen from HST, one per minute, from JPL Horizons.

    Returns an (n, 6) array: JD (UTC), planet RA, Dec (degrees), planet distance (AU),
    moon RA, Dec (degrees; NaN when there is no moon).
    """
    try:
        from astroquery.jplhorizons import Horizons
    except ImportError as exc:
        raise ValueError('planet frames and maps need astroquery (pip install astroquery) '
                         f'or a saved ephemeris file in {config.HORIZONS_DIR.name}/') from exc
    mjd0 = datetime(1858, 11, 17)
    start = (mjd0 + timedelta(days=expstart) - timedelta(minutes=1)).replace(second=0, microsecond=0)
    stop = (mjd0 + timedelta(days=expend) + timedelta(minutes=2)).replace(second=0, microsecond=0)
    epochs = {'start': start.strftime('%Y-%m-%d %H:%M'),
              'stop': stop.strftime('%Y-%m-%d %H:%M'), 'step': '1m'}

    def ask(name):
        return Horizons(id=HORIZONS_IDS[name], id_type=None, location='@hst',
                        epochs=epochs).ephemerides()

    started = time.time()
    try:
        planet_tab = ask(planet)
        moon_tab = ask(moon) if moon else None
        audit('horizons_ok', planet=planet, moon=moon, start=epochs['start'], stop=epochs['stop'],
              seconds=f'{time.time() - started:.1f}')
    except Exception as exc:  # network / service errors surface as many types
        audit('horizons_failed', logging.WARNING, planet=planet, moon=moon, error=exc)
        raise ValueError(f'JPL Horizons query failed: {exc}') from exc
    jd = np.asarray(planet_tab['datetime_jd'], dtype=float)
    try:
        delta = np.asarray(planet_tab['delta'], dtype=float)
    except KeyError as exc:
        raise ValueError(f'JPL Horizons table has no distance column; columns: {planet_tab.colnames}') from exc
    nan = np.full(len(jd), np.nan)
    ra_m, dec_m = nan, nan
    if moon_tab is not None:
        if len(jd) != len(moon_tab['datetime_jd']):
            raise ValueError('JPL Horizons returned tables of different lengths')
        ra_m, dec_m = (np.asarray(moon_tab[c], dtype=float) for c in ('RA', 'DEC'))
    return np.column_stack((jd, np.asarray(planet_tab['RA'], dtype=float),
                            np.asarray(planet_tab['DEC'], dtype=float), delta, ra_m, dec_m))


def sky_offsets_from_positions(positions: np.ndarray, expstart: float):
    """(t [s since EXPSTART], d_xi, d_eta [arcsec]) of planet minus moon.

    d_xi is east (RA difference times cos Dec, wrap-safe) and d_eta is north.
    """
    jd, ra_p, dec_p, _, ra_m, dec_m = positions.T
    d_ra = (ra_p - ra_m + 180.0) % 360.0 - 180.0           # degrees
    d_xi = d_ra * np.cos(np.radians(dec_m)) * 3600.0        # arcsec, east
    d_eta = (dec_p - dec_m) * 3600.0                        # arcsec, north
    return (jd - 2400000.5 - expstart) * 86400.0, d_xi, d_eta


def horizons_path(rootname: str, planet: str, moon: str | None = None) -> Path:
    name = f'{rootname}_{planet.lower()}' + (f'_{moon.lower()}' if moon else '')
    return config.HORIZONS_DIR / f'{name}.csv'


def write_horizons_file(path: Path, positions: np.ndarray, planet: str, moon: str | None,
                        expstart: float, expend: float) -> None:
    """Save the raw query results as plain CSV (positions, not derived quantities, so
    the maths can change later without asking the server again)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    who = f'planet {planet} ({HORIZONS_IDS[planet]})' + (f', moon {moon} ({HORIZONS_IDS[moon]})' if moon else '')
    header = (
        'JPL Horizons positions seen from HST (location @hst, UTC, 1-minute steps), '
        'saved by tag_explorer\n'
        f'{who}; exposure MJD {expstart:.6f} to {expend:.6f}\n'
        'columns: jd, ra_planet_deg, dec_planet_deg, delta_planet_au, ra_moon_deg, dec_moon_deg'
    )
    np.savetxt(path, positions, delimiter=',', header=header, fmt=['%.9f'] + ['%.10f'] * 5)


def read_horizons_file(path: Path) -> np.ndarray:
    positions = np.loadtxt(path, delimiter=',', comments='#', ndmin=2)
    if positions.shape[1] != 6:
        raise ValueError(f'{path.name}: expected 6 columns, found {positions.shape[1]}')
    return positions


def saved_positions(path: Path, expstart: float, expend: float) -> np.ndarray | None:
    """The positions saved in this file if it exists, can be read, and covers MJD expstart..expend; else None
    (the reason is logged unless the file is simply not there)."""
    if not path.is_file():
        return None
    try:
        positions = read_horizons_file(path)
        mjd = positions[:, 0] - 2400000.5
        if mjd[0] > expstart or mjd[-1] < expend:
            logger.warning('%s does not cover the exposure; querying Horizons again', path.name)
            return None
        return positions
    except (OSError, ValueError, IndexError) as exc:
        logger.warning('could not use %s (%s); querying Horizons again', path.name, exc)
        return None


def load_horizons_positions(info: FrameInfo, rootname: str) -> np.ndarray:
    """Planet (and moon) positions for an exposure: from its saved file in horizons/ if
    that exists and covers the exposure, otherwise from JPL Horizons (then saved)."""
    planet = info.planet or info.target
    moon = info.target if info.planet else None
    if planet not in HORIZONS_IDS or (moon is not None and moon not in HORIZONS_IDS):
        raise ValueError(f'no ephemeris id for {planet}/{moon}')
    path = horizons_path(rootname, planet, moon)
    positions = saved_positions(path, info.expstart, info.expend)
    if positions is None:
        positions = query_horizons(planet, moon, info.expstart, info.expend)
        try:
            write_horizons_file(path, positions, planet, moon, info.expstart, info.expend)
        except OSError as exc:
            logger.warning('could not save the Horizons results: %s', exc)
    return positions


def load_planet_track(info: FrameInfo, rootname: str) -> PlanetTrack:
    """The planet's motion relative to the moon, from the saved/queried positions."""
    t, d_xi, d_eta = sky_offsets_from_positions(load_horizons_positions(info, rootname), info.expstart)
    return planet_track_from_offsets(t, d_xi, d_eta, info.cd)
