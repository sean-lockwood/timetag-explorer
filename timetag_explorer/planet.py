"""Jupiter and Saturn: orientation against time (IAU rotation model), and the projections between the sky,
the planet's surface and the longitude / latitude maps."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np




# --------------------------------------------------------------------------
# Planet geometry (Jupiter, Saturn): orientation vs time, and the projection
# between the sky and the planet's surface
# --------------------------------------------------------------------------
#
# Conventions: longitudes are System III, increasing to the WEST (planetographic),
# latitudes are planetographic. "Screen" coordinates are what an observer outside
# the planet sees with north up: x to the right, y up, in units of the 1-bar
# equatorial radius. The sky's east is to the LEFT of that view, and the surface
# moves left to right as a prograde planet rotates, so west longitude at the
# central meridian increases with time.

AU_KM = 149597870.7
LIGHT_TIME_PER_AU = 499.00478384  # seconds
ARCSEC_PER_RAD = 206264.80624709636


@dataclass(frozen=True)
class Body:
    """IAU (WGCCRE) pole and rotation elements and radii; System III rotation."""

    name: str
    radius_eq_km: float
    radius_pol_km: float
    ra0: float        # deg; pole right ascension at J2000 and its rate [deg / Julian century]
    ra0_rate: float
    dec0: float       # deg; pole declination
    dec0_rate: float
    w0: float         # deg; prime-meridian angle at J2000 and its rate [deg / day]
    w_rate: float
    altitude_km: float  # typical height of the auroral emission above the 1-bar level


BODIES = {
    'JUPITER': Body('Jupiter', 71492.0, 66854.0, 268.056595, -0.006499, 64.495303, 0.002413,
                    284.95, 870.5366420, 400.0),
    'SATURN': Body('Saturn', 60268.0, 54364.0, 40.589, -0.036, 83.537, -0.004,
                   38.90, 810.7939024, 1100.0),
}

# TAI - UTC in seconds from the given MJD (the leap seconds since 1996).
_LEAP_MJD = np.array([50083, 50630, 51179, 53736, 54832, 56109, 57204, 57754])
_LEAP_SEC = np.array([30, 31, 32, 33, 34, 35, 36, 37])

POLAR_MAX_COLATITUDE = 50.0  # polar maps extend this far from the pole (degrees)
POLAR_TOP_LONGITUDE = 180.0  # System III longitude at the top of polar maps
MAP_MIN_MU = 0.1             # events/bins seen at a smaller cosine of the emission angle are dropped


def tdb_minus_utc(mjd):
    """TDB - UTC in seconds (leap seconds + 32.184 s; periodic TDB terms, < 2 ms, ignored)."""
    idx = np.clip(np.searchsorted(_LEAP_MJD, mjd, side='right') - 1, 0, len(_LEAP_SEC) - 1)
    return _LEAP_SEC[idx] + 32.184


@dataclass(frozen=True)
class PlanetOrientation:
    """The planet's orientation as seen from HST, sampled at t seconds after EXPSTART."""

    t: np.ndarray
    lon_obs: np.ndarray        # deg, System III west longitude of the sub-observer point (unwrapped)
    lat_obs_c: np.ndarray      # deg, planetocentric latitude of the sub-observer point
    lat_obs_g: np.ndarray      # deg, planetographic latitude of the sub-observer point
    pole_pa: np.ndarray        # deg, position angle of the north pole (east of north)
    radius_arcsec: np.ndarray  # apparent equatorial radius
    ra: np.ndarray             # deg, planet position
    dec: np.ndarray

    def at(self, times) -> dict:
        """Orientation at arbitrary times (linear interpolation)."""
        return {name: np.interp(times, self.t, getattr(self, name))
                for name in ('lon_obs', 'lat_obs_c', 'lat_obs_g', 'pole_pa', 'radius_arcsec')}


def planet_orientation(positions: np.ndarray, body: Body, expstart: float) -> PlanetOrientation:
    """Orientation from the IAU model, the planet's RA/Dec and its distance.

    positions columns: JD (UTC), planet RA, Dec (deg), distance (AU), ... (as saved in horizons/).
    The prime meridian is evaluated at the light-time-corrected epoch in TDB. The
    sub-observer point is the body-fixed direction of the vector from the planet to HST.
    """
    jd, ra, dec, delta = (positions[:, i] for i in range(4))
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    n = np.stack((np.cos(dec_r) * np.cos(ra_r), np.cos(dec_r) * np.sin(ra_r), np.sin(dec_r)), axis=1)

    days = jd + (tdb_minus_utc(jd - 2400000.5) - delta * LIGHT_TIME_PER_AU) / 86400.0 - 2451545.0
    centuries = days / 36525.0
    a0 = np.radians(body.ra0 + body.ra0_rate * centuries)
    d0 = np.radians(body.dec0 + body.dec0_rate * centuries)
    pole = np.stack((np.cos(d0) * np.cos(a0), np.cos(d0) * np.sin(a0), np.sin(d0)), axis=1)
    node = np.stack((-np.sin(a0), np.cos(a0), np.zeros_like(a0)), axis=1)  # where the equator crosses the ICRF equator
    ahead = np.cross(pole, node)
    w = np.radians(body.w0 + body.w_rate * days)[:, None]
    x_axis = node * np.cos(w) + ahead * np.sin(w)    # prime meridian
    y_axis = -node * np.sin(w) + ahead * np.cos(w)   # 90 deg east

    to_observer = -n
    lon_east = np.arctan2((to_observer * y_axis).sum(1), (to_observer * x_axis).sum(1))
    lat_c = np.arcsin(np.clip((to_observer * pole).sum(1), -1.0, 1.0))
    lat_g = np.arctan((body.radius_eq_km / body.radius_pol_km) ** 2 * np.tan(lat_c))

    east = np.stack((-np.sin(ra_r), np.cos(ra_r), np.zeros_like(ra_r)), axis=1)
    north = np.stack((-np.sin(dec_r) * np.cos(ra_r), -np.sin(dec_r) * np.sin(ra_r), np.cos(dec_r)), axis=1)
    return PlanetOrientation(
        t=(jd - 2400000.5 - expstart) * 86400.0,
        lon_obs=np.degrees(np.unwrap(-lon_east)),
        lat_obs_c=np.degrees(lat_c), lat_obs_g=np.degrees(lat_g),
        pole_pa=np.degrees(np.arctan2((pole * east).sum(1), (pole * north).sum(1))),
        radius_arcsec=body.radius_eq_km / (delta * AU_KM) * ARCSEC_PER_RAD,
        ra=ra, dec=dec,
    )


def spheroid_radii(body: Body, altitude_km: float) -> tuple[float, float]:
    """Equatorial and polar semi-axes of the (raised) spheroid, in 1-bar equatorial radii."""
    return ((body.radius_eq_km + altitude_km) / body.radius_eq_km,
            (body.radius_pol_km + altitude_km) / body.radius_eq_km)


def screen_to_surface(x, y, lat_c_obs, lon_obs, a, b):
    """Inverse orthographic projection onto a spheroid with semi-axes a, b.

    x, y: screen coordinates (right, up) in 1-bar equatorial radii. lat_c_obs and lon_obs
    (degrees) are the planetocentric latitude and west longitude of the sub-observer point;
    they may be scalars or arrays. Returns (west longitude, planetographic latitude,
    mu = cos(emission angle), on_planet).
    """
    cp, sp = np.cos(np.radians(lat_c_obs)), np.sin(np.radians(lat_c_obs))
    a2, b2 = a * a, b * b
    # Line of sight z (toward the observer): P = x*e_x + y*e_y + z*e_z lies on the spheroid.
    qa = cp * cp / a2 + sp * sp / b2
    qb = (sp * cp / b2 - sp * cp / a2) * y
    qc = (sp * sp / a2 + cp * cp / b2) * y * y + x * x / a2 - 1.0
    disc = qb * qb - qa * qc
    on = disc >= 0
    z = (-qb + np.sqrt(np.where(on, disc, 0.0))) / qa
    px, py, pz = -y * sp + z * cp, x, y * cp + z * sp
    lon = (lon_obs - np.degrees(np.arctan2(py, px))) % 360.0
    lat_g = np.degrees(np.arctan(np.tan(np.arctan2(pz, np.hypot(px, py))) * (a2 / b2)))
    nx_, ny_, nz_ = px / a2, py / a2, pz / b2
    mu = (nx_ * cp + nz_ * sp) / np.sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_)
    return lon, lat_g, mu, on


def surface_to_screen(lon_w, lat_g, lat_c_obs, lon_obs, a, b):
    """Forward projection: (west longitude, planetographic latitude) -> screen x, y and
    mu = cos(emission angle) (mu > 0: on the visible hemisphere)."""
    cp, sp = np.cos(np.radians(lat_c_obs)), np.sin(np.radians(lat_c_obs))
    a2, b2 = a * a, b * b
    phi = np.radians(lat_g)
    e2 = 1.0 - b2 / a2
    sin_phi = np.sin(phi)
    n_radius = a / np.sqrt(1.0 - e2 * sin_phi * sin_phi)  # prime-vertical radius of curvature
    rho = n_radius * np.cos(phi)
    z_body = n_radius * (1.0 - e2) * sin_phi
    d_east = np.radians(lon_obs - lon_w)                  # east longitude from the sub-observer meridian
    px, py = rho * np.cos(d_east), rho * np.sin(d_east)
    x, y = py, -px * sp + z_body * cp
    nx_, ny_, nz_ = px / a2, py / a2, z_body / b2
    mu = (nx_ * cp + nz_ * sp) / np.sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_)
    return x, y, mu


def pixels_to_screen(u, v, cx, cy, cd, pole_pa, radius_arcsec):
    """Planet-frame detector pixels -> screen coordinates (right, up) in equatorial radii.

    cd is the WCS CD matrix (deg per pixel, rows = east/north), pole_pa the position angle
    of the north pole. Sky east is to the left of the planet's own view, hence the sign."""
    dx, dy = u - cx, v - cy
    xi = (cd[0] * dx + cd[1] * dy) * 3600.0   # arcsec east
    eta = (cd[2] * dx + cd[3] * dy) * 3600.0  # arcsec north
    th = np.radians(pole_pa)
    c, s = np.cos(th), np.sin(th)
    return (-xi * c + eta * s) / radius_arcsec, (xi * s + eta * c) / radius_arcsec


def screen_to_pixels(x, y, cx, cy, cd_inv, pole_pa, radius_arcsec):
    """Inverse of pixels_to_screen."""
    th = np.radians(pole_pa)
    c, s = np.cos(th), np.sin(th)
    xi = (-x * c + y * s) * radius_arcsec / 3600.0   # degrees
    eta = (x * s + y * c) * radius_arcsec / 3600.0
    return cx + cd_inv[0][0] * xi + cd_inv[0][1] * eta, cy + cd_inv[1][0] * xi + cd_inv[1][1] * eta


def limb_outline(lat_c_obs, a, b, n=181):
    """Screen-plane outline of the spheroid (an ellipse) seen from the given sub-observer latitude."""
    psi = np.linspace(0.0, 2.0 * np.pi, n)
    sp, cp = np.sin(np.radians(lat_c_obs)), np.cos(np.radians(lat_c_obs))
    return a * np.cos(psi), np.sqrt(b * b * cp * cp + a * a * sp * sp) * np.sin(psi)


# ---- map projections (the plane the map frames are binned in) --------------------------------

def map_plane_coordinates(lon_w, lat_g, frame, lon_centre):
    """(west longitude, planetographic latitude) -> map plane (u, v, in_map).

    map_rect   u = longitude (unwrapped around lon_centre), v = latitude.
    map_north  azimuthal equidistant polar map: radius = colatitude from the north pole,
               longitude POLAR_TOP_LONGITUDE at the top, increasing clockwise.
    map_south  the same for the south pole, seen from below (mirror image).
    """
    if frame == 'map_rect':
        u = lon_centre + (lon_w - lon_centre + 180.0) % 360.0 - 180.0
        return u, lat_g, np.ones(np.shape(lat_g), dtype=bool)
    sign = 1.0 if frame == 'map_north' else -1.0
    rho = 90.0 - sign * lat_g
    ang = np.radians(lon_w - POLAR_TOP_LONGITUDE)
    return sign * rho * np.sin(ang), rho * np.cos(ang), (rho >= 0) & (rho <= POLAR_MAX_COLATITUDE)


def map_plane_to_surface(u, v, frame):
    """Inverse of map_plane_coordinates (west longitude, planetographic latitude)."""
    if frame == 'map_rect':
        return u % 360.0, v
    sign = 1.0 if frame == 'map_north' else -1.0
    rho = np.hypot(u, v)
    return (POLAR_TOP_LONGITUDE + np.degrees(np.arctan2(sign * u, v))) % 360.0, sign * (90.0 - rho)


def map_bin_area(u, v, frame, binsize, a_km, b_km):
    """Surface area (km^2) of the map bins centred at (u, v) on a spheroid with semi-axes a_km, b_km."""
    _, lat = map_plane_to_surface(u, v, frame)
    phi = np.radians(np.clip(lat, -90.0, 90.0))
    e2 = 1.0 - (b_km / a_km) ** 2
    w = 1.0 - e2 * np.sin(phi) ** 2
    n_rad, m_rad = a_km / np.sqrt(w), a_km * (1.0 - e2) / w ** 1.5
    per_rad2 = m_rad * n_rad * np.cos(phi)             # km^2 per radian^2 of (longitude x latitude)
    deg2 = np.radians(1.0) ** 2
    if frame == 'map_rect':
        return per_rad2 * deg2 * binsize ** 2
    # Polar plane: x, y = rho*(sin, cos)(angle), rho in degrees, angle in radians, so
    # dx dy = rho d(rho) d(angle) and d(lon) d(lat) = (pi/180) dx dy / rho   (cos(phi)/rho is finite at the pole).
    rho = np.maximum(np.hypot(u, v), 1e-6)
    return per_rad2 * np.radians(1.0) * binsize ** 2 / rho
