"""Planet geometry against answers that can be derived independently: round trips, areas of spheres and caps,
the rotation rate, and the leap seconds."""
import support  # noqa: F401
import numpy as np

from timetag_explorer.planet import (AU_KM, BODIES, LIGHT_TIME_PER_AU, map_bin_area, map_plane_coordinates,
                                     map_plane_to_surface, limb_outline, planet_orientation, screen_to_surface,
                                     spheroid_radii, surface_to_screen, tdb_minus_utc)

jup = BODIES['JUPITER']
assert (jup.radius_eq_km, jup.radius_pol_km) == (71492.0, 66854.0)          # IAU 2015 radii
assert spheroid_radii(jup, 0.0) == (1.0, 66854.0 / 71492.0)

# Leap seconds: TAI-UTC was 36 s until the end of 2016 and 37 s from 1 January 2017 (MJD 57754); TDB = TT + ~0.
assert tdb_minus_utc(57753) == 36 + 32.184 and tdb_minus_utc(57754) == 37 + 32.184
assert tdb_minus_utc(50082) == 29 + 32.184 or tdb_minus_utc(50082) == 30 + 32.184   # before the table: first entry

# The sub-observer point is at the centre of the disc, looking straight at it.
x, y, mu = surface_to_screen(100.0, 0.0, 0.0, 100.0, 1.0, 1.0)
assert abs(x) < 1e-12 and abs(y) < 1e-12 and abs(mu - 1.0) < 1e-12
# ... and the limb of a sphere seen from its equator is a circle of radius 1.
lx, ly = limb_outline(0.0, 1.0, 1.0)
assert np.allclose(np.hypot(lx, ly), 1.0)

# Round trip, surface -> screen -> surface, over the visible hemisphere of an oblate planet seen from 3 N.
a, b = spheroid_radii(jup, 400.0)
rng = np.random.default_rng(0)
lon = rng.uniform(0, 360, 2000)
lat = rng.uniform(-85, 85, 2000)
sx, sy, smu = surface_to_screen(lon, lat, 3.0, 120.0, a, b)
vis = smu > 0.05
lon2, lat2, mu2, on = screen_to_surface(sx[vis], sy[vis], 3.0, 120.0, a, b)
assert on.all() and np.allclose(mu2, smu[vis], atol=1e-9)
assert np.allclose(((lon2 - lon[vis] + 180) % 360) - 180, 0, atol=1e-7) and np.allclose(lat2, lat[vis], atol=1e-7)
# a point outside the disc is not on the planet
assert not screen_to_surface(np.array([2.0]), np.array([0.0]), 3.0, 120.0, a, b)[3][0]

# Map planes: round trips, and the area of a whole sphere and of a polar cap.
for frame in ('map_rect', 'map_north', 'map_south'):
    lo = rng.uniform(0, 360, 200)
    la = rng.uniform(45, 85, 200) * (-1 if frame == 'map_south' else 1) if frame != 'map_rect' else rng.uniform(-80, 80, 200)
    u, v, ok = map_plane_coordinates(lo, la, frame, 180.0)
    lo2, la2 = map_plane_to_surface(u, v, frame)
    assert np.allclose(((lo2 - lo + 180) % 360) - 180, 0, atol=1e-7) and np.allclose(la2, la, atol=1e-7), frame
R = 71492.0
step = 0.5
cu = np.arange(0, 360, step) + step / 2
cv = np.arange(-90, 90, step) + step / 2
gu, gv = np.meshgrid(cu, cv)
area = map_bin_area(gu, gv, 'map_rect', step, R, R).sum()
assert abs(area / (4 * np.pi * R ** 2) - 1) < 1e-4, area
# north polar map out to 50 degrees of colatitude: a cap of area 2 pi R^2 (1 - cos 50 deg)
step = 0.25
c = np.arange(-50, 50, step) + step / 2
gu, gv = np.meshgrid(c, c)
inside = np.hypot(gu, gv) < 50.0
area = map_bin_area(gu[inside], gv[inside], 'map_north', step, R, R).sum()
cap = 2 * np.pi * R ** 2 * (1 - np.cos(np.radians(50.0)))
assert abs(area / cap - 1) < 2e-3, (area, cap)

# Orientation: Jupiter 4.5 AU away is 71492 km / 4.5 AU across, and System III turns 870.5366 deg per day.
n = 5
jd = 2460000.5 + np.arange(n) / 24.0
positions = np.column_stack([jd, np.full(n, 150.0), np.full(n, 15.0), np.full(n, 4.5)])
orient = planet_orientation(positions, jup, jd[0] - 2400000.5)
assert np.allclose(orient.radius_arcsec, 71492.0 / (4.5 * AU_KM) * 206264.80624709636)
rate = np.diff(orient.lon_obs) / np.diff(orient.t)                  # deg per second, west longitude
assert np.allclose(np.abs(rate) * 86400, 870.5366420, rtol=1e-3), rate * 86400
assert abs(LIGHT_TIME_PER_AU - 499.0048) < 1e-3
assert np.all(np.abs(orient.lat_obs_c) < 90) and np.all(np.abs(orient.lat_obs_g) < 90)
print('planet ok')
