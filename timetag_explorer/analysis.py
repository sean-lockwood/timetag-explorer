"""Numbers and text that describe a selection: the area and mean brightness of a region (on the detector or on a planet
map), the exposure, and the frame and planet geometry in use."""
from __future__ import annotations

import numpy as np

from .constants import SPATIAL_BIN, TIME_CORR_MIN_FRACTION, TimeSel
from .dataset import TagDataset
from .frames import ARCSEC_PER_LOWRES_PIXEL, Corrections, frame_label, is_map, is_planet
from .planet import map_bin_area
from .regions import Region
from .transforms import counts_label


def _map_region_stats(ds: TagDataset, region: Region, frame: str, corr: Corrections, time_sel):
    """(surface area in km^2, mean brightness) of a region drawn on a planet map.

    Mean brightness is the counts over the pixel-seconds that saw the region's bins
    (counts/s per 4x4 hi-res pixels), so it does not depend on how much of the region was visible."""
    grid = ds.grid(frame, corr)
    inside = ds.region_grid(region, frame, corr)
    geom = ds.planet_geometry()
    alt = geom.altitude_km(corr)
    cu = grid.xlo + (np.arange(grid.nx) + 0.5) * grid.binsize
    cv = grid.ylo + (np.arange(grid.ny) + 0.5) * grid.binsize
    gu, gv = np.meshgrid(cu, cv)
    area = float(map_bin_area(gu[inside], gv[inside], frame, grid.binsize,
                              geom.body.radius_eq_km + alt, geom.body.radius_pol_km + alt).sum())
    img = ds.image(time_sel, frame, corr)
    pixel_seconds, full = ds.map_exposure(time_sel, frame, corr)
    den = pixel_seconds if corr.timecorr else full
    weight = np.where(den >= TIME_CORR_MIN_FRACTION * full, den, 0.0)[inside]
    mean = float((img[inside] * weight).sum() / weight.sum()) if weight.sum() > 0 else float('nan')
    return area, mean


def describe_selection(ds: TagDataset, corr: Corrections, frame: str, time_sel: TimeSel | None, region: Region | None,
                       show_limb: bool = False) -> str:
    """Text summary of the frame, the planet, the time selection and the region."""
    label = frame_label(ds.frame_info, frame)
    lines = []
    if frame != 'image':
        info = ds.frame_info
        if info.slew_rate:
            shift = abs(info.slew_rate) / 86400.0 / ARCSEC_PER_LOWRES_PIXEL * 2.0 * ds.tmax
            detail = (f'slew of {info.slew_rate:g} arcsec/day along {info.slew_angle:.1f} deg '
                      f'removed (about {shift:.0f} hi-res px over the exposure)')
        else:
            detail = 'this observation has no time-dependent slew'
        if is_planet(frame) or is_map(frame) and info.planet:
            detail += (f'; {(info.planet or "planet").title()} motion relative to '
                       f'{(info.target or "the moon").title()} removed')
        lines.append(f'Frame: {label}; {detail}')
    if ds.frame_info.body and (is_map(frame) or show_limb):
        geom = ds.planet_geometry()
        og = {k: float(v) for k, v in geom.orient.at(0.5 * ds.tmax).items()}
        cx, cy = geom.centre(corr)
        lines.append(
            f'{geom.body.name} at mid-exposure: sub-observer longitude {og["lon_obs"] % 360:.1f} W, '
            f'latitude {og["lat_obs_g"]:.1f}, pole PA {og["pole_pa"]:.1f} deg, radius {og["radius_arcsec"]:.2f}"; '
            f'centre ({cx:.1f}, {cy:.1f}) px = header guess ({geom.centre0[0]:.1f}, {geom.centre0[1]:.1f}) '
            f'+ offset ({corr.centre_dx:g}, {corr.centre_dy:g}); emission altitude {geom.altitude_km(corr):.0f} km')

    total = ds.exposure()
    if time_sel is None:
        lines.append(f'All times: {total:.1f} s of good exposure')
    else:
        t0, t1 = time_sel
        lines.append(f'Times {t0:.1f}-{t1:.1f} s: {ds.exposure(t0, t1):.1f} s of good '
                     f'exposure (of {total:.1f} s total)')

    if region is None:
        lines.append('Region: none selected')
    elif is_map(frame):
        area, mean_all = _map_region_stats(ds, region, frame, corr, None)
        means = [f'all times {mean_all:.4g}']
        if time_sel is not None:
            means.append(f'selected times {_map_region_stats(ds, region, frame, corr, time_sel)[1]:.4g}')
        lines.append(f'Region: {region.describe()} (map degrees); surface area {area:,.0f} km\u00b2; '
                     f'mean brightness {", ".join(means)} ({counts_label(ds, corr)}/s per '
                     f'{SPATIAL_BIN}x{SPATIAL_BIN} hi-res px)')
    else:
        grid = ds.region_grid(region, frame, corr)
        area = int(grid.sum()) * SPATIAL_BIN ** 2
        # Summing the rate image over the region's bins gives the region's total rate.
        rates = [f'all times {ds.image(None, frame, corr)[grid].sum():.4g}']
        if time_sel is not None:
            rates.append(f'selected times {ds.image(time_sel, frame, corr)[grid].sum():.4g}')
        lines.append(f'Region: {region.describe()} (hi-res px, {label.lower()}); '
                     f'area {area:,} hi-res px\u00b2; total rate {", ".join(rates)} '
                     f'({counts_label(ds, corr)}/s)')
    return '\n'.join(lines)

