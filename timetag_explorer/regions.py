"""Region shapes on the detector: hit-testing with numpy, and a one-line description for logs and FITS headers.

Coordinates are detector pixels (or, on a planet map, degrees)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np
from matplotlib.path import Path as MplPath



# --------------------------------------------------------------------------
# Region shapes (data coordinates = detector pixels)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Box:
    x0: float
    x1: float
    y0: float
    y1: float
    kind = 'box'

    def describe(self) -> str:
        return f'box x {self.x0:.0f}-{self.x1:.0f}, y {self.y0:.0f}-{self.y1:.0f}'

    def contains(self, x, y):
        return (x >= self.x0) & (x < self.x1) & (y >= self.y0) & (y < self.y1)

@dataclass(frozen=True)
class Circle:
    cx: float
    cy: float
    r: float
    kind = 'circle'

    def describe(self) -> str:
        return f'circle centre ({self.cx:.0f}, {self.cy:.0f}), r = {self.r:.0f}'

    def contains(self, x, y):
        return (x - self.cx) ** 2 + (y - self.cy) ** 2 <= self.r ** 2

@dataclass(frozen=True)
class Polygon:
    points: tuple[tuple[float, float], ...]
    kind = 'free-form'

    def describe(self) -> str:
        xs, ys = zip(*self.points)
        return (f'free-form, {len(self.points)} vertices, bounding box '
                f'x {min(xs):.0f}-{max(xs):.0f}, y {min(ys):.0f}-{max(ys):.0f}')

    def contains(self, x, y):
        out = np.zeros(x.shape, dtype=bool)
        if len(self.points) < 3:
            return out
        pts = np.asarray(self.points, dtype=float)
        (xmin, ymin), (xmax, ymax) = pts.min(axis=0), pts.max(axis=0)
        # Cheap bounding-box prefilter, then exact test on the survivors.
        cand = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
        if cand.any():
            out[cand] = MplPath(pts).contains_points(np.column_stack((x[cand], y[cand])))
        return out

Region = Union[Box, Circle, Polygon]
