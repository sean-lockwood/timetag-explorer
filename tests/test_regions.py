"""Region shapes: which pixels they contain, and how they describe themselves."""
import support  # noqa: F401
import numpy as np

from timetag_explorer.regions import Box, Circle, Polygon

x, y = np.meshgrid(np.arange(0, 10, dtype=float), np.arange(0, 10, dtype=float))

# Box: half-open, so adjacent boxes never count a pixel twice
box = Box(2, 5, 3, 6)
m = box.contains(x, y)
assert m.sum() == 9 and m[3:6, 2:5].all()
assert box.contains(np.array([2.0]), np.array([3.0]))[0] and not box.contains(np.array([5.0]), np.array([3.0]))[0]
assert Box(2, 5, 3, 6).describe() == 'box x 2-5, y 3-6'

# Circle: the boundary counts as inside
c = Circle(5, 5, 2)
assert c.contains(np.array([7.0]), np.array([5.0]))[0] and not c.contains(np.array([7.1]), np.array([5.0]))[0]
assert abs(c.contains(x, y).sum() - np.pi * 4) <= 3          # a pixel grid approximates pi r^2
assert Circle(5, 5, 2).describe().startswith('circle centre (5, 5), r = 2')

# Polygon: a square given as vertices, and a triangle
sq = Polygon(((1, 1), (4, 1), (4, 4), (1, 4)))
inside = sq.contains(np.array([2.5, 0.5, 5.0]), np.array([2.5, 2.5, 2.5]))
assert list(inside) == [True, False, False]
tri = Polygon(((0, 0), (8, 0), (0, 8)))
assert tri.contains(np.array([1.0]), np.array([1.0]))[0] and not tri.contains(np.array([7.0]), np.array([7.0]))[0]
assert not Polygon(((0, 0), (1, 1))).contains(x, y).any()    # fewer than three vertices contain nothing
assert 'free-form, 4 vertices' in sq.describe()

# hashable (they are used in cache keys) and immutable
assert len({Box(0, 1, 0, 1), Box(0, 1, 0, 1), Circle(0, 0, 1)}) == 2
try:
    box.x0 = 9
    raise SystemExit('a Box must be immutable')
except AttributeError:
    pass
print('regions ok')
