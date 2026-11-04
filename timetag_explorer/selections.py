"""Named time ranges and regions, and the bits that mark them in a selection file (see SELECTIONS.md).

A selection file (`<rootname>_tagsel.fits`) flags each event with one bit per time range in `SEL_TIME` and one bit per
region in `SEL_REGION` (16-bit columns). The two tables `TIMESEL` and `REGSEL` say what each bit means. A selection
built from others is an *expression* over their names (`oval - io`) and has no bit of its own.

Everything here is plain Python and NumPy; the FITS tables are written and read in `tagsel`. Regions are held in
**file coordinates**: hi-res pixels counted from 1, like AXIS1 and the position columns (XCORR, XREST, ...). The Explorer
draws on 0-indexed arrays, so `RegionSelection.from_array` shifts a drawn region by +1 and `.for_array()` shifts it back.
"""
from __future__ import annotations

import itertools
import re
import unicodedata
from dataclasses import dataclass, replace

import numpy as np

from .frames import FRAMES, is_map
from .regions import Box, Circle, Polygon, Region, shift_region

MAX_BITS = 16                                        # per column (the flag columns are 16-bit)
NAME_PATTERN = re.compile(r'[A-Za-z][A-Za-z0-9_]*')    # starts with a letter, so a name is never mistaken for a bit
SELECTABLE_FRAMES = ('image', 'rest', 'rest_north', 'planet', 'planet_north')    # frames with position columns
# Colors given to selections that do not name one (the Explorer's own plot colors can replace them).
DEFAULT_COLORS = ('#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#9467bd', '#8c564b', '#e377c2', '#17becf')
SHAPES = {'BOX': Box, 'CIRCLE': Circle, 'POLYGON': Polygon}


def as_ascii(text: str) -> str:
    """FITS text columns are ASCII: accents are dropped (é to e) and anything else becomes '?'."""
    decomposed = unicodedata.normalize('NFKD', str(text))
    return ''.join(c if c.isascii() and c.isprintable() else '?' for c in decomposed if not unicodedata.combining(c))


# ---- expressions ---------------------------------------------------------------------------------------------------

_TOKEN = re.compile(r'\s*(?:([A-Za-z][A-Za-z0-9_]*)|([|&\-~()]))')


def parse_expression(text: str):
    """Parse `A | B`, `A & B`, `A - B` (A and not B), `~A` and parentheses into a tree of tuples:
    ('name', 'A'), ('not', x), ('and', x, y), ('or', x, y), ('sub', x, y).
    Precedence: `~`, then `&` and `-` (left to right), then `|`. Raises ValueError for anything else; never uses eval."""
    tokens, pos = [], 0
    text = str(text)
    while pos < len(text):
        if not text[pos:].strip():
            break
        match = _TOKEN.match(text, pos)
        if match is None:
            raise ValueError(f'cannot read {text[pos:].strip()[:12]!r} in the expression {text!r}')
        tokens.append(match.group(1) or match.group(2))
        pos = match.end()
    if not tokens:
        raise ValueError('the expression is empty')
    at = [0]

    def peek():
        return tokens[at[0]] if at[0] < len(tokens) else None

    def take():
        at[0] += 1
        return tokens[at[0] - 1]

    def or_expr():
        left = and_expr()
        while peek() == '|':
            take()
            left = ('or', left, and_expr())
        return left

    def and_expr():
        left = unary()
        while peek() in ('&', '-'):
            op = take()
            left = ('and' if op == '&' else 'sub', left, unary())
        return left

    def unary():
        if peek() == '~':
            take()
            return ('not', unary())
        return atom()

    def atom():
        token = peek()
        if token == '(':
            take()
            inner = or_expr()
            if peek() != ')':
                raise ValueError(f'a ")" is missing in the expression {text!r}')
            take()
            return inner
        if token is None or not NAME_PATTERN.fullmatch(token):
            raise ValueError(f'a name was expected in the expression {text!r}' + (f', not {token!r}' if token else ''))
        return ('name', take())

    tree = or_expr()
    if at[0] != len(tokens):
        raise ValueError(f'unexpected {tokens[at[0]]!r} in the expression {text!r}')
    return tree


def expression_names(tree) -> list:
    """The names an expression tree refers to, in order of appearance."""
    if tree[0] == 'name':
        return [tree[1]]
    return [name for part in tree[1:] for name in expression_names(part)]


class _Logic:
    """How `&`, `|` and `~` combine masks (boolean arrays)."""
    @staticmethod
    def and_(a, b): return np.logical_and(a, b)
    @staticmethod
    def or_(a, b): return np.logical_or(a, b)
    @staticmethod
    def not_(a): return np.logical_not(a)


class Intervals:
    """A set of time intervals [start, stop), as sorted, separate rows. Combines with `&`, `|` and `~` exactly, so the
    exposure of a derived time range needs no events."""
    def __init__(self, rows=()):
        rows = np.asarray(rows, dtype=float).reshape(-1, 2)
        rows = rows[rows[:, 0] < rows[:, 1]]
        rows = rows[np.argsort(rows[:, 0], kind='stable')]
        merged = []
        for start, stop in rows:
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], stop)
            else:
                merged.append([start, stop])
        self.rows = np.array(merged, dtype=float).reshape(-1, 2)

    def __or__(self, other):
        return Intervals(np.vstack((self.rows, other.rows)))

    def __and__(self, other):
        lo = np.maximum(self.rows[:, None, 0], other.rows[None, :, 0])
        hi = np.minimum(self.rows[:, None, 1], other.rows[None, :, 1])
        keep = lo < hi
        return Intervals(np.column_stack((lo[keep], hi[keep])))

    def __invert__(self):
        starts = np.concatenate(([-np.inf], self.rows[:, 1]))
        stops = np.concatenate((self.rows[:, 0], [np.inf]))
        return Intervals(np.column_stack((starts, stops)))

    def contains(self, t):
        t = np.asarray(t, dtype=float)
        inside = np.zeros(t.shape, dtype=bool)
        for start, stop in self.rows:
            inside |= (t >= start) & (t < stop)
        return inside

    def exposure(self, gti) -> float:
        """Good time (s): the overlap of the good-time intervals `gti` (an (n, 2) array of START, STOP) with the set."""
        gti = np.asarray(gti, dtype=float).reshape(-1, 2)
        return float(np.sum(np.diff((self & Intervals(gti)).rows, axis=1)))


class _IntervalLogic:
    @staticmethod
    def and_(a, b): return a & b
    @staticmethod
    def or_(a, b): return a | b
    @staticmethod
    def not_(a): return ~a


# ---- one time range or region --------------------------------------------------------------------------------------

@dataclass(frozen=True, kw_only=True)
class _Selection:
    name: str = ''           # an expression needs one (other expressions refer to it); otherwise the table names it T1, R1, ...
    note: str = ''
    color: str = ''          # '#rrggbb'; filled in by the table when left empty
    visible: bool = True
    expr: str = ''           # a derived selection: an expression over earlier names; it has no bit
    bit: int = -1            # 0..15; filled in by the table for a selection that is not derived

    @property
    def derived(self) -> bool:
        return bool(self.expr)

    @property
    def mask(self) -> int:
        return 0 if self.bit < 0 else 1 << self.bit


@dataclass(frozen=True, kw_only=True)
class TimeRange(_Selection):
    """Events with start <= TIME < stop, in seconds since the exposure start. (A derived range has neither.)"""
    start: float | None = None
    stop: float | None = None

    def check(self):
        if self.derived:
            if self.start is not None or self.stop is not None:
                raise ValueError(f'time range {self.name!r} is an expression, so it has no start and stop of its own')
        elif self.start is None or self.stop is None or not (np.isfinite(self.start) and np.isfinite(self.stop)
                                                             and self.start < self.stop):
            raise ValueError(f'time range {self.name!r}: needs start < stop (seconds), got {self.start} and {self.stop}')


@dataclass(frozen=True, kw_only=True)
class RegionSelection(_Selection):
    """Events inside `region`, drawn in `frame`; coordinates in file coordinates (hi-res pixels counted from 1)."""
    region: Region | None = None
    frame: str = ''

    @classmethod
    def from_array(cls, region: Region, frame: str, *, name: str = '', **kwargs) -> 'RegionSelection':
        """A region as the Explorer draws it, on 0-indexed arrays, moved to the file's coordinates (+1)."""
        return cls(name=name, region=shift_region(region, 1.0), frame=frame, **kwargs)

    def for_array(self) -> Region:
        """The region on 0-indexed arrays (-1), for code that works on the Explorer's grids."""
        return shift_region(self.region, -1.0)

    def check(self):
        if self.derived:
            if self.region is not None or self.frame:
                raise ValueError(f'region {self.name!r} is an expression, so it has no shape or frame of its own')
            return
        if self.region is None:
            raise ValueError(f'region {self.name!r} has no shape')
        if is_map(self.frame):
            raise ValueError(f'region {self.name!r} was drawn on a planet map ({self.frame}); maps are not supported '
                             'yet, because they need an accurate planet registration')
        if self.frame not in SELECTABLE_FRAMES:
            raise ValueError(f'region {self.name!r}: unknown reference frame {self.frame!r}' if self.frame not in FRAMES
                             else f'region {self.name!r}: frame {self.frame!r} cannot be selected')


def region_params(region: Region) -> tuple[str, list]:
    """(SHAPE, numbers) of a region, as stored in the REGSEL table."""
    if isinstance(region, Box):
        return 'BOX', [region.x0, region.x1, region.y0, region.y1]
    if isinstance(region, Circle):
        return 'CIRCLE', [region.cx, region.cy, region.r]
    return 'POLYGON', [float(v) for point in region.points for v in point]


def region_from_params(shape: str, params) -> Region:
    """The region a SHAPE and its numbers (padding NaN dropped) describe."""
    values = [float(v) for v in np.asarray(params, dtype=float) if not np.isnan(v)]
    shape = str(shape).strip().upper()
    if shape == 'BOX' and len(values) == 4:
        return Box(*values)
    if shape == 'CIRCLE' and len(values) == 3:
        return Circle(*values)
    if shape == 'POLYGON' and len(values) >= 6 and len(values) % 2 == 0:
        return Polygon(tuple(zip(values[0::2], values[1::2])))
    raise ValueError(f'cannot rebuild a {shape!r} region from {len(values)} numbers')


# ---- the table of selections ---------------------------------------------------------------------------------------

class SelectionTable:
    """Time ranges (`TimeTable`) or regions (`RegionTable`) in order, each with its bit, name and expression.

    A selection needs only its numbers: a missing name becomes `T1`, `T2`, ... (`R1`, ... for regions; numbers already
    taken by a name you gave are skipped), a missing bit the lowest free one, a missing color one from a default list.
    Selections are picked by name (any case) or by bit position: `table.row('Oval')`, `table.row(0)`. The bit is the
    position, 0 to 15, not the mask value 1, 2, 4, ...; a selection built from others has no bit and is picked by name.
    """
    kind = ''
    row_class = _Selection
    prefix = ''              # the start of the names given to selections without one: T1, T2, ... or R1, R2, ...

    def __init__(self, rows=()):
        rows = list(rows)
        for row in rows:
            if not isinstance(row, self.row_class):
                raise TypeError(f'a {self.kind} table holds {self.row_class.__name__} rows, not {type(row).__name__}')
            if row.derived and not row.name:
                raise ValueError(f'a {self.kind} selection made from an expression ({row.expr!r}) needs a name')
        taken = {row.name.lower() for row in rows if row.name}
        numbers = iter(n for n in itertools.count(1) if f'{self.prefix}{n}'.lower() not in taken)
        rows = [row if row.name else replace(row, name=f'{self.prefix}{next(numbers)}') for row in rows]
        for row in rows:
            if not NAME_PATTERN.fullmatch(row.name):
                raise ValueError(f'{self.kind} name {row.name!r}: use a letter, then letters, digits or _')
            row.check()
        names = [row.name.lower() for row in rows]
        for name in set(names):
            if names.count(name) > 1:
                raise ValueError(f'two {self.kind} selections are named {name!r} (names ignore case)')
        for row in rows:
            if row.derived and row.bit >= 0:
                raise ValueError(f'{self.kind} {row.name!r} is an expression, so it has no bit')
            if row.bit >= MAX_BITS:
                raise ValueError(f'{self.kind} {row.name!r}: bits run from 0 to {MAX_BITS - 1}, not {row.bit}')
        given = [row.bit for row in rows if not row.derived and row.bit >= 0]
        if len(set(given)) != len(given):
            raise ValueError(f'two {self.kind} selections share a bit')
        taken = set(given)
        primitives = [row for row in rows if not row.derived]
        if len(primitives) > MAX_BITS:
            raise ValueError(f'at most {MAX_BITS} {self.kind} selections can have a bit (derived ones do not count); '
                             f'got {len(primitives)}')
        free = iter(b for b in range(MAX_BITS) if b not in taken)
        done = []
        for i, row in enumerate(rows):
            changes = {}
            if not row.derived and row.bit < 0:
                changes['bit'] = next(free)
            if not row.color:
                changes['color'] = DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
            if row.derived:
                tree = parse_expression(row.expr)   # (a ValueError names the problem)
                earlier = {r.name.lower() for r in rows[:i]}
                for name in expression_names(tree):
                    if name.lower() not in earlier:
                        raise ValueError(f'{self.kind} {row.name!r}: the expression {row.expr!r} refers to {name!r}, which '
                                         'is not a name above it in the table')
            done.append(replace(row, **changes) if changes else row)
        self.rows = tuple(done)

    # -- looking up
    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def __bool__(self):
        return bool(self.rows)

    @property
    def names(self) -> tuple:
        return tuple(row.name for row in self.rows)

    @property
    def primitives(self) -> tuple:
        """The selections that have a bit, in table order."""
        return tuple(row for row in self.rows if not row.derived)

    def row(self, key):
        """The selection with this name (any case) or bit position."""
        if isinstance(key, _Selection):
            key = key.name
        if isinstance(key, (int, np.integer)) and not isinstance(key, (bool, np.bool_)):
            for row in self.rows:
                if not row.derived and row.bit == key:
                    return row
            raise KeyError(f'no {self.kind} selection has bit {int(key)} (bits: '
                           f'{", ".join(str(r.bit) for r in self.primitives) or "none"})')
        if isinstance(key, str):
            for row in self.rows:
                if row.name.lower() == key.strip().lower():
                    return row
            raise KeyError(f'no {self.kind} selection is named {key!r} (names: {", ".join(self.names) or "none"})')
        raise TypeError(f'pick a {self.kind} selection by name or bit position, not {key!r}')

    # -- combining
    def evaluate(self, key, leaf, logic=_Logic):
        """`leaf(row)` for a selection with a bit; for a derived one the expression combined with `logic`."""
        row = self.row(key)
        if not row.derived:
            return leaf(row)

        def walk(tree):
            if tree[0] == 'name':
                return self.evaluate(tree[1], leaf, logic)
            if tree[0] == 'not':
                return logic.not_(walk(tree[1]))
            a, b = walk(tree[1]), walk(tree[2])
            if tree[0] == 'and':
                return logic.and_(a, b)
            if tree[0] == 'or':
                return logic.or_(a, b)
            return logic.and_(a, logic.not_(b))     # 'sub'
        return walk(parse_expression(row.expr))

    def leaves(self, key) -> tuple:
        """The selections with a bit that `key` is made of (itself, if it has one)."""
        found = []

        def collect(row):
            found.append(row)
            return True
        self.evaluate(key, collect, _Any)
        return tuple(dict.fromkeys(found))

    def mask(self, flags, key) -> np.ndarray:
        """Boolean per event: True where the event is in the selection, from a flag column (SEL_TIME or SEL_REGION)."""
        bits = np.asarray(flags).astype(np.uint16, copy=False)
        return self.evaluate(key, lambda row: (bits & np.uint16(row.mask)) != 0)

    # -- the FITS table, as plain arrays
    def to_columns(self) -> dict:
        raise NotImplementedError

    def _common_columns(self) -> dict:
        return {'NAME': [as_ascii(r.name) for r in self.rows], 'NOTE': [as_ascii(r.note) for r in self.rows],
                'BIT': np.array([r.bit for r in self.rows], dtype=np.int16),
                'MASK': np.array([r.mask for r in self.rows], dtype=np.int32)}

    def _tail_columns(self) -> dict:
        return {'EXPR': [as_ascii(r.expr) for r in self.rows], 'COLOR': [as_ascii(r.color) for r in self.rows],
                'VISIBLE': np.array([r.visible for r in self.rows], dtype=bool)}


class _Any:
    """Collects the selections an expression uses (the values are not combined)."""
    @staticmethod
    def and_(a, b): return True
    @staticmethod
    def or_(a, b): return True
    @staticmethod
    def not_(a): return True


class TimeTable(SelectionTable):
    kind = 'time'
    row_class = TimeRange
    prefix = 'T'

    def intervals(self, key) -> Intervals:
        """The set of time intervals (seconds since the exposure start) a time range covers."""
        return self.evaluate(key, lambda r: Intervals([[r.start, r.stop]]), _IntervalLogic)

    def window(self, key) -> tuple[float, float]:
        """(first start, last stop) of the intervals: the span a plot of this selection should cover."""
        rows = self.intervals(key).rows
        if not len(rows):
            raise ValueError(f'time range {self.row(key).name!r} is empty')
        return float(rows[0, 0]), float(rows[-1, 1])

    def exposure(self, key, gti) -> float:
        """Good exposure time (s) in a time range: its overlap with the good-time intervals `gti`, an (n, 2) array of
        START and STOP (a GTI table: `np.column_stack((gti['START'], gti['STOP']))`). Computed when asked, not stored."""
        return self.intervals(key).exposure(gti)

    def to_columns(self) -> dict:
        return {**self._common_columns(),
                'TSTART': np.array([np.nan if r.start is None else r.start for r in self.rows], dtype=float),
                'TSTOP': np.array([np.nan if r.stop is None else r.stop for r in self.rows], dtype=float),
                **self._tail_columns()}

    @classmethod
    def from_columns(cls, cols) -> 'TimeTable':
        rows = []
        for i in range(len(cols['NAME'])):
            expr = str(cols['EXPR'][i]).strip()
            span = {} if expr else {'start': float(cols['TSTART'][i]), 'stop': float(cols['TSTOP'][i])}   # (NaN in the file for derived)
            rows.append(TimeRange(name=str(cols['NAME'][i]).strip(), note=str(cols['NOTE'][i]).strip(),
                                  bit=int(cols['BIT'][i]), expr=expr, color=str(cols['COLOR'][i]).strip(),
                                  visible=bool(cols['VISIBLE'][i]), **span))
        return cls(rows)


class RegionTable(SelectionTable):
    kind = 'region'
    row_class = RegionSelection
    prefix = 'R'

    def frames(self, key) -> tuple:
        """The frames the selection is drawn in (one for a drawn region, several for an expression)."""
        return tuple(dict.fromkeys(row.frame for row in self.leaves(key)))

    @property
    def used_frames(self) -> tuple:
        """Every frame any region is drawn in, in table order."""
        return tuple(dict.fromkeys(row.frame for row in self.primitives))

    def to_columns(self) -> dict:
        shapes = [region_params(r.region) if r.region is not None else ('', []) for r in self.rows]
        width = max([4] + [len(p) for _, p in shapes])
        params = np.full((len(self.rows), width), np.nan)
        for i, (_, p) in enumerate(shapes):
            params[i, :len(p)] = p
        return {**self._common_columns(), 'FRAME': [as_ascii(r.frame) for r in self.rows],
                'SHAPE': [s for s, _ in shapes], 'PARAMS': params,
                'DESCR': [as_ascii(r.region.describe()) if r.region is not None else '' for r in self.rows],
                **self._tail_columns()}

    @classmethod
    def from_columns(cls, cols) -> 'RegionTable':
        n = len(cols['NAME'])
        rows = []
        for i in range(n):
            expr = str(cols['EXPR'][i]).strip()
            shape = str(cols['SHAPE'][i]).strip()
            rows.append(RegionSelection(
                name=str(cols['NAME'][i]).strip(), note=str(cols['NOTE'][i]).strip(), bit=int(cols['BIT'][i]),
                frame=str(cols['FRAME'][i]).strip(), expr=expr, color=str(cols['COLOR'][i]).strip(),
                region=None if expr else region_from_params(shape, cols['PARAMS'][i]),
                visible=bool(cols['VISIBLE'][i])))
        return cls(rows)


@dataclass(frozen=True)
class Selections:
    """The two tables read from a selection file, with the checks a person would otherwise write by hand:

        sel = read_selections('od8k57itq_tagsel.fits')
        keep = sel.in_time(t, 'orbit2') & sel.in_region(t, 'Oval') & t['EDGE_OK']     # a name (any case) or a bit position
    """
    times: TimeTable
    regions: RegionTable

    def in_time(self, t, key) -> np.ndarray:
        """Boolean per event of the events table `t`: in this time range (derived ones included)."""
        return self.times.mask(t['SEL_TIME'], key)

    def in_region(self, t, key) -> np.ndarray:
        """Boolean per event of the events table `t`: in this region (derived ones included)."""
        return self.regions.mask(t['SEL_REGION'], key)

    def time_exposure(self, key, gti) -> float:
        return self.times.exposure(key, gti)


def single_selections(time_sel, region, frame, *, array_coordinates: bool = True):
    """(TimeTable or None, RegionTable or None) for the Explorer's one time window and one region: the time range is
    named `T1` and the region `R1`. A region as drawn (`array_coordinates`) is moved to the file's coordinates."""
    times = TimeTable([TimeRange(name='T1', start=float(time_sel[0]), stop=float(time_sel[1]))]) \
        if time_sel is not None else None
    if region is None:
        return times, None
    selection = RegionSelection.from_array(region, frame, name='R1') if array_coordinates else \
        RegionSelection(name='R1', region=region, frame=frame)
    return times, RegionTable([selection])
