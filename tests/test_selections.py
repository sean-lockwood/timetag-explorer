"""Named time ranges and regions: names and bits, expressions, interval algebra, the table columns, coordinates, and the
exposure and area that are worked out when asked. No FITS file is needed (the files are checked on real data by
validation/check_selection_scripts.py)."""
import support  # noqa: F401
import types

import numpy as np
from support import OFF, dataset

from timetag_explorer import analysis
from timetag_explorer.regions import Box, Circle, Polygon, shift_region
from timetag_explorer.selections import (Intervals, RegionSelection, RegionTable, TimeRange, TimeTable, as_ascii,
                                         expression_names, parse_expression, region_from_params, region_params,
                                         single_selections, Selections)


def refused(function, text):
    try:
        function()
    except (ValueError, KeyError, TypeError) as err:
        assert text in str(err), (text, str(err))
        return
    raise SystemExit(f'expected an error containing {text!r}')


# ---- expressions: precedence ~, then & and - (left to right), then ^, then | --------------------------------------------------
assert parse_expression('A') == ('name', 'A')
assert parse_expression('a | b & c') == ('or', ('name', 'a'), ('and', ('name', 'b'), ('name', 'c')))
assert parse_expression('a & b - c') == ('sub', ('and', ('name', 'a'), ('name', 'b')), ('name', 'c'))
assert parse_expression('a - b - c') == ('sub', ('sub', ('name', 'a'), ('name', 'b')), ('name', 'c'))
assert parse_expression('~a & b') == ('and', ('not', ('name', 'a')), ('name', 'b'))
assert parse_expression('~(a | b)') == ('not', ('or', ('name', 'a'), ('name', 'b')))
assert expression_names(parse_expression('(a | b) - ~c')) == ['a', 'b', 'c']
assert parse_expression('a | b ^ c & d') == ('or', ('name', 'a'), ('xor', ('name', 'b'), ('and', ('name', 'c'), ('name', 'd'))))
assert parse_expression('a ^ b ^ c') == ('xor', ('xor', ('name', 'a'), ('name', 'b')), ('name', 'c'))
assert parse_expression('a - b ^ c') == ('xor', ('sub', ('name', 'a'), ('name', 'b')), ('name', 'c'))
for bad, why in (('', 'empty'), ('a b', 'unexpected'), ('a |', 'name was expected'), ('(a', 'missing'), ('a $ b', 'cannot read'),
                 ('1 | a', 'cannot read'), ('__import__("os")', 'cannot read')):
    refused(lambda: parse_expression(bad), why)

_xor = TimeTable([TimeRange(name='a', start=0.0, stop=60.0), TimeRange(name='b', start=40.0, stop=100.0),
                  TimeRange(name='either', expr='a ^ b')])
assert _xor.intervals('either').rows.tolist() == [[0.0, 40.0], [60.0, 100.0]]

# ---- time ranges: bits, names, picking ------------------------------------------------------------------------------
times = TimeTable([TimeRange(name='orbit1', start=0.0, stop=1500.0, note='first'),
                   TimeRange(name='Orbit2', start=3000.0, stop=4500.0, visible=False),
                   TimeRange(name='both', expr='orbit1 | ORBIT2'),
                   TimeRange(name='gap', expr='~both')])
assert [r.bit for r in times] == [0, 1, -1, -1] and [r.mask for r in times] == [1, 2, 0, 0]
assert times.names == ('orbit1', 'Orbit2', 'both', 'gap') and len(times) == 4 and times.primitives == times.rows[:2]
assert all(r.color.startswith('#') and len(r.color) == 7 for r in times) and len({r.color for r in times}) == 4
assert times.row('ORBIT2').bit == 1 and times.row(1).name == 'Orbit2' and times.row(np.int64(0)).name == 'orbit1'
assert times.row('  Both ').derived and times.row(times.rows[0]).name == 'orbit1'
refused(lambda: times.row('nope'), 'no time selection is named')
refused(lambda: times.row(5), 'has bit 5')
refused(lambda: times.row(2), 'has bit 2')            # a derived selection has no bit
refused(lambda: times.row(1.5), 'by name or bit')
refused(lambda: times.row(True), 'by name or bit')
# the bit is the position, not the mask value: 2 is not "bit 1"
flags = np.array([0, 1, 2, 3, 2 ** 15], dtype=np.uint16)
assert times.mask(flags, 0).tolist() == [False, True, False, True, False]
assert times.mask(flags, 'orbit2').tolist() == [False, False, True, True, False]
assert times.mask(flags, 'both').tolist() == [False, True, True, True, False]
assert times.mask(flags, 'gap').tolist() == [True, False, False, False, True]
assert times.mask(flags.astype(np.int16), 'orbit2').tolist() == times.mask(flags, 'orbit2').tolist()     # signed readers

# interval algebra: the exposure of a derived range needs no events
assert times.intervals('both').rows.tolist() == [[0.0, 1500.0], [3000.0, 4500.0]]
assert times.intervals('gap').rows.tolist() == [[-np.inf, 0.0], [1500.0, 3000.0], [4500.0, np.inf]]
gti = [[0.0, 1000.0], [1400.0, 3100.0], [4000.0, 5000.0]]
assert times.exposure('orbit1', gti) == 1100.0 and times.exposure('Orbit2', gti) == 100.0 + 500.0
assert times.exposure('both', gti) == 1700.0 and times.window('both') == (0.0, 4500.0)
assert times.exposure('gap', gti) == 1500.0 + 500.0           # [1500, 3000) inside the second GTI row, [4500, 5000) in the third
gen = np.random.default_rng(3)
t_ev = gen.uniform(-100, 6000, 20000)
bits = np.zeros(len(t_ev), dtype=np.uint16)
for row in times.primitives:
    bits[(t_ev >= row.start) & (t_ev < row.stop)] |= np.uint16(row.mask)
for name in ('orbit1', 'both', 'gap'):          # event flags and interval membership agree, derived ones included
    assert np.array_equal(times.mask(bits, name), times.intervals(name).contains(t_ev)), name
overlap = TimeTable([TimeRange(name='a', start=0, stop=10), TimeRange(name='b', start=5, stop=20), TimeRange(name='c', expr='a - b')])
assert overlap.intervals('c').rows.tolist() == [[0.0, 5.0]] and overlap.exposure('c', [[0, 100]]) == 5.0
assert Intervals([[0, 5], [5, 8], [20, 21]]).rows.tolist() == [[0.0, 8.0], [20.0, 21.0]]        # touching intervals join
assert Intervals([[3, 3]]).exposure([[0, 9]]) == 0.0 and Intervals().exposure([[0, 9]]) == 0.0
assert Intervals([[0, 4]]).contains(np.array([0.0, 3.999, 4.0])).tolist() == [True, True, False]      # half-open

# ---- the rules a table keeps ------------------------------------------------------------------------------------------
refused(lambda: TimeTable([TimeRange(name='a', start=1, stop=1)]), 'start < stop')
refused(lambda: TimeTable([TimeRange(name='a', start=float('nan'), stop=1)]), 'start < stop')
refused(lambda: TimeTable([TimeRange(name='a')]), 'start < stop')
refused(lambda: TimeTable([TimeRange(name='1a', start=0, stop=1)]), 'use a letter')
refused(lambda: TimeTable([TimeRange(name='a b', start=0, stop=1)]), 'use a letter')
refused(lambda: TimeTable([TimeRange(name='a', start=0, stop=1), TimeRange(name='A', start=2, stop=3)]), 'ignore case')
refused(lambda: TimeTable([TimeRange(name='a', start=0, stop=1), TimeRange(name='c', expr='a | b')]), 'not a name above it')
refused(lambda: TimeTable([TimeRange(name='c', expr='a'), TimeRange(name='a', start=0, stop=1)]), 'not a name above it')
refused(lambda: TimeTable([TimeRange(name='a', start=0, stop=1, bit=16)]), 'bits run from 0 to 15')
refused(lambda: TimeTable([TimeRange(name='a', start=0, stop=1, bit=2), TimeRange(name='b', start=2, stop=3, bit=2)]), 'share a bit')
refused(lambda: TimeTable([TimeRange(name='a', start=0, stop=1), TimeRange(name='b', expr='a', start=0, stop=1)]), 'no start and stop')
refused(lambda: TimeTable([RegionSelection(name='a', region=Box(0, 1, 0, 1), frame='image')]), 'holds TimeRange')
sixteen = TimeTable([TimeRange(name=f't{i}', start=i, stop=i + 0.5) for i in range(16)] + [TimeRange(name='any', expr=' | '.join(f't{i}' for i in range(16)))])
assert [r.bit for r in sixteen.primitives] == list(range(16)) and len(sixteen) == 17            # derived ones do not use bits
refused(lambda: TimeTable([TimeRange(name=f't{i}', start=i, stop=i + 0.5) for i in range(17)]), 'at most 16')
assert sixteen.mask(np.array([2 ** 15], dtype=np.uint16), 15).tolist() == [True]              # bit 15 is the top bit
# an explicit bit is kept, and the others fill the gaps
gaps = TimeTable([TimeRange(name='a', start=0, stop=1, bit=3), TimeRange(name='b', start=2, stop=3), TimeRange(name='c', start=4, stop=5)])
assert [r.bit for r in gaps] == [3, 0, 1]
# a table needs only the numbers: names, bits, colors and visibility are filled in
bare = TimeTable([TimeRange(start=0, stop=1), TimeRange(start=2, stop=3), TimeRange(start=4, stop=5)])
assert bare.names == ('T1', 'T2', 'T3') and [r.bit for r in bare] == [0, 1, 2] and all(r.visible for r in bare)
assert len({r.color for r in bare}) == 3 and all(len(r.color) == 7 and r.color[0] == '#' for r in bare)
assert bare.row('t2').bit == 1 and bare.row(2).name == 'T3'
# names already taken (any case) are skipped; unnamed rows are numbered in order
skip = TimeTable([TimeRange(start=0, stop=1), TimeRange(name='t2', start=2, stop=3), TimeRange(start=4, stop=5)])
assert skip.names == ('T1', 't2', 'T3')
assert TimeTable([TimeRange(start=0, stop=1), TimeRange(name='all', expr='t1 | t1')]).names == ('T1', 'all')
refused(lambda: TimeTable([TimeRange(start=0, stop=1), TimeRange(expr='t1')]), 'needs a name')
assert RegionTable([RegionSelection(region=Box(1, 2, 1, 2), frame='image'),
                    RegionSelection.from_array(Box(0, 1, 0, 1), 'rest')]).names == ('R1', 'R2')

# ---- regions: coordinates, frames, tables ---------------------------------------------------------------------------
drawn = Box(10.0, 20.0, 30.0, 40.0)                           # as the Explorer draws it: 0-indexed arrays
sel = RegionSelection.from_array(drawn, 'rest', name='oval')
assert sel.region == Box(11.0, 21.0, 31.0, 41.0) and sel.for_array() == drawn                # file coordinates are +1
assert shift_region(Circle(5.0, 6.0, 2.0), 1.0) == Circle(6.0, 7.0, 2.0)
assert shift_region(Polygon(((0.0, 0.0), (4.0, 0.0), (0.0, 4.0))), 1.0, 2.0).points == ((1.0, 2.0), (5.0, 2.0), (1.0, 6.0))
# a pixel's event is in the region in both ways of counting
x0, y0 = np.array([10.0, 19.9, 20.0]), np.array([30.0, 39.9, 35.0])
assert np.array_equal(sel.region.contains(x0 + 1, y0 + 1), sel.for_array().contains(x0, y0))
regions = RegionTable([sel, RegionSelection(name='io', region=Circle(100.0, 100.0, 5.0), frame='image', note='Io, east'),
                       RegionSelection(name='outside', expr='oval - io')])
assert [r.bit for r in regions] == [0, 1, -1] and regions.used_frames == ('rest', 'image')
assert regions.frames('oval') == ('rest',) and regions.frames('outside') == ('rest', 'image')
assert [r.name for r in regions.leaves('outside')] == ['oval', 'io'] and regions.leaves(1)[0].name == 'io'
flags = np.array([0, 1, 2, 3], dtype=np.uint16)
assert regions.mask(flags, 'outside').tolist() == [False, True, False, False]                # oval and not io
# a planet-map region is allowed in a table and in a file (its coordinates go in the MAP columns)
maprow = RegionTable([RegionSelection(name='m', region=Box(0, 1, 0, 1), frame='map_rect'), RegionSelection(name='i', region=Box(0, 1, 0, 1), frame='image')])
assert maprow.used_frames == ('map_rect', 'image') and not hasattr(maprow, 'require_file_frames')
refused(lambda: RegionTable([RegionSelection(name='m', region=Box(0, 1, 0, 1), frame='nowhere')]), 'unknown reference frame')
# degrees are not counted from 1: a map region is not shifted, a pixel region is
degrees = RegionSelection.from_array(Box(10.0, 20.0, -5.0, 5.0), 'map_rect', name='d')
assert degrees.region == Box(10.0, 20.0, -5.0, 5.0) and degrees.for_array() == Box(10.0, 20.0, -5.0, 5.0)
pix = RegionSelection.from_array(Box(10.0, 20.0, 30.0, 40.0), 'image', name='p')
assert pix.region == Box(11.0, 21.0, 31.0, 41.0) and pix.for_array() == Box(10.0, 20.0, 30.0, 40.0)
refused(lambda: RegionTable([RegionSelection(name='m', region=Box(0, 1, 0, 1), frame='sky')]), 'unknown reference frame')
refused(lambda: RegionTable([RegionSelection(name='m', frame='image')]), 'has no shape')
refused(lambda: RegionTable([RegionSelection(name='m', region=Box(0, 1, 0, 1))]), 'unknown reference frame')
refused(lambda: RegionTable([RegionSelection(name='a', region=Box(0, 1, 0, 1), frame='image'),
                             RegionSelection(name='b', expr='a', frame='image')]), 'no shape or frame')

# shapes round-trip through the numbers the table stores
for region in (Box(1.5, 2.5, 3.5, 4.5), Circle(1.0, 2.0, 3.0), Polygon(((0.0, 0.0), (4.0, 0.0), (2.0, 3.0), (-1.0, 5.0)))):
    shape, params = region_params(region)
    assert region_from_params(shape, params + [float('nan')] * 3) == region                    # (NaN padding is dropped)
refused(lambda: region_from_params('BOX', [1, 2, 3]), 'cannot rebuild')
refused(lambda: region_from_params('STAR', [1, 2, 3]), 'cannot rebuild')

cols = regions.to_columns()
assert list(cols) == ['NAME', 'NOTE', 'BIT', 'MASK', 'FRAME', 'SHAPE', 'PARAMS', 'DESCR', 'EXPR', 'COLOR', 'VISIBLE', 'CURVE']
assert cols['BIT'].dtype == np.int16 and cols['BIT'].tolist() == [0, 1, -1] and cols['MASK'].tolist() == [1, 2, 0]
assert cols['SHAPE'] == ['BOX', 'CIRCLE', ''] and cols['FRAME'] == ['rest', 'image', ''] and cols['EXPR'] == ['', '', 'oval - io']
assert cols['PARAMS'].shape == (3, 4) and cols['PARAMS'][0].tolist() == [11.0, 21.0, 31.0, 41.0] and np.isnan(cols['PARAMS'][1, 3])
back = RegionTable.from_columns(cols)
assert back.rows == regions.rows
poly = RegionTable([RegionSelection(name='p', frame='image', region=Polygon(tuple((float(i), float(i * i)) for i in range(1, 8))))])
assert poly.to_columns()['PARAMS'].shape == (1, 14) and RegionTable.from_columns(poly.to_columns()).rows == poly.rows
tcols = times.to_columns()
assert list(tcols) == ['NAME', 'NOTE', 'BIT', 'MASK', 'TSTART', 'TSTOP', 'EXPR', 'COLOR', 'VISIBLE']
assert TimeTable.from_columns(tcols).rows == times.rows and tcols['VISIBLE'].tolist() == [True, False, True, True]
assert np.isnan(tcols['TSTART'][2]) and 'EXPTIME' not in tcols and 'AREA' not in cols            # worked out when asked
# FITS text is ASCII
assert as_ascii('Io, east µ') == 'Io, east ?'      # (the micro sign decomposes to a Greek letter, which becomes ?)
assert as_ascii('café') == 'cafe' and as_ascii('line\nbreak') == 'line?break'
assert RegionTable([RegionSelection(name='x', region=Box(0, 1, 0, 1), frame='image', note='naïve')]).to_columns()['NOTE'] == ['naive']

# the Explorer's single window and region
t1, r1 = single_selections((10.0, 20.0), drawn, 'rest')
assert t1.names == ('T1',) and r1.names == ('R1',) and r1.rows[0].region == Box(11.0, 21.0, 31.0, 41.0)
assert single_selections(None, None, 'image') == (None, None) and single_selections((1.0, 2.0), None, 'image')[1] is None

# reading by name or bit from an events table
sel_all = Selections(times, regions)
events = {'SEL_TIME': np.array([1, 2, 0], dtype=np.uint16), 'SEL_REGION': np.array([0, 1, 3], dtype=np.uint16)}
assert sel_all.in_time(events, 'Orbit1').tolist() == [True, False, False] and sel_all.in_time(events, 1).tolist() == [False, True, False]
assert sel_all.in_region(events, 'OVAL').tolist() == [False, True, True] and sel_all.in_region(events, 'outside').tolist() == [False, True, False]
assert sel_all.time_exposure('orbit1', gti) == 1100.0

# ---- exposure and area, when asked -------------------------------------------------------------------------------
span = np.linspace(1, 2048, 400).round().astype(int)         # events across the detector, so the Explorer's bin grid covers it
ds = dataset(np.linspace(0, 99, 400), span, span[::-1], gti=((0.0, 40.0), (50.0, 99.0)))
whole = RegionTable([RegionSelection.from_array(Box(0.0, 2048.0, 0.0, 2048.0), 'image', name='big'),
                     RegionSelection.from_array(Box(400.0, 600.0, 400.0, 800.0), 'image', name='mid'),
                     RegionSelection.from_array(Circle(1000.0, 1000.0, 100.0), 'image', name='ring'),
                     RegionSelection(name='rest', frame='rest', region=Box(1.0, 2.0, 1.0, 2.0)),
                     RegionSelection(name='mixed', expr='mid | rest'),
                     RegionSelection(name='midonly', expr='mid & big')])
tt = TimeTable([TimeRange(name='early', start=10.0, stop=30.0), TimeRange(name='late', start=45.0, stop=95.0),
                TimeRange(name='either', expr='early | late')])
stats = analysis.region_stats(ds, whole, 'mid', OFF)
assert stats.frame == 'image' and stats.area == 200 * 400                                        # hi-res px^2, on 4x4 bins
assert stats.mean_exposure == 89.0 and stats.min_exposure == 89.0 and stats.max_exposure == 89.0   # all the good time
stats = analysis.region_stats(ds, whole, 'mid', OFF, tt, 'early')
assert stats.mean_exposure == 20.0 and stats.area == 200 * 400
stats = analysis.region_stats(ds, whole, 'mid', OFF, tt, 'either')
assert stats.mean_exposure == 20.0 + 45.0                                                        # [10, 30) and [50, 95) are good
assert analysis.region_stats(ds, whole, 'midonly', OFF).area == 200 * 400                      # a derived region: its bins
assert analysis.region_stats(ds, whole, 'big', OFF).area == 2048 * 2048
ring = analysis.region_stats(ds, whole, 'ring', OFF)
assert abs(ring.area - np.pi * 100 ** 2) < 0.06 * np.pi * 100 ** 2                              # bins approximate the disc
refused(lambda: analysis.region_stats(ds, whole, 'mixed', OFF), 'more than one frame')
# the exposure-time correction on: the Explorer's per-bin exposure (the same field of view here, so the same time)
on = OFF.__class__(randcorr=False, geocorr=False, flatcorr=False, timecorr=True)
stats = analysis.region_stats(ds, whole, 'mid', on, tt, 'early')
assert abs(stats.mean_exposure - 20.0) < 1e-6 and stats.min_exposure <= stats.mean_exposure <= stats.max_exposure
assert analysis.time_exposure(ds, tt, 'early') == 20.0 and analysis.time_exposure(ds, tt, 'late') == 45.0
assert analysis.time_exposure(ds, tt, 'either') == 65.0
print('selections ok')
