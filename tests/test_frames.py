"""Header parsing and the small rules around frames and light-curve bins."""
import support  # noqa: F401

from timetag_explorer.constants import DEFAULT_TIME_BIN, MAX_LIGHT_CURVE_BINS, MIN_TIME_BIN
from timetag_explorer.frames import (FRAMES, FrameInfo, clamp_binsize, is_map, is_north, is_planet, names_from_levels,
                                     parse_mt_levels)

# Motion-tracking levels: parts of one level are joined in order
header = {'MT_LV1_1': 'TYPE=STD,NAME=JUPI', 'MT_LV1_2': 'TER', 'MT_LV2_1': 'TYPE=POS_ANGLE,ANG=44.0,R_RAD=0.5',
          'OTHER': 'x'}
levels = parse_mt_levels(header)
assert levels[0] == {'TYPE': 'STD', 'NAME': 'JUPITER'}
assert levels[1] == {'TYPE': 'POS_ANGLE', 'ANG': '44.0', 'R_RAD': '0.5'}
assert parse_mt_levels({}) == []
assert isinstance(names_from_levels(levels), tuple)

# Light-curve bins are always usable, whatever was typed
assert clamp_binsize(3) == 3.0
assert clamp_binsize('nan') == DEFAULT_TIME_BIN and clamp_binsize(None) == DEFAULT_TIME_BIN
assert clamp_binsize(0) == MIN_TIME_BIN and clamp_binsize(-5) == MIN_TIME_BIN
assert clamp_binsize(1e-9, tmax=1000.0) == 1000.0 / MAX_LIGHT_CURVE_BINS      # not billions of bins

# Frame names
assert is_map('map_rect') and not is_map('planet') and is_planet('planet_north') and is_north('rest_north')
assert 'image' in FRAMES and all(isinstance(v, str) for v in FRAMES.values())
info = FrameInfo(orientat=0.0, slew_angle=0.0, slew_rate=0.0)
assert info.unavailable_reason('image') is None and info.unavailable_reason('map_rect') is not None
print('frames ok')
