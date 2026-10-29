"""TIME-TAG Explorer science code: the calculations behind the web app for HST/STIS FUV TIME-TAG data.

This package has no web-server code. It reads a TIME-TAG (and SPT) file, applies the processing steps, transforms
positions between reference frames and planet maps, and writes the standalone selection scripts.

    from timetag_explorer import TagDataset, Corrections, Box
"""
from __future__ import annotations

__version__ = '1.0.1'

from .analysis import describe_selection
from .constants import SELECTION_SCRIPT_VERSION
from .dataset import TagDataset, load_dataset
from .frames import DEFAULT_CORRECTIONS, FRAMES, Corrections, FrameInfo
from .regions import Box, Circle, Polygon, Region
from .selection_script import SelectionScript, build_selection_script

__all__ = ['__version__', 'SELECTION_SCRIPT_VERSION', 'TagDataset', 'load_dataset', 'describe_selection',
           'Corrections', 'DEFAULT_CORRECTIONS', 'FRAMES', 'FrameInfo', 'Box', 'Circle', 'Polygon', 'Region',
           'SelectionScript', 'build_selection_script']
