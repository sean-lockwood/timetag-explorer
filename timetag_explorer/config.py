"""Where the package finds its files, and the address it writes into the scripts it generates.

The defaults come from environment variables (or, failing that, from the current folder). A program that uses the
package can also assign to these names, which are read each time they are needed:

    from timetag_explorer import config
    config.FLAT_FILE = Path('/data/combined_flat.fits')
    config.HORIZONS_DIR = Path('/data/horizons')
"""
from __future__ import annotations

import os
from pathlib import Path


def _env_path(name: str, default: str) -> Path:
    return Path(os.environ.get(name) or default)


# The flat field: a 2048 x 2048 image of multiplicative weights (optional: without it, counts are unweighted).
FLAT_FILE = _env_path('TAG_EXPLORER_FLAT', 'combined_flat.fits')
# Where the JPL Horizons ephemerides are saved and reused (one small CSV per exposure and body).
HORIZONS_DIR = _env_path('TAG_EXPLORER_HORIZONS', 'horizons')
# Web address of the tool that is written into generated scripts and FITS headers (empty: the default address).
PUBLIC_URL = os.environ.get('TAG_EXPLORER_URL', '').strip().rstrip('/')
