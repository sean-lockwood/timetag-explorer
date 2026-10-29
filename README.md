# timetag-explorer

[![Tests](https://github.com/sean-lockwood/timetag-explorer/actions/workflows/ci.yml/badge.svg)](https://github.com/sean-lockwood/timetag-explorer/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/timetag-explorer)](https://pypi.org/project/timetag-explorer/)
[![Documentation](https://app.readthedocs.org/projects/timetag-explorer/badge/?version=latest)](https://timetag-explorer.readthedocs.io)

The science code behind [TIME-TAG Explorer](https://explorer.planetarylightshow.com), a web app for exploring
HST/STIS FUV TIME-TAG observations (aurorae on Jupiter and Saturn, and other targets that move or vary).

This repository holds **only the calculations**: reading a TIME-TAG file, the processing steps, the transformations
between reference frames and planet maps, the exposure-time maps, and the generator of the standalone "selection
scripts". The web server around it (sign-in, pages, plots and mouse selection) is a separate, private project, and
this package knows nothing about it. Everything the app computes about the data comes from here, so this is the part
that scientific review should look at.

> **Status: first release, not yet validated.** [VALIDATION.md](https://github.com/sean-lockwood/timetag-explorer/blob/main/VALIDATION.md) lists each step with its reference and
> its review status. The geometric distortion correction is **disabled** in this release: `gdcorr.py` holds
> placeholder (zero) coefficients until values from a public reference file replace them.

## Install

Python 3.11 or newer.

```
pip install timetag-explorer
```

To work on the code, clone the repository and `pip install -e .` in it. Documentation:
<https://timetag-explorer.readthedocs.io>.

## A first look

```python
from timetag_explorer import load_dataset, Corrections, Box, describe_selection

ds = load_dataset('od8k57itq_tag.fits')       # a TIME-TAG file from MAST; its _spt.fits file, if next to it, is used too
corr = Corrections(randcorr=True, flatcorr=False, timecorr=True)   # which processing steps are on

edges, rate = ds.light_curve(3.0)                                    # count rate in 3 s bins
image = ds.image(time_sel=(100.0, 250.0), frame='image', corr=corr)     # counts/s per 4x4 hi-res-pixel bin
region = Box(900, 1100, 900, 1100)                                   # detector pixels
edges, region_rate = ds.light_curve(3.0, region=region, corr=corr)
print(describe_selection(ds, corr, 'image', (100.0, 250.0), region))
```

Rest frames and planet maps (`frame='rest'`, `'planet'`, `'map_rect'`, ...) need the exposure's SPT file next to the
TIME-TAG file; the planet ones also ask JPL Horizons for positions (through `astroquery`) and keep the answers in
`config.HORIZONS_DIR`. The flat field is optional: put `combined_flat.fits`
([download](https://www.planetarylightshow.com/static/data/combined_flat.fits);
[how it was made](https://www.planetarylightshow.com/aurora_processing.html#Flat-Fielding)) where `config.FLAT_FILE` points,
or set the environment variables `TAG_EXPLORER_FLAT` and `TAG_EXPLORER_HORIZONS`.

`build_selection_script(ds, time_sel, region, frame, corr)` returns a standalone Python script that reproduces a
selection on the original TIME-TAG file and writes flag columns into a copy of it.

## What is where

| Module | Contents |
|---|---|
| `regions` | `Box`, `Circle`, `Polygon`: which detector pixels a region holds |
| `frames` | the processing switches (`Corrections`), reference frames, header parsing (pointing, slew, WCS) |
| `horizons` | planet and moon positions from JPL Horizons, and the local cache of them |
| `planet` | Jupiter and Saturn orientation (IAU rotation model, light time) and the sky / surface / map projections |
| `transforms` | detector pixels to frames and maps, exposure-time maps, the distortion lookup table |
| `dataset` | `TagDataset`: events, slicing in time, binning, light curves, count-rate images |
| `analysis` | areas, mean brightness and the exposure for a selection |
| `selection_script` | the generator of selection scripts (and its version, `SELVER`) |
| `gdcorr` | the geometric distortion polynomial (placeholder coefficients for now) |
| `constants`, `config`, `log` | shared numbers; file locations; logging |

## Tests and validation

```
pip install -e .
for t in tests/test_*.py; do python "$t"; done            # quick: synthetic data, no network
python validation/check_selection_scripts.py OFAP44G5Q    # slow: runs generated scripts on a real exposure (downloads it)
```

The first set runs on every pull request. The second runs every week on GitHub Actions. See
[VALIDATION.md](https://github.com/sean-lockwood/timetag-explorer/blob/main/VALIDATION.md) for what each step is checked against.

## Contributing

Questions and concerns about the science are welcome as [issues](https://github.com/sean-lockwood/timetag-explorer/issues); proposed changes as pull requests.
Please read [CONTRIBUTING.md](https://github.com/sean-lockwood/timetag-explorer/blob/main/CONTRIBUTING.md) first. Requests for access to the web app go to the address on the
[app's front page](https://explorer.planetarylightshow.com/about), not here.

## Citing

Use the "Cite this repository" button on GitHub, or [CITATION.cff](https://github.com/sean-lockwood/timetag-explorer/blob/main/CITATION.cff). A DOI will be added once the first
release is archived.

## License

BSD 3-Clause. See [LICENSE](https://github.com/sean-lockwood/timetag-explorer/blob/main/LICENSE).
