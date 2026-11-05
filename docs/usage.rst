A first look
============

.. code-block:: python

   from timetag_explorer import load_dataset, Corrections, Box, describe_selection

   ds = load_dataset('od8k57itq_tag.fits')       # a TIME-TAG file from MAST; its _spt.fits file, if next to it, is used too
   corr = Corrections(randcorr=True, flatcorr=False, timecorr=True)   # which processing steps are on

   edges, rate = ds.light_curve(3.0)                                    # count rate in 3 s bins
   image = ds.image(time_sel=(100.0, 250.0), frame='image', corr=corr)     # counts/s per 4x4 hi-res-pixel bin
   region = Box(900, 1100, 900, 1100)                                   # detector pixels
   edges, region_rate = ds.light_curve(3.0, region=region, corr=corr)
   print(describe_selection(ds, corr, 'image', (100.0, 250.0), region))

Rest frames and planet maps (``frame='rest'``, ``'planet'``, ``'map_rect'``, ...) need the exposure's SPT file next to
the TIME-TAG file; the planet ones also ask JPL Horizons for positions.

Selection scripts
-----------------

:func:`timetag_explorer.build_selection_script` returns a standalone Python script that reproduces a selection on the
original TIME-TAG file. The script writes a copy of the file, ``<rootname>_tagsel.fits``, with columns added. The
original tables and ``AXIS1`` and ``AXIS2`` are left as recorded. A selection is optional: without one, the
corrections are still applied and the selection columns are simply absent.

The columns added are:

* the selection flags: ``SEL_TIME`` and ``SEL_REGION`` (one bit per named time range or region) and ``EDGE_OK``;
* the corrected image-frame position, ``XCORR`` and ``YCORR``;
* the position in each reference frame, ``XREST``/``YREST`` and ``XPLAN``/``YPLAN``. The north-up versions,
  ``XRESTN``/``YRESTN`` and ``XPLANN``/``YPLANN``, make the file much larger and are written only if the script sets
  ``INCLUDE_NORTH_UP = True``;
* the weights, ``FLATW`` and ``VIGW``. An event counts as their product; there is no combined column.

The primary header records what was done, in one block under the label ``/ TIMETAG-EXPLORER KEYWORDS``: which steps
ran (``TTERAND``, ``TTEGEO``, ``TTEFLAT``, ``TTEVIGN``, ``TTEEDGE``: ``DONE`` or ``OMIT``), the flat-field file used
(``TTEFFILE``), and whether a time or space selection was made (``TTETSEL``, ``TTERSEL``). These names are deliberately
not the calibration pipeline's step keywords (``RANDCORR``, ``FLATCORR`` and so on, whose ``PERFORM`` means something
else to ``calstis`` and ``stistools``), which are left as they were.

The script is short. It holds the choices and calls :func:`timetag_explorer.tagsel.make_tagsel`, and it names the exact
version of this package that the Explorer ran, with the ``pip install`` line for it. The region is tested at each
event's own position, not on the Explorer's 4-pixel bins.

The version of the script generator, :data:`timetag_explorer.SELECTION_SCRIPT_VERSION`, is written into every script
and into the ``SELVER`` keyword of the FITS files they write. It is raised whenever a generated script would compute or
write something different.

Several named selections
~~~~~~~~~~~~~~~~~~~~~~~~

:func:`timetag_explorer.build_selections_script` takes any number of named time ranges and regions
(:class:`~timetag_explorer.selections.TimeRange`, :class:`~timetag_explorer.selections.RegionSelection`; the format is in ``SELECTIONS.md``).
The ``TIMESEL`` and ``REGSEL`` extensions of the output name the bits, and :func:`timetag_explorer.read_selections` picks
events by name (any case) or bit position, including selections written as expressions of others::

   sel = read_selections('od8k57itq_tagsel.fits')
   t = Table.read('od8k57itq_tagsel.fits', hdu=1)
   keep = sel.in_time(t, 'orbit2') & sel.in_region(t, 'Oval') & t['EDGE_OK']

Without a selection script
~~~~~~~~~~~~~~~~~~~~~~~~~~

The script is only a convenience: the same file can be made by calling the package. A selection needs just its numbers;
a missing name becomes ``T1``, ``T2``, ... (``R1``, ``R2``, ... for regions), a missing bit the lowest free one, a missing
color one from a default list, and it is visible. Name them, or give a ``bit``, when you want to find them again in a
way that does not depend on their order::

   from timetag_explorer import Corrections, Box
   from timetag_explorer.selections import TimeRange, RegionSelection
   from timetag_explorer.tagsel import make_tagsel

   path = make_tagsel('od8k57itq',
                      corr=Corrections(randcorr=True, flatcorr=False, timecorr=True),
                      times=[TimeRange(start=0, stop=1500), TimeRange(name='orbit2', start=3000, stop=4500)],
                      regions=[RegionSelection(region=Box(900, 1100, 900, 1100), frame='image')])    # hi-res pixels from 1

The first range is ``T1`` (bit 0), the second ``orbit2`` (bit 1) and the region ``R1`` (bit 0 of ``SEL_REGION``).
A selection made from others (``RegionSelection(name='both', expr='R1 | R2')``) must have a name. A region drawn on
0-indexed arrays, as the Explorer draws them, goes through ``RegionSelection.from_array(region, frame)``, which moves it
to the file's coordinates.

Regions are in the file's coordinates, hi-res pixels counted from 1 like ``AXIS1``. Area and exposure are not stored:
:func:`timetag_explorer.analysis.region_stats` and :func:`timetag_explorer.analysis.time_exposure` work them out when asked.

What is where
-------------

=========================  ==================================================================================
Module                     Contents
=========================  ==================================================================================
``regions``                ``Box``, ``Circle``, ``Polygon``: which detector pixels a region holds
``frames``                 the processing switches (``Corrections``), reference frames, header parsing
``horizons``               planet and moon positions from JPL Horizons, and the local cache of them
``planet``                 Jupiter and Saturn orientation and the sky / surface / map projections
``transforms``             detector pixels to frames and maps, exposure-time maps, the distortion lookup table
``dataset``                ``TagDataset``: events, slicing in time, binning, light curves, count-rate images
``analysis``               areas, mean brightness and the exposure for a selection
``selection_script``       the generator of selection scripts (and its version, ``SELVER``)
``selections``             named time ranges and regions: bits, names, expressions, the ``TIMESEL`` / ``REGSEL`` rows
``tagsel``                 what a selection script runs: flags, corrected positions, weights and selection tables in ``_tagsel.fits``
``gdcorr``                 the geometric distortion polynomial (placeholder coefficients for now)
``vignette``               the experimental vignetting model (F25SRF2 only; off by default)
``constants``, ``config``  shared numbers; file locations; logging
=========================  ==================================================================================
