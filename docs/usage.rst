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
original TIME-TAG file and writes a copy of it, ``<rootname>_tagsel.fits``, with columns added: the selection flags
(``SEL_TIME``, ``SEL_REGION``, ``EDGE_OK``), the corrected image-frame position (``XCORR``, ``YCORR``) and the position
in each reference frame (``XREST``/``YREST``, ``XRESTN``/``YRESTN``, ``XPLAN``/``YPLAN``, ``XPLANN``/``YPLANN``), and the
weights (``FLATW``, ``VIGW``, ``WEIGHT``); ``AXIS1`` and ``AXIS2`` are left as recorded. A selection is optional: without
one the corrections are still applied and the selection columns are simply absent. The primary header says which steps
were performed (``RANDCORR``, ``GEOCORR``, ``FLATCORR``, ``VIGNCORR``, ``EXPTCORR``: ``PERFORM`` or ``OMIT``) and whether
a time or space selection was made (``TIMESEL``, ``SPCSEL``). The script is short: it holds
the choices and calls :func:`timetag_explorer.tagsel.make_tagsel`, and it names the exact version of this package that
the Explorer ran (and the ``pip install`` line for it). The region is tested at each event's own position, not on the
Explorer's 4-pixel bins. The version of the script generator,
:data:`timetag_explorer.SELECTION_SCRIPT_VERSION`, is written into every script and into the ``SELVER`` keyword of the
FITS files they write; it is raised whenever a generated script would compute or write something different.

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
``tagsel``                 what a selection script runs: flags, corrected positions and weights in ``_tagsel.fits``
``gdcorr``                 the geometric distortion polynomial (placeholder coefficients for now)
``vignette``               the experimental vignetting model (F25SRF2 only; off by default)
``constants``, ``config``  shared numbers; file locations; logging
=========================  ==================================================================================
