TIME-TAG Explorer |release|
===========================

The science code behind `TIME-TAG Explorer <https://explorer.planetarylightshow.com>`_, a web app for exploring
HST/STIS FUV TIME-TAG observations (aurorae on Jupiter and Saturn, and other targets that move or vary).

This package holds **only the calculations**: reading a TIME-TAG file, the processing steps, the transformations
between reference frames and planet maps, the exposure-time maps, and the generator of the standalone "selection
scripts". The web server around it is a separate, private project, and this package knows nothing about it. Everything
the app computes about the data comes from here, so this is the part that scientific review should look at.

.. warning::

   First releases, not yet validated. The
   `validation table <https://github.com/sean-lockwood/timetag-explorer/blob/main/VALIDATION.md>`_ lists each step with
   its reference and review status. The geometric distortion correction uses the STIS IDC table ``1771408ro_idc.fits``
   (bundled with the package; see :mod:`timetag_explorer.gdcorr`), adapted to the 2048 x 2048 hi-res pixel grid.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   usage
   api
   changelog

Source code, issues and pull requests: https://github.com/sean-lockwood/timetag-explorer
