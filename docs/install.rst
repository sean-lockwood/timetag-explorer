Installation
============

Python 3.11 or newer.

.. code-block:: console

   pip install timetag-explorer

To work on the code, clone the repository and install it in editable mode:

.. code-block:: console

   git clone https://github.com/sean-lockwood/timetag-explorer.git
   cd timetag-explorer
   pip install -e .

Running the tests
-----------------

.. code-block:: console

   for t in tests/test_*.py; do python "$t"; done            # quick: synthetic data, no network
   python validation/check_selection_scripts.py OFAP44G5Q    # slow: runs generated scripts on a real exposure

Files the package uses
----------------------

* The exposure's TIME-TAG file (``*_tag.fits``) from MAST, and, for the rest frames and planet maps, its SPT file
  (``*_spt.fits``) next to it.
* Optionally the flat field, ``combined_flat.fits``
  (`download <https://www.planetarylightshow.com/static/data/combined_flat.fits>`_;
  `how it was made <https://www.planetarylightshow.com/aurora_processing.html#Flat-Fielding>`_). Put it where
  :data:`timetag_explorer.config.FLAT_FILE` points, or set the environment variable ``TAG_EXPLORER_FLAT``.
* A folder where JPL Horizons positions are kept (``TAG_EXPLORER_HORIZONS`` or
  :data:`timetag_explorer.config.HORIZONS_DIR`). The planet frames ask Horizons for positions through ``astroquery``
  the first time and reuse the saved answers afterwards.
