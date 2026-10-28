# Contributing

Thank you for looking at the science. Reports of mistakes are the most useful contribution of all.

## What belongs here

This repository is the science code only: calculations on TIME-TAG events and their geometry. The web server that hosts
the app is not part of it, and changes to it cannot be made here. The package must stay free of web-server code; a test
checks that.

## Raising a concern

Open an issue with the "Science question or concern" form. It helps to say which step is involved (see
[VALIDATION.md](VALIDATION.md)), what you expected and why (an equation, a handbook, a reference file), and the
exposure you saw it on (the rootname is enough).

## Proposing a change

1. Fork, branch, and make the change with a test that would have failed before it. Small synthetic datasets are
   enough for most things (see `tests/support.py`).
2. Run the tests: `for t in tests/test_*.py; do python "$t"; done`.
3. If the change can alter a number the code produces, say so in `CHANGELOG.md` in plain words, raise the minor
   version (`timetag_explorer/__init__.py` and `CITATION.cff`, which a test keeps in step), and cite the source of
   the new equation or value in the pull request.
4. If a generated selection script would compute or write something different, raise `SELECTION_SCRIPT_VERSION` in
   `timetag_explorer/constants.py` too.
5. Open the pull request. Changes to results are reviewed with more care than anything else, so expect questions.

The maintainer merges after review. Your commits stay in the history under your name, and you are credited in the
release notes.

By contributing you agree that your contribution is released under the project's BSD 3-Clause license.

## Releases

Each release raises `__version__`, updates `CHANGELOG.md` and `CITATION.cff`, and is tagged `vX.Y.Z`. Publishing the
GitHub release archives it on Zenodo (once connected), which is where the DOI comes from. The web app installs a
specific tag, so a release is only used by the app after the maintainer moves the app to it.
