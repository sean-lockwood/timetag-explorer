Changelog
=========

The changelog is kept in the repository, where each entry says whether a computed number can have changed:
https://github.com/sean-lockwood/timetag-explorer/blob/main/CHANGELOG.md

Version numbers
---------------

* **patch** (1.0.x): no result changes (fixes to wording, tests, documentation, performance).
* **minor** (1.x.0): results can change (a corrected equation, a new step, new reference data), or features are added.
* **major** (x.0.0): the Python interface changes incompatibly.

``SELECTION_SCRIPT_VERSION`` is raised separately, whenever a generated script would compute or write something
different.
