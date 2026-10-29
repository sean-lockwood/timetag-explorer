"""Sphinx configuration for the timetag-explorer documentation (built by Read the Docs)."""
import timetag_explorer

project = 'timetag-explorer'
author = 'Sean Lockwood'
copyright = '2026, Sean Lockwood'
release = timetag_explorer.__version__
version = '.'.join(release.split('.')[:2])

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',      # the docstrings use INPUTS:/OUTPUTS: style sections in places
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',      # "[source]" links: the code is what reviewers read
]
autodoc_member_order = 'bysource'
autodoc_default_options = {'members': True, 'show-inheritance': True}
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
}

html_theme = 'sphinx_rtd_theme'
html_title = f'timetag-explorer {release}'
exclude_patterns = ['_build']
