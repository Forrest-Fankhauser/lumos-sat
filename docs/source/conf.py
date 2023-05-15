project = 'lumos-sat'
copyright = '2023, Forrest Fankhauser'
author = 'Forrest Fankhauser'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
   'sphinx.ext.autodoc',
   'sphinx.ext.autosummary',
   'sphinx.ext.viewcode',
   'sphinx_copybutton'
]

autosummary_typehints = None
templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_context = { 'display_github': False }
html_show_sourcelink = False
