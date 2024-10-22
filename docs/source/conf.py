# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# import m2r
sys.path.insert(0, os.path.abspath("../.."))
import LoopStructural
from pyvista.plotting.utilities.sphinx_gallery import DynamicScraper
from pyvista.core.utilities.docs import linkcode_resolve  # noqa: F401
from pyvista.core.utilities.docs import pv_html_page_context

##Pyvista stuff

# Manage errors
import pyvista
from pathlib import Path

pyvista.set_error_output_file('errors.txt')
# Ensure that offscreen rendering is used for docs generation
pyvista.OFF_SCREEN = True  # Not necessary - simply an insurance policy
# Preferred plotting style for documentation
pyvista.set_plot_theme('document')
pyvista.global_theme.window_size = [1024, 768]
pyvista.global_theme.font.size = 22
pyvista.global_theme.font.label_size = 22
pyvista.global_theme.font.title_size = 22
pyvista.global_theme.return_cpos = False
pyvista.set_jupyter_backend(None)
# Save figures in specified directory
pyvista.FIGURE_PATH = str(Path('./images/').resolve() / 'auto-generated/')
if not Path(pyvista.FIGURE_PATH).exists():
    Path(pyvista.FIGURE_PATH).mkdir()

# necessary when building the sphinx gallery
pyvista.BUILDING_GALLERY = True
os.environ['PYVISTA_BUILDING_GALLERY'] = 'true'

###

# -- Project information -----------------------------------------------------

project = "LoopStructural"
copyright = "2020, Lachlan Grose"
author = "Lachlan Grose"

# The full version, including alpha/beta/rc tags
release = LoopStructural.__version__

# -- Internationalization ----------------------------------------------------

# specifying the natural language populates some key tags
language = "en"

# -- General configuration ---------------------------------------------------
autoclass_content = "both"  # include both class docstring and __init__
autodoc_default_flags = [
    # Make sure that any autodoc declarations show the right members
    "members",
    "inherited-members",
    "private-members",
    "show-inheritance",
]
autosummary_generate = True  # Make _autosummary files and include them
napoleon_numpy_docstring = True  # False  # Force consistency, leave only Google
napoleon_use_rtype = False  # More legible
autosummary_imported_members = True
autosummary_ignore_module_all = False
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # Need the autodoc and autosummary packages to generate our docs.
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    # The Napoleon extension allows for nicer argument formatting.
    "sphinx.ext.napoleon",
    # add sphinx gallery
    "sphinx_gallery.gen_gallery",
    # citations
    "myst_parser",
    "sphinxcontrib.bibtex",
    # 'enum_tools.autoenum',
    'jupyter_sphinx',
    # 'notfound.extension',
    'numpydoc',
    'pyvista.ext.coverage',
    'pyvista.ext.plot_directive',
    'pyvista.ext.viewer_directive',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.linkcode',  # This adds the button ``[Source]`` to each Python API site by calling ``linkcode_resolve``
    'sphinx.ext.extlinks',
    'sphinx.ext.intersphinx',
    # 'sphinx_copybutton',
    'sphinx_design',  # tabs for gallery
    'sphinx_gallery.gen_gallery',
    # 'sphinxcontrib.asciinema',
    # 'sphinx_tags',
    # 'sphinx_toolbox.more_autodoc.overloads',
    # 'sphinx_toolbox.more_autodoc.typevars',
    # 'sphinx_toolbox.more_autodoc.autonamedtuple',
    # 'sphinxext.opengraph',
    # 'sphinx_sitemap',
]
bibtex_bibfiles = ["docs_references.bib"]
bibtex_default_style = "plain"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

html_theme_options = {
    "icon_links": [
        {
            "name": "Twitter",
            "url": "https://twitter.com/loop3d",
            "icon": "fab fa-twitter",
        },
        {
            "name": "GitHub",
            "url": "https://github.com/loop3d/LoopStructural",
            "icon": "fab fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/LoopStructural",
            "icon": "fa-custom fa-pypi",
        },
    ],
    #     "navbar_start": ["navbar-logo", "navbar-version"],
    #     "use_edit_page_button": True,
    "collapse_navigation": True,
    "external_links": [
        {"name": "Loop3d", "url": "https://www.loop3d.org"},
    ],
    "header_links_before_dropdown": 4,
    "logo": {
        "text": "LoopStructural - {}".format(release),
        "image_light": "_static/infinity_loop_icon.svg",
        "image_dark": "_static/infinity_loop_icon.svg",
    },
}
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"
html_sourcelink_suffix = ""
html_last_updated_fmt = ""  # to reveal the build date in the pages meta

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_js_files = ["custom-icon.js"]
todo_include_todos = True

autosummary_mock_imports = [
    "LoopStructural.interpolators._cython",
]
# Sphinx gallery examples
# from LoopStructural.visualisation.sphinx_scraper import Scraper as LoopScraper
from sphinx_gallery.sorting import FileNameSortKey
import pyvista


pyvista.BUILDING_GALLERY = True


class ResetPyVista:
    """Reset pyvista module to default settings."""

    def __call__(self, gallery_conf, fname):
        """Reset pyvista module to default settings.

        If default documentation settings are modified in any example, reset here.
        """
        import pyvista

        pyvista._wrappers['vtkPolyData'] = pyvista.PolyData
        pyvista.set_plot_theme('document')

    def __repr__(self):
        return 'ResetPyVista'


reset_pyvista = ResetPyVista()
sphinx_gallery_conf = {
    'pypandoc': True,
    "examples_dirs": ["../../examples/"],
    "gallery_dirs": ["_auto_examples/"],  # path to where to save gallery generated output
    "within_subsection_order": FileNameSortKey,
    "reference_url": {"LoopStructural": None},
    'doc_module': 'pyvista',
    'image_scrapers': (DynamicScraper(), 'matplotlib'),
    'first_notebook_cell': '%matplotlib inline',
    'reset_modules': (reset_pyvista,),
    'reset_modules_order': 'both',
}

suppress_warnings = ['config.cache']

import re

# -- .. pyvista-plot:: directive ----------------------------------------------
from numpydoc.docscrape_sphinx import SphinxDocString

IMPORT_PYVISTA_RE = r'\b(import +pyvista|from +pyvista +import)\b'
IMPORT_MATPLOTLIB_RE = r'\b(import +matplotlib|from +matplotlib +import)\b'

plot_setup = """
from pyvista import set_plot_theme as __s_p_t
__s_p_t('document')
del __s_p_t
"""
plot_cleanup = plot_setup


def _str_examples(self):
    examples_str = '\n'.join(self['Examples'])

    if (
        self.use_plots
        and re.search(IMPORT_MATPLOTLIB_RE, examples_str)
        and 'plot::' not in examples_str
    ):
        out = []
        out += self._str_header('Examples')
        out += ['.. plot::', '']
        out += self._str_indent(self['Examples'])
        out += ['']
        return out
    elif re.search(IMPORT_PYVISTA_RE, examples_str) and 'plot-pyvista::' not in examples_str:
        out = []
        out += self._str_header('Examples')
        out += ['.. pyvista-plot::', '']
        out += self._str_indent(self['Examples'])
        out += ['']
        return out
    else:
        return self._str_section('Examples')


SphinxDocString._str_examples = _str_examples


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
import sphinx_book_theme  # noqa: F401

html_theme = 'sphinx_book_theme'
html_context = {
    'github_user': 'pyvista',
    'github_repo': 'pyvista',
    'github_version': 'main',
    'doc_path': 'doc/source',
    'examples_path': 'examples',
}
html_show_sourcelink = False
html_copy_source = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = False


def get_version_match(semver):
    """Evaluate the version match for the multi-documentation."""
    if semver.endswith('dev0'):
        return 'dev'
    major, minor, _ = semver.split('.')
    return f'{major}.{minor}'


# def setup(app):
#     app.add_stylesheet('custom.css')
def setup(app):  # noqa: D103
    app.connect('html-page-context', pv_html_page_context)
    app.add_css_file('copybutton.css')
    app.add_css_file('no_search_highlight.css')
