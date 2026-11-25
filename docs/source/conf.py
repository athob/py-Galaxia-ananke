# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

package_name = 'galaxia_ananke'

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

from src import __metadata__ as md

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = md.__project__
copyright = md.__copyright_short__.replace('Copyright', '').strip(' ')
author = md.__author__
release = md.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_math_dollar',
    'sphinx.ext.mathjax',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'sphinx_copybutton',
    'myst_nb',
]

exclude_patterns = ['.ipynb_checkpoints/*']

napoleon_use_param = False

# Two pages can have a section of the same name: apidoc writes a "Submodules" section in the page of a package and in the one of
# each of its subpackages. A label made of the name of the page and the title of the section keeps them apart.
autosectionlabel_prefix_document = True

myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

nb_execution_timeout = 100

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'

# The theme prints it in the header, and Sphinx's own is '{project} {release} documentation': without the last word it is shorter.
html_title = f'{project} {release}'

html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": md.__url__,
            "icon": "fab fa-github-square",
            "type": "fontawesome",
        },
    ],
}

html_sidebars = {
    "**": [],
}

# -- Options for automatic API doc

autodoc_default_options = {
    'member-order': 'bysource',
    'special-members': '__init__',
    'ignore-module-all': True
}

def run_apidoc(_):
    try:
        from sphinx.ext.apidoc import main
    except ImportError:
        from sphinx.apidoc import main

    sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
    cur_dir = os.path.abspath(os.path.dirname(__file__))

    api_doc_dir = os.path.join(cur_dir, "modules")
    module = os.path.join(cur_dir, "../..", f"src/{package_name}")
    ignore = [
    ]

    main(["-M", "-f", "-e", "-T", "-d 0", "-o", api_doc_dir, module, *ignore])


def setup(app):
    app.connect("builder-inited", run_apidoc)
