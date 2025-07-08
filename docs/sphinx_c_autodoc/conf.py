# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
# sys.path.insert(0, os.path.abspath('.'))
#sys.path.append(os.path.abspath('../../source'))
#cautodoc_root = os.path.abspath('../../source')

# -- Project information -----------------------------------------------------

project = 'Solve-MOND-MM'
#copyright = '2020, d'
author = 'Claudio Llinares'

# The full version, including alpha/beta/rc tags
release = 'e'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
#"sphinx.ext.autodoc", 
#"sphinx.ext.napoleon", 
extensions = ["sphinx.ext.viewcode", 
              "sphinx.ext.napoleon",
              "sphinx.ext.autodoc", 
              "sphinx_c_autodoc",
              "sphinx_c_autodoc.napoleon"] 

#              "sphinx_c_autodoc.viewcode"]
#c_autodoc_roots = ['../../source/']
# This is a local copy of only the files we want to document:
c_autodoc_roots = ['./sphinx-c-apidoc/source/']
autodoc_default_options = {
    "members": True, 
    "private-members": True
}
primary_domain = "c"

# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
#exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
exclude_patterns = [ 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'classic'
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']
