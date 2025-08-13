# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------
project = 'Sanger DNA Damage Analysis Pipeline'
copyright = '2025, Allysson Allan'
author = 'Allysson Allan'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.githubpages',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# Use Furo theme for modern, minimalist design
html_theme = 'furo'

html_static_path = ['_static']

# -- Furo theme options (minimalist and professional) -----------------------
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#2563eb",  # Professional blue
        "color-brand-content": "#1e40af",
        "color-admonition-background": "#f8fafc",
        "font-stack": "'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif",
        "font-stack--monospace": "'JetBrains Mono', 'Fira Code', Consolas, 'Liberation Mono', monospace",
        "font-size--small": "0.875rem",
        "font-size--normal": "1rem",
        "font-size--large": "1.125rem",
        "line-height": "1.6",
        "sidebar-width": "280px",
        "content-padding": "3rem",
    },
    "dark_css_variables": {
        "color-brand-primary": "#60a5fa",
        "color-brand-content": "#93c5fd",
        "color-admonition-background": "#1e293b",
        "font-stack": "'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif",
        "font-stack--monospace": "'JetBrains Mono', 'Fira Code', Consolas, 'Liberation Mono', monospace",
    },
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "top_of_page_button": "edit",
    "source_repository": "https://github.com/allyssonallan/sanger_adna_damage",
    "source_branch": "main",
    "source_directory": "docs/source/",
}

# -- Extension configuration -------------------------------------------------

# Napoleon settings for clean docstring parsing
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'special-members': '__init__',
}

# Intersphinx mapping for clean cross-references
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}

# MyST settings for Markdown support
myst_enable_extensions = [
    "deflist",
    "tasklist",
    "attrs_inline",
    "colon_fence",
]

# HTML context
html_context = {
    'display_github': True,
    'github_user': 'allyssonallan',
    'github_repo': 'sanger_adna_damage',
    'github_version': 'main',
    'conf_py_path': '/docs/source/',
}

# Custom CSS for professional typography
html_css_files = [
    'https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap',
    'https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;600&display=swap',
    'custom.css',
]

# Custom JavaScript
html_js_files = [
    'custom.js',
]

# Clean up some settings
html_copy_source = False
html_show_sourcelink = False
html_show_sphinx = False
html_show_copyright = True

# Source file suffixes
source_suffix = {
    '.rst': None,
    '.md': 'myst_parser',
}

# Suppress certain warnings for cleaner builds
suppress_warnings = [
    'image.nonlocal_uri',
    'ref.doc',
]

# Clean navigation
html_sidebars = {
    "**": [
        "sidebar/scroll-start.html",
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/navigation.html",
        "sidebar/ethical-ads.html",
        "sidebar/scroll-end.html",
    ]
}
