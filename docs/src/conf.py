# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information

project = 'microcosm'
copyright = '2026, Minh-Quan Ton-Ngoc, Si-Nguyen Mai'

release = '1.0'
version = '1.0.1'

# -- General configuration
html_context = {
    "doc_author": "Si-Nguyen Mai"
}

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    "sphinx.ext.extlinks",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx_tabs.tabs",
    "sphinx_copybutton",
    "sphinx_design",
    "myst_parser",
    "jupyter_sphinx",
    "sphinx_togglebutton",
    "nbsphinx",
    "numpydoc",
    # "sphinx_sitemap",
    "sphinxcontrib.mermaid",
    "sphinxcontrib.video",
    "sphinxcontrib.youtube",
    "sphinx_click",
    "sphinx_contributors",
    # "sphinx_iconify",
    "click_extra.sphinx",
    "sphinx_datatables",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
    "numpy": ("https://numpy.org/devdocs/", None),
}
intersphinx_disabled_domains = ['std']
todo_include_todos = True
myst_enable_extensions = ["colon_fence"]
# jupyter_sphinx_thebelab_config = {
#     "requestKernel": True,
# }
jupyter_sphinx_require_url = ""
# iconify_script_url = ""
# nbsphinx_requirejs_path = ""
# sitemap_excludes = []

# The master toctree document
# master_doc = 'index'

templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = [
    'css/custom.css',
]
html_extra_path = ["_public"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_title = 'microcosm'
html_theme = 'shibuya'
# html_baseurl = ""
# sitemap_url_scheme = "{link}"

html_copy_source = False
html_show_sourcelink = False

# html_additional_pages = {
#     "molepi": "molepi.html",
# }

html_logo = "_static/metaflow_logo.png"
# html_favicon = "_static/icon.png"

html_theme_options = {
    "accent_color": "brown", 
    "color_mode": "light", 
    "logo_target": "/",
    # "light-logo": "_static/logo-light.svg",
    # "dark_logo": "_static/logo-dark.svg",
    # "og_image_url": "https://shibuya.lepture.com/icon.png",
    # "twitter_creator": "lepture",
    # "twitter_site": "lepture",
    # "discussion_url": "https://github.com/lepture/shibuya/discussions",
    # "twitter_url": "https://twitter.com/lepture",
    "github_url": "https://github.com/tnmquann/metaflow",
    "globaltoc_expand_depth": 1,
    "open_in_perplexity": True,
    "nav_links": [
        {
            "title": "Quick Start",
            "url": "/quickstart",
        },
        {
            "title": "Installation",
            "url": "/installation/",
        },
        # {
        #     "title": "Usage",
        #     "url": "/usage/",
        # },
        # {
        #     "title": "Tutorial",
        #     "url": "/tutorials/",
        # },
        # {
        #     "title": "Contributors",
        #     "url": "/contributors",
        # },
    ],
}

html_sidebars = {
  "**": [
    "localtoc.html",
  ]
}
# -- Options for PDF output

