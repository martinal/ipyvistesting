jupyter-volren-widget
===============================

A Jupyter Widget for Volume Rendering with WebGL

Installation
------------

To install use pip:

    $ pip install volrenwidget
    $ jupyter nbextension enable --py --sys-prefix volrenwidget


For a development installation (requires npm),

    $ git clone https://github.com/martinal/jupyter-volren-widget.git
    $ cd jupyter-volren-widget
    $ pip install -e .
    $ jupyter nbextension install --py --symlink --sys-prefix volrenwidget
    $ jupyter nbextension enable --py --sys-prefix volrenwidget
