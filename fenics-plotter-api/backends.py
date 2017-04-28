# -*- coding: utf-8 -*-


from __future__ import print_function
from __future__ import division

import numpy as np


class PlotterBackend:
    """Base class for plotter backends, taking only numpy arrays."""
    def __init__(self):
        pass


    # TODO: Implement scatter plots through the entire
    # pipeline with some backend and Jupyter notebooks,
    # to test the design along the entire stack.
    def scatter(self, points, point_values, **kwargs):
        """Scatter plot.

        @param points
        @param point_values
        """
        raise NotImplementedError("Missing implementation in plotting backend.")


    def wireframe(self, vertices, edges, vertex_values=None, edge_values=None, **kwargs):
        """Wireframe plot.

        @param vertices
        @param edges
        @param vertex_values
        @param edge_values
        """
        raise NotImplementedError("Missing implementation in plotting backend.")


    def surface(self, vertices, faces, vertex_values=None, face_values=None, **kwargs):
        """Surface plot.

        @param vertices
        @param faces
        @param vertex_values
        @param face_values
        """
        raise NotImplementedError("Missing implementation in plotting backend.")


    def displacement(self, vertices, cells, vertex_values, **kwargs):
        """Displacement plot.

        @param vertices
        @param cells
        @param vertex_values
        """
        raise NotImplementedError("Missing implementation in plotting backend.")


# FIXME: Implement some backends here or in separate files,
# just add the functions as they are implemented and let the
# missing functions fail with NotImplementedError in the base class


class MatplotlibBackend(PlotterBackend):
    pass  # FIXME: Convert some code to implement the API for 2D data here


class BqplotBackend(PlotterBackend):
    pass  # FIXME: Convert some code to implement the API for 2D data here


class ScivijsBackend(PlotterBackend):
    pass  # FIXME: Convert some code to implement the API for 2D data here


class K3dBackend(PlotterBackend):
    pass  # FIXME: Convert some code to implement the API for 2D data here


class VTKBackend(PlotterBackend):
    pass  # FIXME: Convert some code to implement the API for 2D data here


class VTKBackend(PlotterBackend):
    pass  # FIXME: Convert some code to implement the API for 2D data here

