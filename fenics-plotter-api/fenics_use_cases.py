# -*- coding: utf-8 -*-


from __future__ import print_function
from __future__ import division

import numpy as np

from dolfin import *

# FIXME: Convert this to standalone test cases
# as the various functionality is implemented
def fenics_use_cases(backend, **kwargs):
    """Compact form of imagined fenics plotting use cases.

    To make this happen, we should create a self contained
    unit test for each plotting command here, and implement
    it for each wanted plotting backend.

    This is intended as a high level overview of relevant
    plotting types, and Details of **kwargs are deliberately
    skipped here. For each plotting function, a different set
    of parameters will be relevant for kwargs, although these
    parameters should be chosen consistently and preferrably
    consistently with existing plotting frameworks for user
    familiarity.
    """

    # --- Setup FEniCS object types that we want to plot

    # Tetrahedron mesh
    mesh = UnitCubeMesh(1, 1, 1)

    # Integer marker for each tetrahedron cell
    cell_markers = CellFunctionSizet(mesh)

    # Integer marker for each triangular facet
    face_markers = FacetFunctionSizet(mesh)

    # Integer marker for each edge
    edge_markers = EdgeFunctionSizet(mesh)

    # Integer marker for each vertex
    vertex_markers = VertexFunctionSizet(mesh)

    # Discrete function composed of discontinuous Lagrange polynomials of degree 0 over mesh
    V0 = FunctionSpace(mesh, "DP", 0)
    f0 = Function(V0)
    e0 = Expression("x[0]", degree=0)

    # Discrete function composed of continuous Lagrange polynomials of degree 1 over mesh
    V1 = FunctionSpace(mesh, "P", 1)
    f1 = Function(V1)
    e1 = Expression("x[0]", degree=1)

    # Discrete function composed of discontinuous Lagrange polynomials of degree 0 over mesh
    VD1 = FunctionSpace(mesh, "DP", 1)
    fd1 = Function(VD1)


    # --- Create a plotter object configured for a specific plotting backend

    p = DolfinPlotter(backend)


    # --- Examples of plotting commands

    # Note: A high level interface can automatically convert
    # e.g. a function to its mesh where that is expected,
    # but to keep things simple I suggest we leave that
    # for a later iteration.


    # --- Scatter plots (0D topology)

    # Scatter plot showing mesh vertices
    p.scatter(mesh, **kwargs)

    # Scatter plot placing values in cell midpoints
    p.scatter(f0, **kwargs)

    # Scatter plot placing values in vertices
    p.scatter(f1, **kwargs)


    # --- Wireframe plots (1D topology)

    # Wireframe showing mesh edges
    p.wireframe(mesh, **kwargs)

    # Wireframe showing mesh edges colored by marker values
    p.wireframe(edge_function, **kwargs)


    # TODO: Streamlines? Streaklines?


    # --- Surface plots (2D topology)

    # Surface showing mesh faces
    p.surface(mesh, **kwargs)

    # Surface showing face markers
    p.surface(face_markers, **kwargs)

    # Surface showing piecewise constant function on faces
    p.surface(f0, **kwargs)

    # Surface showing continuous function on faces
    p.surface(f1, **kwargs)

    # Surface showing discontinuous function on faces
    p.surface(fd1, **kwargs)


    # --- Isosurface plots (2D topology approximation of 3D data)

    # Isosurfaces of function in certain values
    p.isosurface(f1, **kwargs)

    # Note: consider projecting f0 and fd1 to V1,
    # at least as a first approximation


    # --- Glyph rendering of function

    # Instanced glyphs in set of points transformed via values of function
    # (which points must be configurable: vertices, cells, random)

    # Scalar fields can use e.g. spherical glyphs
    # Default to vertex coordinates for DG0 functions
    p.glyphs(e1, **kwargs)
    p.glyphs(f1, **kwargs)

    # Default to cell midpoints for DG0 functions
    p.glyphs(e0, **kwargs)
    p.glyphs(f0, **kwargs)

    # Vector fields can use e.g. line or array glyphs
    p.glyphs(v1, **kwargs)
    p.glyphs(v0, **kwargs)

    # Tensor fields can also be shown with some special glyphs
    p.glyphs(t1, **kwargs)
    p.glyphs(t0, **kwargs)


    # --- Displacement rendering (3D topology)

    # Displacement rendering of vector valued function
    p.displacement(v1, **kwargs)


    # --- Direct volume rendering (3D topology)

    # Direct volume rendering of mesh with fixed density
    #p.volume(mesh, **kwargs)

    # Direct volume rendering of function
    p.volume(f1, **kwargs)

    # Direct volume rendering of function
    #p.volume(f0, **kwargs)

    # Direct volume rendering of cell markers
    #p.volume(cell_markers, **kwargs)


    # --- Annotation plots (showing indices as text, useful for debugging)

    # Annotate cell indices of mesh
    p.annotate(mesh, dim=3, **kwargs)

    # Annotate face indices of mesh
    p.annotate(mesh, dim=2, **kwargs)

    # Annotate edge indices of mesh
    p.annotate(mesh, dim=1, **kwargs)

    # Annotate vertex indices of mesh
    p.annotate(mesh, dim=0, **kwargs)

    # Annotate marker values for each cell
    p.annotate(cell_markers, **kwargs)

    # Annotate marker values for each face
    p.annotate(face_markers, **kwargs)

    # Annotate marker values for each edge
    p.annotate(edge_markers, **kwargs)

    # Annotate marker values for each vertex
    p.annotate(vertex_markers, **kwargs)



def test():
    # TODO: Use generic mock framework to test the DolfinPlotter parts
    kwargs = {}
    backend = DummyPlotterBackend()
    fenics_use_cases(backend, **kwargs)
