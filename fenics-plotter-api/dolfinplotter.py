# -*- coding: utf-8 -*-


from __future__ import print_function
from __future__ import division

import numpy as np

from dolfin import *
from dolfin import cpp

_meshfunction_types = (
    cpp.MeshFunction,
    cpp.MeshFunctionBool,
    cpp.MeshFunctionInt,
    cpp.MeshFunctionDouble,
    cpp.MeshFunctionSizet
    )

_plottable_types = (
    cpp.Function,
    cpp.Expression,
    cpp.Mesh,
    cpp.MultiMesh,
    cpp.DirichletBC
    ) + _meshfunction_types


def has_continuous_element(obj):
    elm = obj.ufl_element()
    return elm.degree() > 0 and elm.family() == "Lagrange"


class DolfinPlotDataExtractor:

    @classmethod
    def _extract_mesh(cls, obj, mesh=None):
        mesharg = mesh
        if isinstance(obj, cpp.Mesh):
            mesh = obj
        elif isinstance(obj, cpp.Function):
            mesh = obj.function_space().mesh()
        elif hasattr(obj, "mesh"):
            mesh = obj.mesh()
        else:
            mesh = mesharg

        if mesharg is not None and mesharg.id() != mesh.id():
            raise ValueError("Passed Mesh is not the same as object mesh.")
        if not isinstance(mesh, cpp.Mesh):
            raise ValueError("Expecting a Mesh argument or object having a Mesh.")

        return mesh


    @classmethod
    def _try_project(cls, obj, mesh):
        from dolfin.fem.projection import project
        try:
            cpp.info("Object cannot be plotted directly, projecting to"
                     " piecewise linear function space.")
            obj = project(obj, mesh=mesh)
            mesh = obj.function_space().mesh()
        except Exception as e:
            msg = ("Don't know how to plot given object:\n  %s\n"
                   "and projection failed:\n  %s")
            raise ValueError(msg % (str(object), str(e)))
        return obj, mesh


    # FIXME: Figure out how much the decision of which data to extract is plot
    # mode specific or if it can be covered by a separate plot mode agnostic API,
    # something like:
    def extract_entity_coords_and_values(self, obj, entity_dim):
        """
        
        FIXME

        What would make this easier are some functions:

            coordinates = mesh.compute_entity_midpoint_coordinates(entity_dim)
            values = function.compute_entity_midpoint_values(entity_dim[, mesh])

        which for entity_dim == 0 falls back to vertices,
        and entity_dim == tdim falls back to cell midpoints.
        I think vertices, cells, and facets will all have
        fairly common use cases in visualization.

        Maybe back it up with something like this:

            ufc::coordinate_mapping::tabulate_reference_entity_midpoint_coordinates(X, entity_dim)
            ufc::coordinate_mapping::compute_entity_midpoint_geometry(X, J, entity_dim)
            ufc::finite_element::evaluate_reference_entity_midpoint_basis(X, entity_dim)
            ufc::finite_element::evaluate_entity_midpoint_basis(X, entity_dim)

        """

    def extract_scatter_data(self, obj, **kwargs):
        """Extract arrays with data for a scatter plot.

        @param obj
            If obj is a Mesh, shows mesh vertices with constant values.
            If obj is a Function, shows function values in mesh vertices or cell midpoints.
            If obj is a MeshFunction, shows marker values on mesh entity midpoints.

        @return points
        @return point_values
        """
        mesh = self._extract_mesh(obj, mesh=kwargs.get("mesh"))
        coordinates = mesh.coordinates()
        cells = mesh.cells()

        if isinstance(obj, cpp.Mesh):
            points = coordinates
            point_values = None
        elif isinstance(obj, (cpp.Expression, cpp.Function)):
            if has_continuous_element(obj):
                points = coordinates
                point_values = obj.compute_vertex_values(mesh)
            else:
                #points = mesh.cell_midpoints()
                #point_values = obj.compute_cell_midpoint_values(mesh)
                raise NotImplementedError("Need cell midpoint computation in dolfin to cover this case.")
        elif isinstance(obj, _meshfunction_types):
            dim = obj.dim()
            point_values = obj.array()
            if dim == 0:
                points = coordinates
            #elif obj.dim() == 1:
            #    points = FIXME  # edge midpoints
            #elif obj.dim() == 2:
            #    points = FIXME  # face midpoints
            #elif obj.dim() == 3:
            #    points = FIXME  # volume cell midpoints
            else:
                raise NotImplementedError("Scatter only implemented for vertex based meshfunctions.")
        else:
            raise NotImplementedError("Scatter not implemented for %s." % (obj.__class__,))


class DolfinPlotter:
    """Initial draft of a plotter object for dolfin object.

    To be used as design input for wanted features for a Jupyter plotting backend.
    """
    def __init__(self, backend, extractor):
        # Injected dependency on plotting backend
        self.backend = backend
        self.extractor = extractor


    def plot(self, obj, plottype=None, **kwargs):
        "Generic plotting pattern."
        data = self.extractor.extract(obj, plottype, **kwargs)
        return self.backend.plot(data, plottype, **kwargs)


    def scatter(self, obj, **kwargs):
        """Scatter plot.

        @param obj
            If obj is a Mesh, shows mesh vertices with constant values.
            If obj is a Function, shows function values in mesh vertices or cell midpoints.
            If obj is a MeshFunction, shows marker values on mesh entity midpoints.
        """
        points, point_values = self.extract_scatter_data(obj, **kwargs)
        return self.backend.scatter(points, point_values, **kwargs)


    def wireframe(self, obj, **kwargs):
        """Wireframe plot.

        @param obj
            If obj is a Mesh, shows wireframe of mesh edges.
            If obj is a Function, colors mesh edges by function values in mesh vertices.
            If obj is a MeshFunction, colors mesh edges by marker values on mesh entities.
        """
        return self.plot(obj, "wireframe", **kwargs)


    def surface(self, obj, **kwargs):
        """Surface plot.

        @param obj
            If obj is a Mesh, shows mesh faces.
            If obj is a Function, colors mesh faces by function values in mesh vertices or cells.
            If obj is a MeshFunction, colors mesh faces by marker values on mesh faces or cells.
        """
        data = self.extract("surface", obj, **kwargs)
        return self.backend.surface(data.coordinates, data.triangles, **kwargs)


    def isosurface(self, obj, **kwargs):
        """Isosurface plot.

        @param obj
            If obj is a Function, shows isosurfaces of that function.
        @param values
            List of isovalues.
        """
        raise NotImplementedError("FIXME")


    def quiver(self, obj, **kwargs):
        """Quiver plot, showing instanced glyphs in set of points transformed via values of function.

        TODO: Which points to select should be configurable: vertices, cells, random.

        Glyph representation is left up to backend.

        @param obj
            A Function, scalar, vector or tensor valued:
            - scalar fields can use e.g. spherical glyphs
            - vector fields can use e.g. line or array glyphs
            - tensor fields can also be shown with some special glyphs
        """
        raise NotImplementedError("FIXME")


    def extract_displacement_data(self, obj, **kwargs):
        """Displacement rendering of vector field.

        @param obj
            A vector valued Function of same value size
            as the geometric dimension of the mesh.
        """
        if not isinstance(obj, (cpp.Expression, cpp.Function)):
            raise NotImplementedError("Displacement not implemented for %s." % (obj.__class__,))
        elif not has_continuous_element(obj):
            raise NotImplementedError("Invalid element for displacement function.")

        mesh = self.extract_mesh(obj, mesh=kwargs.get("mesh"))
        vertices = mesh.coordinates()
        cells = mesh.cells()
        vertex_values = obj.compute_vertex_values(mesh)

        # TODO: Extract boundary mesh instead
        return vertices, cells, vertex_values


    def displacement(self, obj, **kwargs):
        """Displacement rendering of vector field.

        @param obj
            A vector valued Function of same value size
            as the geometric dimension of the mesh.
        """
        vertices, cells, vertex_values = self.extract_displacement_data(obj, **kwargs)
        return self.backend.displacement(vertices, cells, vertex_values, **kwargs)


    def warp(self, obj, **kwargs):
        """Displacement rendering of scalar height field in 2D.

        @param obj
            A scalar valued Function which will be used as the height.
        """
        raise NotImplementedError("FIXME")


    def volume(self, obj, **kwargs):
        """Direct volume rendering.

        @param obj
            If obj is a Mesh, render using a fixed density.
            If obj is a Function, render using density from function values.
        """
        raise NotImplementedError("FIXME")


    def annotate(self, markers, **kwargs):
        """Annotate mesh entity indices or marker values with numbers.

        @param markers
            A dolfin.MeshFunction or dolfin.Mesh

        @param dim
            If markers is a dolfin.Mesh, dim is the entity dimension to annotate
        """
        raise NotImplementedError("FIXME")
