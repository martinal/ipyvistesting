# -*- coding: utf-8 -*-


from __future__ import print_function
from __future__ import division

import numpy as np


class MatplotlibBackend:
    def __init__(self):
        import matplotlib
        import matplotlib.pyplot
        self._matplotlib = matplotlib

    def plot(self, data, **kwargs):
        "Generic plotting pattern."
        plottype = data["plottype"]
        return getattr(self, plottype)(data, **kwargs)

    def scatter(self, data, **kwargs):
        pass

    def wireframe(self, data, **kwargs):
        pass

    def surface(self, data, **kwargs):
        plt = self._matplotlib.pyplot
        Triangulation = self._matplotlib.tri.Triangulation
        ax = plt.gca()

        gdim = data["gdim"]
        tdim = data.get("tdim", gdim)
        coordinates = data["coordinates"]
        cells = data["cells"]

        if tdim != 2:
            raise ValueError("Expecting 2D data in surface plot.")

        point_values = data["point_values"]
        face_values = data["face_values"]

        if point_values:
            fixme
        elif face_values:
            ax.tripcolor(triang, facecolors=face_values, **kwargs)
        else:
            if gdim == 2:
                color = kwargs.pop("color", '#808080')
                triangulation = Triangulation(coordinates[:, 0], coordinates[:, 1], cells)
                return ax.triplot(triangulation, color=color, **kwargs)
            elif gdim == 3:
                return ax.plot_trisurf(*[coordinates[:,i] for i in range(gdim)],
                                       triangles=cells, **kwargs)
            else:
                raise ValueError("%dD not supported." % (gdim,))

    def interactive(self):
        plt = self._matplotlib.pyplot
        plt.show()


from pprint import pprint

class PrintBackend:

    def plot(self, data, **kwargs):
        "Generic plotting pattern."
        plottype = data["plottype"]
        print("****** Plotting data with plottype '%s':" % (plottype,))
        pprint(data)
        if kwargs:
            print("using kwargs:")
            pprint(kwargs)
        if hasattr(self, plottype):
            return getattr(self, plottype)(data, **kwargs)

    def scatter(self, data, **kwargs):
        print("Calling scatter.")

    def interactive(self):
        print("No interaction in print backend.")


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
    d = elm.degree()
    f = elm.family()
    print(d, f)
    return elm.degree() > 0 and elm.family() == "Lagrange"


class Extractor:

    def _extract_mesh(self, obj, mesh=None):
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

    def _try_project(self, obj, mesh):
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

    def determine_plottype(self, obj, plottype, **kwargs):
        """Automatically detect or override plottype to be compatible with this data.
        
        If no plottype can be determined, raise ValueError.
        """
        if plottype is None:
            plottype = "scatter"
        return plottype

    def extract(self, obj, plottype, **kwargs):
        "Extract data for object."
        plottype = self.determine_plottype(obj, plottype)
        if not hasattr(self, plottype):
            raise ValueError("Missing implementation of plottype '%s'." % (plottype,))
        data = getattr(self, plottype)(obj, **kwargs)
        data["plottype"] = plottype
        return data

    def scatter(self, obj, **kwargs):
        """Extract data for a scatter plot.

        @param obj
            If obj is a Mesh, shows mesh vertices with constant values.
            If obj is a Function, shows function values in mesh vertices or cell midpoints.
            If obj is a MeshFunction, shows marker values on mesh entity midpoints.

        @return data
            data["points"]
            data["point_values"]
        """
        mesh = self._extract_mesh(obj, mesh=kwargs.get("mesh"))
        gdim = mesh.geometry().dim()
        tdim = mesh.topology().dim()
        coordinates = mesh.coordinates()
        cells = mesh.cells()

        if isinstance(obj, cpp.Mesh):
            points = coordinates
            point_values = None
            name = "Mesh with id=%d" % (obj.id(),)

        elif isinstance(obj, (cpp.Expression, cpp.Function)):
            name = "Function with id=%d" % (obj.id(),)
            if has_continuous_element(obj):
                points = coordinates
                point_values = obj.compute_vertex_values(mesh)
            else:
                raise NotImplementedError("Need cell midpoint computation in dolfin to cover this case.")
                points = mesh.cell_midpoints()
                point_values = obj.compute_cell_midpoint_values(mesh)

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
                raise NotImplementedError("Scatter plot data only implemented for vertex based meshfunctions.")
        else:
            raise NotImplementedError("Scatter plot data not implemented for %s." % (obj.__class__,))

        # This data is what the plotting backend gets
        data = dict(
            name=name,
            tdim=tdim,
            gdim=gdim,
            points=points,
            point_values=point_values
        )
        return data

    def wireframe(self, obj, **kwargs):
        """Extract data for a wireframe plot.
        """
        # FIXME
        data = dict(
        )
        return data

    def surface(self, obj, **kwargs):
        """Extract data for a surface plot.
        """
        mesh = self._extract_mesh(obj, mesh=kwargs.get("mesh"))
        gdim = mesh.geometry().dim()
        tdim = mesh.topology().dim()

        if isinstance(obj, cpp.Mesh):
            name = "Mesh with id %d" % (mesh.id(),)
            if tdim == 3:
                assert gdim == 3
                name = "Boundary surface of " + name
                mesh = BoundaryMesh(mesh, "exterior", order=False)
                tdim = 2
            coordinates = mesh.coordinates()
            cells = mesh.cells()
            point_values = None
            face_values = None

        elif isinstance(obj, (cpp.Expression, cpp.Function)):
            name = "Function with id=%d" % (obj.id(),)
            if tdim > 2:
                raise NotImplementedError("Missing implementation of surface of tdim == 3.")

            if has_continuous_element(obj):
                point_values = obj.compute_vertex_values(mesh)
                face_values = None
            else:
                point_values = None
                face_values = obj.compute_facet_midpoint_values(mesh)  # Missing function
                raise NotImplementedError("Need facet midpoint computation in dolfin to cover this case.")
        else:
            raise NotImplementedError("Surface plot data not implemented for %s." % (obj.__class__,))

        # This data is what the plotting backend gets
        data = dict(
            name=name,
            tdim=tdim,
            gdim=gdim,
            coordinates=coordinates,
            cells=cells,
            point_values=point_values,
            face_values=face_values,
        )
        return data

    # Etc. for each supported plottype


class Plotter:
    """Initial draft of a plotter object for dolfin object.

    To be used as design input for wanted features for a Jupyter plotting backend.
    """
    def __init__(self, backend, extractor, default_kwargs=None):
        # Injected dependency on plotting backend
        self.backend = backend
        self.extractor = extractor
        self.default_kwargs = {} if default_kwargs is None else dict(default_kwargs)

    def interactive(self):
        return self.backend.interactive()

    def plot(self, obj, plottype=None, **kwargs):
        "Generic plotting pattern."
        completed_kwargs = dict(self.default_kwargs)
        completed_kwargs.update(kwargs)
        data = self.extractor.extract(obj, plottype, **completed_kwargs)
        assert "plottype" in data
        return self.backend.plot(data, **completed_kwargs)

    def __call__(self, obj, plottype=None, **kwargs):
        return self.plot(obj, plottype, **kwargs)

    def scatter(self, obj, **kwargs):
        return self.plot(obj, "scatter", **kwargs)

    def wireframe(self, obj, **kwargs):
        return self.plot(obj, "wireframe", **kwargs)

    def surface(self, obj, **kwargs):
        return self.plot(obj, "surface", **kwargs)


def main():
    backend = PrintBackend()
    backend = MatplotlibBackend()
    extractor = Extractor()
    plotter = Plotter(backend, extractor, {})

    #plot = plotter.plot
    #interactive = plotter.interactive

    mesh = UnitSquareMesh(5, 3)
    #plotter.plot(mesh)
    plotter.surface(mesh)
    #plotter.scatter(mesh)
    #plotter.wireframe(mesh)

    if 0:
        V0 = FunctionSpace(mesh, "DG", 0)
        f0 = Function(V0)
        plotter.plot(f0)

    if 0:
        V1 = FunctionSpace(mesh, "CG", 1)
        f1 = Function(V1)
        plotter.plot(f1)

    if 0:
        V1d = FunctionSpace(mesh, "DG", 1)
        f1d = Function(V1d)
        plotter.plot(f1d)

    if 0:
        V2 = FunctionSpace(mesh, "CG", 2)
        f2 = Function(V2)
        plotter.plot(f2)

    plotter.interactive()


if __name__ == "__main__":
    main()
