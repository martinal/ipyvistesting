
# Sketch of possible interface for PT code
# with FEniCS on one side and a Jupyter widget
# using THREE.js/WebGL on the other side.
# Most names here are quickly invented without
# too much thought, so feel free to use other
# names or modify the interface


# Some initial imports
from IPython.display import display
import ipywidgets as widgets
import numpy as np


"""
# Setup a fixed dolfin function
from dolfin import *
mesh = UnitCubeMesh(1, 1, 1)
V = FunctionSpace(mesh, "Lagrange", 1)
e = Expression("(x[0]-x0) *(x[1]-x1) * (x[2]-x2)",
               x0=0.5, x1=0.5, x2=0.5, degree=1)
f = Function(V)
f.interpolate(e)
"""


"""
# Compute sorted Projected Tetrahedra triangles
from dolfinpt import *
pt = ProjectedTetrahedraBuilder()
pt.set_mesh(mesh)
pt.set_function(f)
pt.set_mvp(MVP)


# Get shaders from pt object
vertex_shader = pt.get_vertex_shader();
fragment_shader = pt.get_fragment_shader();

# Get coordinates, triangle indices, attributes from pt object
coordinates = pt.get_coordinates()
triangles = pt.get_triangles()
names = pt.get_attribute_names()
attributes = { pt.get_attribute(name) for name in pt.get_attribute_names() }
uniforms = { pt.get_uniform(name) for name in pt.get_uniform_names() }
"""


# Define shaders directly
vertex_shader = """
#version 130

in vec3 vertexPosition;
in vec4 vertexColor;

out vec4 color;

void main(){
    color = vertexColor;
    gl_Position = vec4(vertexPosition, 1.0);
}
"""

fragment_shader = """
#version 130

in vec4 color;

out vec4 fragmentColor;

void main(){
    fragmentColor = color;
}
"""


# Define mesh and attributes directly
coordinates = np.asarray([
    (0.0, 0.0),
    (1.0, 0.0),
    (1.0, 1.1),
    ], dtype="float32")
num_vertices = coordinates.shape[0]

triangles = np.asarray([
    (0, 1, 2),
    ], dtype="int32")
num_triangles = triangles.shape[0]

attributes = {}
attributes["f_front"] = np.asarray([0.0, 0.1, 0.2], dtype="float32")
attributes["f_back"]  = np.asarray([0.9, 0.8, 0.7], dtype="float32")
attributes["s_front"] = np.asarray([1.0, 0.5, 0.5], dtype="float32")
attributes["s_back"]  = np.asarray([0.0, 0.5, 0.5], dtype="float32")

uniforms = {
    "f_min": 0.0,
    "f_max": 1.0,
    "color": [1.0, 1.0, 1.0],
    }

# TODO: Preintegrated textures


# Setup widget using THREE.js and WebGL backend for rendering
from volrenwidget import *
vrw = VolRenWidget()
vrw.set_shaders(vertex_shader, fragment_shader)
vrw.set_triangles(coordinates, triangles)
vrw.set_uniforms(uniforms)
vrw.set_attributes(attributes)


# Finally display widget in notebook
display(vrw)


# Then loop over timesteps and update
# display with new function values!
for t in np.arange(0.0, 1.0, 10):
    e.x0 = 0.5 * t
    f.interpolate(e)
    pt.update_function(f)

    attributes["f_front"][:] = fixme
    attributes["f_back"][:] = fixme

    vrw.resend(["f_front", "f_back"])
