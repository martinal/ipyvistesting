Visualization types
===================

The following is a brainstorming of 3d visualization types we may want to support.


Data types
----------

- Mesh (points and cells)
  - gdim = geometric dimension
  - tdim = topological dimension
  - tdim <= gdim
  - Point cloud (tdim = 0)
  - Wireframe or other set of line segments (tdim = 1)
  - Surface or other set of planar polygons (tdim = 2)
  - Volume (tdim = 3)

- Value:
  - Single data point, scalar, vector or tensor valued
    - dtype (int, float, etc)
    - shape
    - rank = len(shape)

- Mesh values:
  - Value assiciated with each entity in mesh of specific tdim
  - Most common cases:
    - Point values:
      - data assiciated with each point in mesh
    - Cell values:
      - data assiciated with each cell in mesh


Organized by plot type
----------------------

- Single glyph: value
  Scalars map to points with directionless properties such as color and/or size.
  Vectors map to arrows or lines.

- Scatter plot: points, optionally point values.
  Point positions are an important part of what's being communicated.
  Each value is separately translated to a glyph.

- Glyph plot: points, point values.
  Point positions are only samples of a continuous domain.

- Line plot (wireframe etc): points, edges, [point values | edge values]

- Height plot: mesh (points, cells), point scalar values
  - Regular line plot
  - Height surface

- Displacement plot: mesh (points, cells), point vector values
  - Surface displacement
  - 3D displacement plot


Organized by field type
-----------------------

- Point clouds:
  - Basic:
    - Glyphs

- Mesh:
  - Basic:
    - Scatter plot of vertex positions
    - Wireframe

- Scalar fields:
  - Basic:
    - Glyphs
    - Surface color
  - Intermediate:
    - Surface color plot after some postprocessing:
      - Boundary
      - Clip plane
      - Isosurface
      - Any kind of surface that postprocessing can produce
  - Advanced:
    - Volume rendering variants:
      - Full emission-absorption model
      - Splatting
      - Maximum intensity projection
      - Arbitrary number of isosurfaces
      - Stochastic scatter plots

- Vector fields:
  - Basic:
    - Glyphs
    - Displacement
    - Scalar field visualizations of derived quantities (magnitude, directional components)
  - Intermediate:
    - Streamlines
    - Streaklines
  - Advanced:
    - LIC
    - FTLE

- Tensor fields:
  - Basic:
    - Glyphs
    - Scalar field visualizations of derived quantities
