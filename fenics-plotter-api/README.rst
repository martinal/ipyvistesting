Visualization types
===================
The following is a brainstorming of 3d visualization types we may
want to support


Organized by plot type
----------------------

- Scatter plot: points, scalar/vector/tensor values.
  Point positions are an important part of what's being communicated.

- Glyph plot: points, scalar/vector/tensor values.
  Point positions are only samples of a continuous domain.
  Vectors map to arrows.


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
