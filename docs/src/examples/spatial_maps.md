# Spatial maps

When connecting multiple spatial models, we may be required to define a mapping function from the coordinate system of 1 component to another.
For example, a PDE component is solved through [MethodOfLines](@extref MethodOfLines index) which will have some specified discretization, but a connected Agent-based model may have continuous positions, or discrete positions on a different resolution.
We provide a handful of tools for mapping between two different spaces.

TODO
