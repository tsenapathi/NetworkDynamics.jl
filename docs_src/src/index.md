# NetworkDynamics

# Overview

This package implements functions for defining and studying dynamics on networks.
The key construction is a callable struct compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = network_dynamics(vertices!, edges!, g)
nd(dx, x, p, t)
```

The first two parameters are the functions, or function arrays from which a network dynamics is
built. The last parameter g is a graph encoding the network constructed with
the LightGraphs.jl package.
