# NetworkDynamics.jl

## Overview

This package implements functions for defining and studying dynamics on networks.
The key construction is a callable struct compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = network_dynamics(nodes!, lines!, s_e, t_e, dim_n, dim_l; symbols_n=nothing, symbols_l=nothing)
nd(dx, x, p, t)
```

The first two parameters are the functions, or function arrays from which a network dynamics is
constructed:

``$nodes![n](dx, x, [l]_s, [l]_t, p, t)$``
``$lines![e](dl, l, x_s, x_t, p, t)$``

The arrays dim_n and dim_l encode the dimensionality of $x$ and $l$ variables.
The arrays s_e and t_e encode the network structure by giving the source and
target of each edge.

Optionally we can also specify an array of symbols per edge and node that allow
convenience access to the different nodal and line dimensions.

Given edges $e$, and nodes $n$, as well as an orientation encoded by
the source function $s(e)$ and the target function $t(e)$
this implements the system of ODEs:

``$\frac{dx_n}{dt} = dx_n$``

``$\frac{dl_e}{dt} = dl_e$``

with $dx$ and $dl$ calculated by

``$[l]_s = [l_e \text{ if } s(e) = n]$``

``$[l]_t = [l_e \text{ if } t(e) = n]$``

``$nodes![n](dx_n, x_n, [l]_s, [l]_t, p_n, t)$``

``$lines![e](dl_e, l_e, x_{s(e)}, x_{t(e)}, p_e, t)$``

Something that relaxes to a diffusive network would for example be
implemented by


```julia
lines = (dl, l, x_1, x_2) -> dl .= 1000. * ((x_1 - x_2) - l)
nodes = (dx_n, x_n, l_s, l_t, p_n, t) -> dx_n .= f(x_n) - (sum(l_s) - sum(l_t))
```

The package also supplies a node and a line type that combine the the node
function, the dimensionality and the symbols, as well as a constructor for the
network dynamics from these types. Further Constructors are provided for
LightGraphs:

```julia
nd=network_dynamics(nodes::Array{ODE_Node}, lines::Array{ODE_Line}, G::AbstractGraph)
```

## Static lines

For static line relations we similarly have:

```julia
sl_nd = static_lines_network_dynamics(nodes!, lines!, s_e, t_e, ...)
sl_nd(dx, x, p, t)
```

With the convention for lines given by:

``$lines![e](l, x_s, x_t, p, t)$``

and otherwise as above. This implements the system of ODEs:

``$\frac{dx_n}{dt} = dx_n$``

with $dx$ calculated by

``$lines![e](l_e, x_{s(e)}, x_{t(e)}, p_e, t)$``

``$[l]_s = [l_e \text{ if } s(e) = n]$``

``$[l]_t = [l_e \text{ if } t(e) = n]$``

``$nodes![n](dx_n, x_n, [l]_s, [l]_t, p_n, t)$``

A diffusive network would be implemented by

```julia
lines = (l, x_1, x_2) -> l .= x_1 - x_2
nodes = (dx_n, x_n, l_s, l_t, p_n, t) -> dx_n .= f(x_n) - (sum(l_s) - sum(l_t))
```

The alternative constructor is given by:

```julia
sl_nd=static_lines_network_dynamics(nodes::Array{ODE_Node}, lines::Array{Static_Line}, G::AbstractGraph)
```

## Network DAEs

Design question: Don't do implicit DAEs, support mass matrices everywhere by
default. This adds two more (optional) arrays of vectors to the signature of
the constructor and the ODE_Node and the Static_Line types.

The advantage is that we can then deal with only three types of dynamics, ODE,
SDE and DDE. Leaning towards yes on this. The important special case of all ODE
can then be done via a performance optimization.

This would suggest splitting this into three sub packages, ODE, SDE and DDE.
Each sub package should define its own XDE_Node and XDE_line as well as
promotion rules from other XDE_Node and XDE_line types.

## Network SDEs
## Network DDEs

## Convenience functions for symbolic access to node variables

## API

```@autodocs
Modules = [NetworkDynamics]
```
