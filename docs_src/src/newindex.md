# NetworkDynamics.jl

## Overview

This package implements functions for defining and studying dynamics on networks.
The key construction is a callable struct compatible with the
DifferentialEquations.jl calling syntax.

```julia
nd = network_dynamics(nodes!, lines!, g)
nd(dx, x, p, t)
```

The first two parameters are the functions, or function arrays from which a network dynamics is
built. The last parameter is a graph encoding the network, constructed with
the LightGraphs.jl package.
Note that the type network_dynamics is only a placeholder for the dynamics types we will now specify.

##Static lines

A dynamical network with static lines (static meaning that the current on an edge depends solely on the
values on the nodes it connects) is created via the struct static_line_network_dynamics:

```julia
slnd = static_line_network_dynamics(nodes!, lines!, g)
```
The functions nodes! and lines! are of the form:

```julia
nodes![n](dx[n],x[n],l_s[n],l_t[n],p,t)
lines![e](l[e],x_s,x_t,p,t)  
```

Specifically, the given variables are:

```julia
l_s[n] = [l[e] if s[e] == n for e in 1:length(lines!)]
l_t[n] = [l[e] if t[e] == n for e in 1:length(lines!)]
x_s=x[s[e]]
x_t=x[t[e]]
```
The vectors s and t contain the information about the source and target of each
edge, i.e. s[1] == 2 -> The source of edge 1 is node 2. The vectors l_s[n] and l_t[n] are
basically containing the in- and outgoing currents of node n in the form of an array. Thus,
when describing the dynamics of the node variables, we will need to sum over these.

For example, a system of equations describing a diffusive network would be:

```julia
using LightGraphs
g= barabasi_albert(10,5)
nodes! = [(dx,x,l_s,l_t,p,t) -> dx .= sum(l_s) - sum(l_t) for n in nodes(g)]
lines! = [(l,x_s,x_t,p,t) -> l .= x_s - x_t for e in edges(g)]
```

Here, the diffusiveness lies within the lines! function. It states that there is only
a current between two nodes if these nodes have a non-equal value. This current then ultimatively
leads to an equilibrium in which the value on any connected node is equal.

Note that one needs to put a dot before the equal signs. This is due to the use of views
in the internals of the static_line_network_dynamics function.

We finally want to solve the defined diffusive system. This we do by using the magnificient
package DifferentialEquations.jl. We also need to specify a set of initial values x0 as well as a time
interval for which we are solving the problem:

```julia
using DifferentialEquations
x0 = rand(10)
slnd = static_line_network_dynamics(nodes!,lines!,g)
slnd_prob = ODEProblem(slnd,x0,(0.,2.))
sol = solve(slnd_prob)
```
