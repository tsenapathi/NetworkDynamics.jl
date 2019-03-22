# Functions

The Dynamics for the whole Network is constructed from functions for the single vertices and edges. There are several types:

```@docs
ODEVertex(vertexfunction!, dimension, massmatrix, sym)
StaticEdge(edgefunction!, dimension)
ODEEdge(edgefunction!, dimension, massmatrix, sym)
```


### ODEVertex

The arguments mean the following: vertexfunction! is catching the dynamics of a single vertex depending on the vertex value itself as well as in- and outgoing currents (or edges). An example for such a function would be:

```julia
function vertexfunction!(dv, v, e_s, e_d, p, t)
  dv .= 0
  for e in e_s
    dv .-= e
  end
  for e in e_d
    dv .+= e
  end
end
```

The e_s and e_d are arrays containing the edges that have the decribed vertex as source and destination. Other arguments coincide with the usual ODE function arguments. The vertexfunction given to ODEVertex always needs to have the shown argument structure. Note the importance of the broadcast structure of the equations (the dot before every operator), this is necessary due to the use of views in the internal functions.

dimension is the number of Variables on the Vertex.

massmatrix is the mass matrix M, i.e.

```@docs
M*dv = vertexfunction!
```

sym are the symbols of the Vertex. If one had for example a vertex with a frequency and some angle, one would construct sym via:

```@docs
sym = [:omega, :phi]
```

This makes it easier to later fish out the interesting variables one wants to look at.

One may also call ODEVertex via:

```@docs
ODEVertex(vertexfunction!, dimension)
```

The function then defaults to using the identity as mass matrix and [:v for i in 1:dimension] as symbols.


### StaticEdge

Static here means, that the edge value described by edgefunction! solely depends on the vertex values the edge connects. One very simple and natural example is a diffusive system:

```@julia
edgefunction! = (e, v_s, v_d, p, t) -> e .= v_s .- v_d
```

v_s and v_d are the vertex values of the edges source and destination. There is no derivative of the edge value involved, hence we call these problems static.

dimension: see ODEVertex

### ODEEdge

For Problems where the derivative of the edge value indeed depends on the edge value itself, we use the ODEEdge function. Another simple and natural example is a system that quickly diffuses to the static case:

```@julia
edgefunction! = (de, e, v_s, v_d, p, t) -> de .= 1000 * (v_s .- v_d .- e)
```

dimension: see ODEVertex

massmatrix: see ODEVertex

sym: see ODEVertex

Also, one can construct an ODEEdge by only giving the first two arguments:

```@docs
ODEEdge(edgefunction!, dimension)
```

Then the function defaults to using the identity as mass matrix as well as using [:e for in 1:dimension] as sym.




## Constructor

The central constructor of this package is network_dynamics(), this function demands an array of VertexFunction and EdgeFunction as well as a graph (see LightGraphs), and returns an ODEFunction which one can easily solve via the tools given in DifferentialEquations.jl. One calls it via:

```@docs
network_dynamics(Array{VertexFunction}, Array{EdgeFunction}, graph)
```
