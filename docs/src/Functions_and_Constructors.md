# Functions

The key constructor [`network_dynamics`](@ref) assembles the dynamics of the whole network from functions for the single vertices and edges of the graph `g`.
Since the equations describing the local dynamics may differ strongly from each
other, the types `VertexFunction` and `EdgeFunction` are introduced. They
provide a unifying interface between different classes of nodes and edges. Both
have several subtypes that account for the different types of equations that may
represent the local dynamics. At the moment algebraic (static) equations and ordinary differential equations (ODEs) are supported:


```julia
# VertexFunctions
StaticVertex(f!, dimension, symbol)
ODEVertex(f!, dimension, mass_matrix, symbol)

# EdgeFunctions
StaticEdge(f!, dimension, symbol)
ODEEdge(f!, dimension, mass_matrix, symbol)
```


# VertexFunctions

Given a set of (algebraic or differential) equations describing a node or an edge
the first step is to turn them into a **mutating** function `f!`. Depending on the class of the function `f!`, the constructors `StaticVertex` or `ODEVertex` are called in order to turn `f!` into a `VertexFunction` object compatible with [`network_dynamics`](@ref).


Since in general the state of a vertex depends on the vertex value itself as well as on the in- and outgoing edges, the function `f!`
has to respect on of the following calling syntaxes.

```julia
# For static nodes
function f!(v, e_s, e_d, p, t) end
# For dynamics nodes
function f!(dv, v, e_s, e_d, p, t) end
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while `e_s` and `e_d` are arrays containing the edges for which the described vertex is the source or the destination respectively. The typical case of diffusive coupling on a directed graph could be described as

```julia
function vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0.
    for e in e_s
        dv .-= e
    end
    for e in e_d
        dv .+= e
    end
    nothing
end
```
!!! warning
    The arguments `e_s` and `e_d` are **obligatory** even if the graph is undirected and no distinction between source and destination can be made. This is necessary since [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) implements an undirected graph in the same way as a directed graphs, but ignores the directionality information. Therefore some care has
    to be taken when dealing with assymetric coupling terms. A detailed example
    can be found [here](missing).

### [StaticVertex](@ref)

If a vertex is described by an algebraic equation  `f!(v, e_s, e_d, p, t)`, i.e. `dv = 0` the `VertexFunction` is constructed as

```julia
StaticVertex(f!, dim, sym)
```

Here, **dim** is the number of independent variables in the vertex equations and **sym** is an array of symbols of these variables. For example, if a node
models a constant input ``I = p``, then `dim = 1` and `sym = [:I]`. For more details on the use of symbols, check out the [example](missing) section or the Julia [documentation](https://docs.julialang.org/en/v1/manual/metaprogramming/). The use of symbols makes ODEVit easier to later fish out the interesting variables one wants to look at.


### [ODEVertex](@ref)

If a vertex has local dynamics `f!(dv, v, e_s, e_d, p, t)` described by an ODE
the `VertexFunction` is contructed as

```julia
ODEVertex(f!, dim, mass_matrix, sym)
```

As above, **dim** is the number of independent variables in the vertex equations and **sym** corresponds to the symbols of these variables.

**mass_matrix** is an optional argument that defaults to the identity matrix `I`. If a mass matrix M is given the system ``M * dv = f!`` will be solved.


One may also call ODEVertex with keyword arguments:

```julia
ODEVertex(f! = f!, dim = dim)
```

The function then defaults to using the identity as mass matrix and `[:v for i in 1:dimension]` as symbols.


### StaticEdge

Static here means, that the edge value described by **edgefunction!** solely depends on the vertex values the edge connects. One very simple and natural example is a diffusive system:

```@julia
edgefunction! = (e, v_s, v_d, p, t) -> e .= v_s .- v_d
```

v_s and v_d are the vertex values of the edges source and destination. There is no derivative of the edge value involved, hence we call these problems static.

**dimension**: see ODEVertex

### ODEEdge

For Problems where **edgefunction** also contains the differential of an edge value , we use the ODEEdge function. Another simple and natural example for such a system is one that quickly diffuses to the static case:

```@julia
edgefunction! = (de, e, v_s, v_d, p, t) -> de .= 1000 * (v_s .- v_d .- e)
```

**dimension**: see ODEVertex

**mass_matrix**: see ODEVertex

**sym**: see ODEVertex

Also, one can construct an ODEEdge by only giving the first two arguments:

```julia
ODEEdge(edgefunction!, dimension)
```

Then the function defaults to using the identity as mass matrix as well as using [:e for in 1:dimension] as sym.




## Constructor

The central constructor of this package is network_dynamics(), this function demands an array of VertexFunction and EdgeFunction as well as a graph (see LightGraphs), and returns an ODEFunction which one can easily solve via the tools given in DifferentialEquations.jl. One calls it via:

```julia
network_dynamics(Array{VertexFunction}, Array{EdgeFunction}, graph)
```

VertexFunction and EdgeFunction are the Unions of all the Vertex and Edge Functions we specified in the previous section. Let's look at an example. First we define our graph as well as the differential systems connected to its vertices and edges:

```julia

using LightGraphs

g = barabasi_albert(10,5) # The graph is a random graph with 10 vertices and 25 Edges.

function vertexfunction!(dv, v, e_s, e_d, p, t)
  dv .= 0
  for e in e_s
    dv .-= e
  end
  for e in e_d
    dv .+= e
  end
end

function edgefunction! = (de, e, v_s, v_d, p, t) -> de .= 1000 .*(v_s .- v_d .- e)

vertex = ODEVertex(vertexfunction!, 1)
vertexarr = [vertex for v in vertices(g)]

edge = ODEEdge(edgefunction!, 1)
edgearr = [edge for e in edges(g)]

nd = network_dynamics(vertexarr, edgearr, g)
```

Now we have an ODEFunction nd that we can solve with well-known tools from DifferentialEquations. To solve the defined system,
we further need an array with initial values x0 as well as a time span tspan in which we solve the problem:

```julia

using DifferentialEquations

x0 = rand(10 + 25) #10 for the vertices and 25 for the edges
tspan = (0.,2.)

prob = ODEProblem(nd,x0,tspan)
sol = solve(prob)

using Plots
plot(sol, legend = false, vars = 1:10) # vars gives us x[1:10] in the plot
```

The Plot shows the classic diffusive behaviour.

### Mass Matrix

One thing one has to know when working with **mass matrices** is best described via an example, let's consider
the same problem as before with solely changed edge and vertex:

```julia
vertex = ODEVertex(vertexfunction!, 2, [2 1; -1 1], nothing)
edge = ODEEdge(edgefunction!, 2)
```

We now have two dimensional vertex and edge variables, we additionally added a mass matrix for every vertex. The Constructor builds one
big mass matrix from all the given ones. If one now wants to solve the problem, one has to specify the solving algorithm for the solver as the
default solver can't handle mass matrices. The DAE solvers are fit for these kind of problems. One has to be especially aware of putting the variable autodiff inside the algorithm to false, hence one has to write the solver like this:

```julia
sol = solve(prob, Rodas4(autodiff = false)) # Rodas4 is just an exemplary DAE solving algorithm, there are many more.#
```

With that, everything works just fine. One has to put autodiff to false, because the structure of the lastly given equations is not of the standard form that the DAE solvers can handle just like that.
