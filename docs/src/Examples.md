# Examples

To solve the defined system,
we further need an array with initial values x0 as well as a time span tspan in which we solve the problem:

```julia

using OrdinaryDiffEq

x0 = rand(10 + 25) #10 for the vertices and 25 for the edges
tspan = (0.,2.)

prob = ODEProblem(nd,x0,tspan)
sol = solve(prob)

using Plots
plot(sol, legend = false, vars = 1:10) # vars gives us x[1:10] in the plot
```

The Plot shows the classic diffusive behaviour.

### Mass Matrix - OUTDATED

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
