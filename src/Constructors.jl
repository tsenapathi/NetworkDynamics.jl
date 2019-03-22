using LightGraphs
using DifferentialEquations

include("nd_ODE_ODE_scalar.jl")
using .nd_ODE_ODE_scalar_mod

include("nd_ODE_Static_scalar.jl")
using .nd_ODE_Static_scalar_mod

include("nd_ODE_ODE.jl")
using .nd_ODE_ODE_mod

include("nd_ODE_Static.jl")
using .nd_ODE_Static_mod

include("Functions.jl")
using .NDFunctions

#=
The internal constructors are named as Vertex type, edge type and scalarity:

_network_dynamics_ODE_ODE
_network_dynamics_ODE_Static
_network_dynamics_ODE_ODE_scalar
_network_dynamics_ODE_Static_scalar
=#

function network_dynamics(vertices!,  edges!, graph)
    try
        va = Array{VertexFunction}(vertices!)
    catch err
        println("Cannot convert the Vertices to an Array{VertexFunction}!")
        println(err)
        return nothing
    end

    ve = Array{EdgeFunction}(edges!)
    # Here goes the logic that tries to create a homogeneous array of Vertex and
    # Edge types...
    nothing
end

function network_dynamics(vertices!::Array{ODEVertex}, edges!::Array{StaticEdge}, graph)
    @assert length(vertices!) == length(nodes(graph))
    @assert length(edges!) == length(edges(graph))
    massmatrix = nothing # Construct Mass Matrix from vertices! and edges!
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = dim_e + dim_v
    massmatrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!

    if all(dim_v .== 1) && all(dim_e .== 1)
        nd! = nd_ODE_Static_scalar(vertex_functions, edge_functions, graph)
    else
        nd! = nd_ODE_Static(vertex_functions, edge_functions, graph, dim_v, dim_e)
    end
    #=
    for i, idx in enumerate(nd!.e_idx)
        massmatrix[idx,idx] = edges![i].massmatrix
    end
    for i, idx in enumerate(nd!.e_idx)
        massmatrix[idx,idx] = edges![i].massmatrix
    end
    =#
    # Test here if massmatrix is proporitonal to the identity and if so drop it.
    ODEFunction(nd!, massmatrix)
end

@doc doc"""
Write docstrings
"""
function network_dynamics(vertices!::Array{ODEVertex}, edges!::Array{ODEEdge}, graph)
    @assert length(vertices!) == length(nodes(graph))
    @assert length(edges!) == length(edges(graph))
    massmatrix = nothing # Construct Mass Matrix from vertices! and edges!
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = dim_e + dim_v
    massmatrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!

    if all(dim_v .== 1) && all(dim_e .== 1)
        nd! = nd_ODE_ODE_scalar(vertex_functions, edge_functions, graph)
    else
        nd! = nd_ODE_ODE(vertex_functions, edge_functions, graph, dim_v, dim_e)
    end

    for i, idx in enumerate(nd!.e_idx)
        massmatrix[idx,idx] = edges![i].massmatrix
    end
    for i, idx in enumerate(nd!.e_idx)
        massmatrix[idx,idx] = edges![i].massmatrix
    end
    # Test here if massmatrix is proporitonal to the identity and if so drop it.
    ODEFunction(nd!, massmatrix)
end

#=
Future work
=#

@doc doc"""
Write docstrings
"""
function network_dynamics(vertex!::ODEVertex, edge!::ODEEdge, graph)
    throw("unimplemented")
    @assert length(vertices!) == length(nodes(graph))
    massmatrix = nothing # Construct Mass Matrix from vertices! and edges!
    vertex_functions = [v.f! for v in vertex!]
    dim_v = [v.dim for v in vertices!]
    if all(dim_v .== 1) && all(dim_e .== 1)
        nd! = scalar_homogeneous_network_dynamics(vertex_functions, ...)
    else
        nd! = full_homogeneous_network_dynamics(vertex_functions, dim_v, ...)
    end
    ODEFunction(nd!, massmatrix)
end
