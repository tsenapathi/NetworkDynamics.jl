using LightGraphs
using DifferentialEquations

function network_dynamics(Array{VertexFunction}, Array{EdgeFunction}, graph)
    nothing
end

function network_dynamics(Array{ODEVertex}, Array{StaticEdge}, graph)
    nothing
end

@doc doc"""
Write docstrings
"""
function network_dynamics(vertices!::Array{ODEVertex}, edges!::Array{ODEEdge}, graph)
    @assert length(vertices!) == length(nodes(graph))
    massmatrix = nothing # Construct Mass Matrix from vertices! and edges!
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    dim_nd = dim_e + dim_v
    massmatrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!

    if all(dim_v .== 1) && all(dim_e .== 1)
        nd! = scalar_network_dynamics(vertex_functions, ...)
    else
        nd! = full_network_dynamics(vertex_functions, dim_v, ...)
    end

    for i, idx in enumerate(nd!.e_idx)
        massmatrix[idx,idx] = edges![i].massmatrix
    end
    for i, idx in enumerate(nd!.e_idx)
        massmatrix[idx,idx] = edges![i].massmatrix
    end

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
