module Constructors



# Not in use right now. Shifted into NetworkDynamics.jl for reasons.

using LightGraphs
using DifferentialEquations
using SparseArrays
using LinearAlgebra

include("nd_ODE_ODE_scalar.jl")
using ...nd_ODE_ODE_scalar_mod

include("nd_ODE_Static_scalar.jl")
using ...nd_ODE_Static_scalar_mod

include("nd_ODE_ODE.jl")
using ...nd_ODE_ODE_mod

include("nd_ODE_Static.jl")
using ...nd_ODE_Static_mod

include("Functions.jl")
using ...NDFunctions

#=
The internal constructors are named by vertex type and edge type:

nd_ODE_Static
nd_ODE_ODE
=#

export network_dynamics

function network_dynamic(vertices!,  edges!, graph)
    try
        va! = Array{NDFunctions.VertexFunction}(vertices!)
    catch err
        println("Cannot convert the vertices to an Array{VertexFunction}!")
        println(err)
        return nothing
    end
    try
        ea! = Array{NDFunctions.EdgeFunction}(edges!)
    catch err
        println("Cannot convert the edges to an Array{EdgeFunction}!")
        println(err)
        return nothing
    end
    network_dynamic(va!,  ea!, graph)
end

function network_dynamics(vertices!::Array{ODEVertex,1}, edges!::Array{StaticEdge,1}, graph)
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_v)

    nd! = nd_ODE_Static(vertex_functions, edge_functions, graph, dim_v, dim_e)

    if all([v.mass_matrix == I for v in vertices!])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!
        for i, idx in enumerate(nd!.v_idx)
            if vertices![i].mass_matrix != I
                mass_matrix[idx,idx] = vertices![i].mass_matrix
            end
        end
    end

    ODEFunction(nd!,mass_matrix = mass_matrix)
end



function network_dynamics(vertices!::Array{ODEVertex}, edges!::Array{ODEEdge}, graph)
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))
    vertex_functions = [v.f! for v in vertices!]
    dim_v = [v.dim for v in vertices!]
    edge_functions = [e.f! for e in edges!]
    dim_e = [e.dim for e in edges!]
    dim_nd = sum(dim_e) + sum(dim_v)

    nd! = nd_ODE_ODE(vertex_functions, edge_functions, graph, dim_v, dim_e)

    if all([v.mass_matrix == I for v in vertices!])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,dim_nd,dim_nd) # Construct Mass Matrix from vertices! and edges!
        for i, idx in enumerate(nd!.v_idx)
            if vertices![i].mass_matrix != I
                mass_matrix[idx,idx] = vertices![i].mass_matrix
            end
        end
        for i, idx in enumerate(nd!.e_x_idx)
            if edges![i].mass_matrix != I
                mass_matrix[idx,idx] = edges![i].mass_matrix
            end
        end
    end

    ODEFunction(nd!,mass_matrix = mass_matrix)
end

end # module
