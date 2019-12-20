
# New design, this is a self contained file defining everything for the new
# design of ND.jl

module ND2

using LightGraphs
using DiffEqOperators
using DiffEqBase
using LinearAlgebra


@Base.kwdef struct StaticEdge{T}
    f!::T # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct ODEVertex{T}
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

# We need rather complicated sets of indices into the arrays that hold the
# vertex and the edge variables. We precompute everything we can and store it
# in GraphStruct.

const Idx = UnitRange{Int}

function create_idxs(dims; counter=1)::Array{Idx, 1}
    idxs = [1:1 for dim in dims]
    for (i, dim) in enumerate(dims)
        idxs[i] = counter:(counter + dim - 1)
        counter += dim
    end
    idxs
end

function create_offsets(dims; counter=0)::Array{Int, 1}
    offs = [1 for dim in dims]
    for (i, dim) in enumerate(dims)
        offs[i] = counter
        counter += dim
    end
    offs
end

function create_idxs(offs, dims)::Array{Idx, 1}
    idxs = [1+off:off+dim for (off, dim) in zip(offs, dims)]
end


struct GraphStruct
    num_v::Int
    num_e::Int
    v_dims::Array{Int, 1}
    e_dims::Array{Int, 1}
    s_e::Array{Int, 1}
    d_e::Array{Int, 1}
    v_offs::Array{Int, 1}
    e_offs::Array{Int, 1}
    v_idx::Array{Idx, 1}
    e_idx::Array{Idx, 1}
    s_e_offs::Array{Int, 1}
    d_e_offs::Array{Int, 1}
    s_e_idx::Array{Idx, 1}
    d_e_idx::Array{Idx, 1}
    e_s_v_dat::Array{Array{Tuple{Int,Int}, 1}}
    e_d_v_dat::Array{Array{Tuple{Int,Int}, 1}}
end
function GraphStruct(g, v_dims, e_dims)
    num_v = nv(g)
    num_e = ne(g)

    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]

    v_offs = create_offsets(v_dims)
    e_offs = create_offsets(e_dims)

    v_idx = create_idxs(v_offs, v_dims)
    e_idx = create_idxs(e_offs, e_dims)

    s_e_offs = [v_offs[s_e[i_e]] for i_e in 1:num_e]
    d_e_offs = [v_offs[d_e[i_e]] for i_e in 1:num_e]

    s_e_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_e_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]

    e_s_v_dat = [[(offset, dim) for (i_e, (offset, dim)) in enumerate(zip(e_offs, e_dims)) if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d_v_dat = [[(offset, dim) for (i_e, (offset, dim)) in enumerate(zip(e_offs, e_dims)) if i_v == d_e[i_e]] for i_v in 1:num_v]

    GraphStruct(
    num_v,
    num_e,
    v_dims,
    e_dims,
    s_e,
    d_e,
    v_offs,
    e_offs,
    v_idx,
    e_idx,
    s_e_offs,
    d_e_offs,
    s_e_idx,
    d_e_idx,
    e_s_v_dat,
    e_d_v_dat)
end

# In order to access the data in the arrays efficiently we create views that
# allow us to efficiently index into the underlying arrays.

import Base.getindex, Base.setindex!, Base.length

struct ConstArray{T}
    ca_val::T
end

@inline function getindex(ca::ConstArray, idx)
    ca.ca_val
end


struct EdgeData{G}
    gd::G
    idx_offset::Int
    len::Int
end

@inline function getindex(e_dat::EdgeData, idx)
    e_dat.gd.e_array[idx + e_dat.idx_offset]
end

@inline function setindex!(e_dat::EdgeData, x, idx)
    e_dat.gd.e_array[idx + e_dat.idx_offset] = x
    nothing
end

@inline function length(e_dat::EdgeData)
    e_dat.len
end



struct VertexData{G}
    gd::G
    idx_offset::Int
    len::Int
end

@inline function getindex(v_dat::VertexData, idx)
    v_dat.gd.v_array[idx + v_dat.idx_offset]
end

@inline function setindex!(v_dat::VertexData, x, idx)
    v_dat.gd.v_array[idx + v_dat.idx_offset] = x
    nothing
end

@inline function length(v_dat::VertexData)
    v_dat.len
end

# Putting the above together we create a GraphData object:

# An alternative design that needs to be evaluated for performance is to create
# only one array of VertexData and EdgeData and index into that, possibly with a
# new set of access types...

mutable struct GraphData{T}
    v_array::T
    e_array::T
    v::Array{VertexData{GraphData{T}}, 1}
    e::Array{EdgeData{GraphData{T}}, 1}
    v_s_e::Array{VertexData{GraphData{T}}, 1} # the vertex that is the source of e
    v_d_e::Array{VertexData{GraphData{T}}, 1} # the vertex that is the destination of e
    e_s_v::Array{Array{EdgeData{GraphData{T}}, 1}, 1} # the edges that have v as source
    e_d_v::Array{Array{EdgeData{GraphData{T}}, 1}, 1} # the edges that have v as destination
    function GraphData{T}(v_array::T, e_array::T, gs::GraphStruct) where T
        gd = new{T}(v_array, e_array)
        gd.v = [VertexData{GraphData{T}}(gd, offset, dim) for (offset,dim) in zip(gs.v_offs, gs.v_dims)]
        gd.e = [EdgeData{GraphData{T}}(gd, offset, dim) for (offset,dim) in zip(gs.e_offs, gs.e_dims)]
        gd.v_s_e = [VertexData{GraphData{T}}(gd, offset, dim) for (offset,dim) in zip(gs.s_e_offs, gs.v_dims[gs.s_e])]
        gd.v_d_e = [VertexData{GraphData{T}}(gd, offset, dim) for (offset,dim) in zip(gs.d_e_offs, gs.v_dims[gs.d_e])]
        gd.e_s_v = [[EdgeData{GraphData{T}}(gd, offset, dim) for (offset,dim) in e_s_v] for e_s_v in gs.e_s_v_dat]
        gd.e_d_v = [[EdgeData{GraphData{T}}(gd, offset, dim) for (offset,dim) in e_d_v] for e_d_v in gs.e_d_v_dat]
        gd
    end
end

function GraphData(v_array, e_array, gs)
    GraphData{typeof(v_array)}(v_array, e_array, gs)
end

function construct_mass_matrix(mmv_array, dim_nd, gs::GraphStruct)
    if all([mm == I for mm in mmv_array])
        mass_matrix = I
    else
        mass_matrix = sparse(1.0I,dim_nd,dim_nd)
        for (i, mm) in enumerate(mmv_array)
            if mm != I
                mass_matrix[gs.v_idx[i],gs.v_idx[i]] .= mm
            end
        end
    end
    mass_matrix
end


@Base.kwdef struct nd_ODE_Static_2{G,T,TF1,TF2}
    tup_vertices!::TF1
    v_type_idx::Array{Array{Int, 1}, 1}
    tup_edges!::TF2
    e_type_idx::Array{Array{Int, 1}, 1}
    graph::G
    graph_structure::GraphStruct
    graph_data::GraphData{T}
    dgraph_data::GraphData{T}
end

function (d::nd_ODE_Static_2{G,T,TF1,TF2})(dx::T, x::T, p, t) where G where T where TF1 where TF2
    # print("Type stable version")
    gd = d.graph_data
    gd.v_array = x

    dgd = d.dgraph_data
    dgd.v_array = dx

    for (i_e, edges!) in enumerate(d.tup_edges!)
    for (idx_e, edge) in enumerate(edges!) # This is a typestable loop with the same f! type throughout
        i = d.e_type_idx[i_e][idx_e]
        edge.f!(gd.e[i], gd.v_s_e[i], gd.v_d_e[i], p[i+d.graph_structure.num_v], t)
    end
    end

    for (i_v, vertices!) in enumerate(d.tup_vertices!)
    for (idx_v, vertex) in enumerate(vertices!)
        i = d.v_type_idx[i_v][idx_v]
        vertex.f!(dgd.v[i], gd.v[i], gd.e_s_v[i], gd.e_d_v[i], p[i], t)
    end
    end

    nothing
end

function (d::nd_ODE_Static_2)(dx, x, p, t)
    # print("Type unstable version")

    e_array = similar(x, d.graph_structure.num_e)
    gd = GraphData(x, e_array, d.graph_structure)

    # de_array = similar(dx, d.graph_structure.num_e)
    # dgd = GraphData(dx, de_array, d.graph_structure)

    for i in 1:d.graph_structure.num_e
        d.edges![i].f!(gd.e[i], gd.v_s_e[i], gd.v_d_e[i], p[i+d.graph_structure.num_v], t)
    end

    for i in 1:d.graph_structure.num_v
        # d.vertices![i].f!(dgd.v[i], gd.v[i], gd.e_s_v[i], gd.e_d_v[i], p[i], t)
        d.vertices![i].f!(view(dx,d.graph_structure.v_idx[i]), gd.v[i], gd.e_s_v[i], gd.e_d_v[i], p[i], t)
    end

    nothing
end


function network_dynamics_2(vertices!, edges!, graph, p; x_prototype=zeros(1))
    @assert length(vertices!) == length(vertices(graph))
    @assert length(edges!) == length(edges(graph))

    v_dims = [v.dim for v in vertices!]
    e_dims = [e.dim for e in edges!]

    v_array = similar(x_prototype, sum(v_dims))
    e_array = similar(x_prototype, sum(e_dims))

    v_type_array = vertices! .|> typeof
    tup_v_types = Tuple(v_type_array |> unique)
    v_type_idx = [[i for (i, T) in enumerate(v_type_array) if T == Tv] for Tv in tup_v_types]
    v_tup = Tuple([[vertices![i] for i in v_type_i] for v_type_i in v_type_idx])

    e_type_array = edges! .|> typeof
    tup_e_types = Tuple(e_type_array |> unique)
    e_type_idx = [[i for (i, T) in enumerate(e_type_array) if T == Te] for Te in tup_e_types]
    e_tup = Tuple([[edges![i] for i in e_type_i] for e_type_i in e_type_idx])

    graph_stucture = GraphStruct(graph, v_dims, e_dims)

    graph_data = GraphData{typeof(v_array)}(v_array, e_array, graph_stucture)
    dgraph_data = GraphData{typeof(v_array)}(v_array, e_array, graph_stucture)

    nd! = nd_ODE_Static_2(v_tup, v_type_idx, e_tup, e_type_idx, graph, graph_stucture, graph_data, dgraph_data)

    Jv = JacVecOperator(nd!, v_array, p, 0.0)

    # Construct mass matrix
    mass_matrix = construct_mass_matrix([v.mass_matrix for v in vertices!], sum(v_dims), graph_stucture)

    symbols = [Symbol(vertices![i].sym[j],"_",i) for i in 1:length(vertices!) for j in 1:v_dims[i]]

    ODEFunction(nd!; jac_prototype=Jv, mass_matrix = mass_matrix, syms=symbols)
end


struct SolutionWrapper
    sol
    graph_data
end
function (sw::SolutionWrapper)(t)
    sw.graph_data.v_array = sw.sol(t)
    sw.graph_data
end
function (sw::SolutionWrapper)(s::Symbol, t)
    sw.graph_data.v_array = sw.sol(t)
    if   s == :v
        return sw.graph_data.v
    elseif s == :e
        return sw.graph_data.e
    else
        return nothing
    end
end
function (sw::SolutionWrapper)(s::Symbol, i, t)
    sw.graph_data.v_array = sw.sol(t)
    if   s == :v
        return sw.graph_data.v[i]
    elseif s == :e
        return sw.graph_data.e[i]
    else
        return nothing
    end
end

end
