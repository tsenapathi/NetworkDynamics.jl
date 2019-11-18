using Pkg
Pkg.activate(".")

import NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using DifferentialEquations
using Plots
using BenchmarkTools

const ND = NetworkDynamics

# Get us a graph
g = barabasi_albert(10^2,5)

# First we construct the dynamical equations using sparse matrix multiplications
# We can not expect to reach this level of performance, nor do we need to.

B = incidence_matrix(g, oriented=true)

struct kuramoto_dyn{T, T2, U}
    B::T
    B_trans::T2
    ω::U
end
function (dd::kuramoto_dyn)(dx, x, p, t)
    dx .= dd.ω .- dd.B * sin.(dd.B_trans * x)
    nothing
end

ω = randn(nv(g))

kura_network_L = ODEFunction(kuramoto_dyn(B, transpose(B), ω))

# Next we construct a version that mimics what network dynamics does.

const Idx = UnitRange{Int}

struct kuramoto_dyn_loop
    num_v::Int # Number of vertices
    num_e::Int # Number of edges
    e_int::Array{Float64, 1} # Variables living on edges
    e_idx::Array{Idx, 1} # Array of indices of variables in e_int belonging to edges
    s_idx::Array{Idx, 1} # Array of indices of variables in x belonging to source vertex of edge
    d_idx::Array{Idx, 1} # Array of indices of variables in x belonging to destination vertex of edge
    v_idx::Array{Idx, 1} # Array of indices of variables in x belonging to vertex
    e_s::Array{Array{SubArray{Float64,1,Array{Float64,1},Tuple{Idx},true},1},1}
    e_d::Array{Array{SubArray{Float64,1,Array{Float64,1},Tuple{Idx},true},1},1}
end
function (dd::kuramoto_dyn_loop)(dx, x, p, t)
    @views begin
    for i in 1:dd.num_e
        dd.e_int[dd.e_idx[i]][1] = sin(x[dd.s_idx[i]][1] - x[dd.d_idx[i]][1])
    end
    for i in 1:dd.num_v
        dx[dd.v_idx[i]][1] = p[i].ω
        for e in dd.e_s[i]
            dx[dd.v_idx[i]][1] -= e[1]
        end
        for e in dd.e_d[i]
            dx[dd.v_idx[i]][1] += e[1]
        end
    end
    end #views
    nothing
end

e_int = zeros(ne(g))
gs = ND.create_graph_structure(g, ones(Int, nv(g)), ones(Int, ne(g)), e_int)

kura_network_loop = kuramoto_dyn_loop(nv(g), ne(g), gs.e_int, gs.e_idx, gs.s_idx, gs.d_idx, gs.v_idx, gs.e_s, gs.e_d)

#Now we use the network dynamics construction:

@inline function kuramoto_edge!(e,v_s,v_d,p,t)
    e[1] = sin(v_s[1] - v_d[1])
    nothing
end

@inline function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
    dv[1] = p.ω
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end

@inline function kuramoto_vertex2!(dv, v, e_s, e_d, p, t)
    dv[1] = p.ω
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end


# The Array of ODEVertex and StaticEdge contain the dynamics and other
# information, most importantly the number of variables on each vertex/edge

odevertex = ND.ODEVertex(f! = kuramoto_vertex!, dim = 1)
staticedge = ND.StaticEdge(f! = kuramoto_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

vertex_list2 = [ND.ODEVertex(f! = kuramoto_vertex2!, dim = 1), odevertex]
append!(vertex_list2, [odevertex for i in 1:(nv(g) - 2)])

struct kura_params
    ω::Float64
end

pars = [kura_params(ω_i) for ω_i in ω]
append!(pars, [kura_params(0.) for i in 1:ne(g)])

# @code_warntype ND.network_dynamics(vertex_list,edge_list,g)

kura_network_st = ND.network_dynamics(vertex_list,edge_list,g)
kura_network_st2 = ND.network_dynamics(vertex_list2,edge_list,g)

x0 = rand(nv(g))
dx_L = similar(x0)
dx_st = similar(x0)
dx_st2 = similar(x0)
dx_loop = similar(x0)

# @code_warntype kura_network_st(dx_st, x0, nothing, 0.)
# @code_warntype kura_network_loop(dx_loop, x0, nothing, 0.)

kura_network_st(dx_st, x0, pars, 0.)
kura_network_st2(dx_st2, x0, pars, 0.)
kura_network_loop(dx_loop, x0, pars, 0.)
kura_network_L(dx_L, x0, nothing, 0.)

println("Benachmarking")

println("Network Dynamics")
@btime kura_network_st(dx_st, x0, pars, 0.)
println("Network Dynamics with two functions")
@btime kura_network_st2(dx_st, x0, pars, 0.)
println("Network Dynamics 'by hand'")
@btime kura_network_loop(dx_L, x0, pars, 0.)
println("Sparse network multiplication")
@btime kura_network_L(dx_L, x0, nothing, 0.)
