import NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using Plots
using BenchmarkTools

const ND = NetworkDynamics

include("/home/hellmann/git/NetworkDynamics/dev_new_design.jl")
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

kur_network_L = ODEFunction(kuramoto_dyn(B, transpose(B), ω))

#Now for NetworkDynamics

@inline function kuramoto_edge!(e,v_s,v_d,p,t)
    e[1] = sin(v_s[1] - v_d[1])
    nothing
end

@inline function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
    dv[1] = p
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end

odevertex = ND.ODEVertex(f! = kuramoto_vertex!,dim = 1)
staticedge = ND.StaticEdge(f! = kuramoto_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

kur_network_st = ND.network_dynamics(vertex_list,edge_list,g)

odevertex = ODEVertex(f! = kuramoto_vertex!,dim = 1)
staticedge = StaticEdge(f! = kuramoto_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

kur_network_st2 = network_dynamics_2(vertex_list,edge_list,g)

p = vcat(ω, zeros(ne(g)))

x0 = rand(nv(g))
dx_L = similar(x0)
dx_st = similar(x0)
dx_st2 = similar(x0)

kur_network_st(dx_st, x0, p, 0.)
kur_network_st2(dx_st2, x0, p, 0.)
kur_network_L(dx_L, x0, nothing, 0.)

println("Benchmarking")

@btime kur_network_st(dx_st, x0, p, 0.)
@btime kur_network_st2(dx_st2, x0, p, 0.)
@btime kur_network_L(dx_L, x0, nothing, 0.)
