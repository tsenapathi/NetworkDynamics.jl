import NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using Plots
using BenchmarkTools

const ND = NetworkDynamics

include("/home/hellmann/git/NetworkDynamics/dev_new_design.jl")
g = barabasi_albert(10^3,5)

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

p = vcat(ω, zeros(ne(g)))

kur_network_st2 = network_dynamics_2(vertex_list, edge_list, g, p; x_prototype=zeros(nv(g)))

x0 = rand(nv(g))
dx_L = similar(x0)
dx_st = similar(x0)
dx_st2 = similar(x0)
dx_st3 = similar(x0)

x1 = zeros(Float32, nv(g))
x1 .= x0

kur_network_st(dx_st, x0, p, 0.)
kur_network_L(dx_L, x0, nothing, 0.)
isapprox(dx_st, dx_L)

kur_network_st2(dx_st2, x0, p, 0.) # This calls the type stable method
kur_network_st2(dx_st3, x1, p, 0.) # This calls the non-type stable method

dx_st - dx_st2 .|> abs |> maximum
dx_st - dx_st3 .|> abs |> maximum

println("Benchmarking")

@btime kur_network_st(dx_st, x0, p, 0.)
@btime kur_network_st2(dx_st2, x0, p, 0.)
@btime kur_network_L(dx_L, x0, nothing, 0.)

prob_st = ODEProblem(kur_network_st,x0,(0.,5.),p)
prob_st2 = ODEProblem(kur_network_st2,x0,(0.,5.),p)
prob_L = ODEProblem(kur_network_L,x0,(0.,5.))

println("Benchmarking Solves")
# @time solve(prob_st, Rodas4(autodiff=false))
@time sol=solve(prob_st2, TRBDF2(linsolve=LinSolveGMRES()))
@time sol=solve(prob_L, TRBDF2(linsolve=LinSolveGMRES()))
@time solve(prob_st, Rodas4())

# @time solve(prob_st2, Rodas4(autodiff=false))
@time solve(prob_st2, Rodas4())
@time solve(prob_L, Rodas4())

using Sundials

@time solve(prob_st, CVODE_BDF())
@time solve(prob_st2, CVODE_BDF())
@time solve(prob_L, CVODE_BDF())

# Paul suggested to use GraphData objects to access the solution object more
# elegantly:

sw = SolutionWrapper(sol, kur_network_st2.f.graph_data)

b = t -> [v[1] for v in sw(:v, t)]
sw(:v, 4.)[1]
sw(:v, 1, 4.)[1]

b(4.)
