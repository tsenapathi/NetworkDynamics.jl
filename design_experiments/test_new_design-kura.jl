import NetworkDynamics
using LightGraphs
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using Plots
using BenchmarkTools
using DiffEqOperators

const ND = NetworkDynamics

g = barabasi_albert(10^3,50)

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

knL = kuramoto_dyn(B, transpose(B), ω)
Jv = JacVecOperator(knL, randn(nv(g)), nothing, 0.0)
kur_network_L = ODEFunction(knL, jac_prototype=Jv)

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

@inline function kuramoto_vertex2!(dv, v, e_s, e_d, p, t)
    dv[1] = p
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end

p = vcat(ω, zeros(ne(g)))


nd_odevertex = ND.ODEVertex(f! = kuramoto_vertex!,dim = 1)
nd_odevertex2 = ND.ODEVertex(f! = kuramoto_vertex2!,dim = 1)
nd_staticedge = ND.StaticEdge(f! = kuramoto_edge!, dim = 1)

nd_vertex_list = [nd_odevertex for v in vertices(g)]
nd_vertex_list_in = [v == 2 ? nd_odevertex : nd_odevertex2 for v in vertices(g)]
nd_edge_list = [nd_staticedge for e in edges(g)]

kur_network_nd = ND.network_dynamics(nd_vertex_list,nd_edge_list,g)
kur_network_nd_in = ND.network_dynamics(nd_vertex_list_in,nd_edge_list,g)

include("/home/hellmann/git/NetworkDynamics/dev_new_design.jl")

nd1_odevertex = ND1.ODEVertex(f! = kuramoto_vertex!,dim = 1)
nd1_odevertex2 = ND1.ODEVertex(f! = kuramoto_vertex2!,dim = 1)
nd1_staticedge = ND1.StaticEdge(f! = kuramoto_edge!, dim = 1)

nd1_vertex_list = [nd1_odevertex for v in vertices(g)]
nd1_vertex_list_in = [v == 2 ? nd1_odevertex : nd1_odevertex2 for v in vertices(g)]
nd1_edge_list = [nd1_staticedge for e in edges(g)]

kur_network_nd1 = ND1.network_dynamics_2(nd1_vertex_list, nd1_edge_list, g, p; x_prototype=zeros(nv(g)))
kur_network_nd1_in = ND1.network_dynamics_2(nd1_vertex_list_in, nd1_edge_list, g, p; x_prototype=zeros(nv(g)))

x0 = rand(nv(g))
dx_L = similar(x0)
dx_nd = similar(x0)
dx_nd1 = similar(x0)

kur_network_L(dx_L, x0, nothing, 0.)
kur_network_nd(dx_nd, x0, p, 0.)
kur_network_nd1(dx_nd1, x0, p, 0.)

println("Homogeneous accuracy:")
println(dx_L - dx_nd .|> abs |> maximum)
println(dx_L - dx_nd1 .|> abs |> maximum)

println("Inhomogeneous accuracy:")
kur_network_nd_in(dx_nd, x0, p, 0.)
kur_network_nd1_in(dx_nd1, x0, p, 0.)

println(dx_L - dx_nd .|> abs |> maximum)
println(dx_L - dx_nd1 .|> abs |> maximum)

# Wrong input type triggers the non-cache version

x1 = zeros(Float32, nv(g))
x1 .= x0

println("Homogeneous accuracy, type mismatch:")
kur_network_nd1(dx_nd1, x1, p, 0.) # This calls the type unstable method
println(dx_L - dx_nd1 .|> abs |> maximum)

println("Inhomogeneous accuracy, type mismatch:")
kur_network_nd1_in(dx_nd1, x1, p, 0.) # This calls the type unstable method
println(dx_L - dx_nd1 .|> abs |> maximum)

println("Benchmarking")

println("Sparse Matrix version:")
@btime kur_network_L(dx_L, x0, nothing, 0.)
println("Homogeneous:")
print("Old:")
@btime kur_network_nd(dx_nd, x0, p, 0.)
print("New:")
@btime kur_network_nd1(dx_nd1, x0, p, 0.)

println("Inomogeneous:")
print("Old:")
@btime kur_network_nd_in(dx_nd, x0, p, 0.)
print("New:")
@btime kur_network_nd1_in(dx_nd1, x0, p, 0.)

prob_L = ODEProblem(kur_network_L, x0, (0.,5.))
prob_nd = ODEProblem(kur_network_nd, x0, (0.,5.), p)
prob_nd1 = ODEProblem(kur_network_nd1, x0, (0.,5.), p)
#
# println("Benchmarking Solves (Rodas4)")
# println("Sparse Matrix version:")
# @time solve(prob_L, Rodas4(linsolve=LinSolveGMRES())) # This crashes currently
# println("Old ND:")
# @time solve(prob_nd, Rodas4())
# println("New ND:")
# @time solve(prob_nd1, Rodas4(linsolve=LinSolveGMRES(tol = 1e-6))) # This crashes currently
#
# using Sundials
#
# println("Benchmarking Solves (CVODE_BDF)")
# println("Sparse Matrix version:")
# @time solve(prob_L, CVODE_BDF()) # This crashes currently
# println("Old ND:")
# @time solve(prob_nd, CVODE_BDF())
# println("New ND:")
# @time solve(prob_nd1, CVODE_BDF()) # This crashes currently

println("Benchmarking Solves (TRBDF2)")
println("Sparse Matrix version:")
@btime solve(prob_L, TRBDF2(linsolve=LinSolveGMRES()))
println("Old ND:")
@btime solve(prob_nd, TRBDF2(linsolve=LinSolveGMRES()))
println("New ND:")
@btime solve(prob_nd1, TRBDF2(linsolve=LinSolveGMRES()))

using Sundials
println("Old ND, CVODE_BDF:")
@btime sol_BDF = solve(prob_nd, CVODE_BDF())


sol_L = solve(prob_L, TRBDF2(linsolve=LinSolveGMRES()))
sol_nd = solve(prob_nd, TRBDF2(linsolve=LinSolveGMRES()))
sol_nd1 = solve(prob_nd1, TRBDF2(linsolve=LinSolveGMRES()))
sol_BDF = solve(prob_nd, CVODE_BDF())

plot(sol_L, vars=collect(1:10))
plot(sol_nd, vars=collect(1:10))
plot(sol_nd1, vars=collect(1:10))
plot(sol_BDF, vars=collect(1:10))
