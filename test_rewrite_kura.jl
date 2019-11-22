using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using DiffEqOperators
using BenchmarkTools

N = 10^3
g = barabasi_albert(N,5)

B = incidence_matrix(g, oriented=true)

struct kuramoto_dyn{T, T2, U}
    B::T
    B_trans::T2
    ω::U
    N::Int
end
function (dd::kuramoto_dyn)(dx, x, p, t)
    dx[1:N] .= dd.ω .- x[1:N] .- dd.B * sin.(dd.B_trans * x[N+1:2N])
    dx[N+1:2N] .= x[1:N]
    nothing
end

ω = randn(nv(g))

kn = kuramoto_dyn(B, transpose(B), ω, N)
Jv = JacVecOperator(kn, randn(nv(g)), nothing, 0.0)
kur_network_L = ODEFunction(kn, jac_prototype=Jv)

#Now for NetworkDynamics

@inline Base.@propagate_inbounds function kuramoto_edge!(e,v_s,v_d,p,t)
    e[1] = sin(v_s[2] - v_d[2])
    nothing
end

# @inline Base.@propagate_inbounds

@inline Base.@propagate_inbounds function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
    dv[1] = p - v[1]
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    dv[2] = v[1]
    nothing
end


odevertex = ODEVertex(f! = kuramoto_vertex!,dim = 2)
staticedge = StaticEdge(f! = kuramoto_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]
p = vcat(ω, zeros(ne(g)))

kur_network_nd = network_dynamics(vertex_list,edge_list,g, p)


x0 = rand(2N)
dx_L = similar(x0)
dx_nd = similar(x0)

kur_network_nd(dx_nd, x0, p, 0.)
kur_network_L(dx_L, x0, nothing, 0.)

dx_nd[1:2:2N] .- dx_L[1:N] .|> abs |> maximum

println("Benchmarking")

@btime kur_network_nd(dx_nd, x0, p, 0.)
@btime kur_network_L(dx_L, x0, nothing, 0.)
