Pkg.activate(".")
using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
using DiffEqOperators
using BenchmarkTools

N = 100
g = random_regular_graph(N,3)

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
# Jv = JacVecOperator(kn, randn(nv(g)), nothing, 0.0)
kur_network_L = ODEFunction(kn) #, jac_prototype=Jv)

#Now for NetworkDynamics

@inline Base.@propagate_inbounds function kuramoto_edge!(e,v_s,v_d,p,t)
    e[1] = sin(v_s[2] - v_d[2])
    nothing
end

@inline Base.@propagate_inbounds function kuramoto_dedge!(de, e,v_s,v_d,p,t)
    de[1] = 100. * (sin(v_s[2] - v_d[2]) - e[1])
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


odevertex = ODEVertex(f! = kuramoto_vertex!, dim = 2, sym=[:ω, :ϕ])
staticedge = StaticEdge(f! = kuramoto_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]
ode_sd_edge_list = [ODEEdge(se) for se in edge_list]
ode_edge_list = [ODEEdge(f! = kuramoto_dedge!, dim = 1) for e in edges(g)]
p = vcat(ω, zeros(ne(g)))

kur_network_nd = network_dynamics(vertex_list,edge_list,g, p)
kur_network_eode = network_dynamics(vertex_list,ode_edge_list,g, p)


x0 = rand(2N)
dx_L = zeros(2N)
dx_nd = similar(x0)

x_ode = rand(2N+ne(g))
dx_ode = similar(x_ode)

println("A new way to interact with the state of the system is to call the
network rhs with the signature (x, p, t, GetGD). This will return a GraphData
object built on x, that allows reading and manipulating the undelrying data.")

println("gd_x0 = kur_network_nd(x0, p, 0., GetGD)")
gd_x0 = kur_network_nd(x0, p, 0., GetGD) # This provides us with the graph view of the x0 data
gd_x0.e[1]
println("x0[1] == $(x0[1])")
println("gd_x0.v[1][1] == $(gd_x0.v[1][1])")
gd_x0.v[1][1] = 5.
println("x0[1] == $(x0[1])")
println("gd_x0.v[1][1] == $(gd_x0.v[1][1])")


kur_network_nd(dx_nd, x0, p, 0.)
kur_network_eode(dx_ode, x_ode, p, 0.)
kur_network_L(dx_L, x0, nothing, 0.)

kur_network_L.f.ω .- p[1:N]

ω_idx = idx_containing(kur_network_nd, :ω)

dx_nd[ω_idx] .- dx_L[1:N] .|> abs |> maximum

dgd = kur_network_nd(dx_nd, p, 0., GetGD)
dgd.v[1]

gd = kur_network_nd(x0, p, 0., GetGD)

gd.v[1]

gd.e_s_v[1]
gd.e_d_v[1]

p[1] - gd.v[1][1] - gd.e_s_v[1][1][1] - gd.e_s_v[2][1][1] - gd.e_s_v[3][1][1]

sin.(kn.B_trans * x0[N+1:2N])
kn.B[1,:]

dx_L[1]

ω[1]
x0[N+1]

println("Benchmarking")

@btime kur_network_nd(dx_nd, x0, p, 0.)
@btime kur_network_L(dx_L, x0, nothing, 0.)
@btime kur_network_eode(dx_ode, x_ode, p, 0.)
