using NetworkDynamics
using LightGraphs
# using LinearAlgebra
# using SparseArrays
using OrdinaryDiffEq
# using Plots
using BenchmarkTools


g = barabasi_albert(10^4,5)
L = laplacian_matrix(g)

function diffusion_dyn(dx, x, p, t)
    dx .= - L * x
    nothing
end

diff_network_L = ODEFunction(diffusion_dyn)

x0 = rand(nv(g))

# prob_L = ODEProblem(diff_network_L,x0,(0.,5.))
# sol_L = solve(prob_L)

#Now for NetworkDynamics

@inline Base.@propagate_inbounds  function diffusion_edge!(e,v_s,v_d,p,t)
    e[1] = v_s[1] - v_d[1]
    nothing
end

@inline Base.@propagate_inbounds  function diffusion_vertex!(dv, v, e_s, e_d, p, t)
    dv[1] = 0.
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] += e[1]
    end
    nothing
end

odevertex = ODEVertex(f! = diffusion_vertex!,dim = 1)
staticedge = StaticEdge(f! = diffusion_edge!, dim = 1)

vertex_list = [odevertex for v in vertices(g)]
edge_list = [staticedge for e in edges(g)]

diff_network_st = network_dynamics(vertex_list,edge_list,g, nothing)


x0 = rand(nv(g))
dx_L = similar(x0)
dx_st = similar(x0)

diff_network_st(dx_st, x0, nothing, 0.)
diff_network_L(dx_L, x0, nothing, 0.)

dx_st .- dx_L .|> abs |> maximum

println("Benchmarking")

@btime diff_network_st(dx_st, x0, nothing, 0.)
@btime diff_network_L(dx_L, x0, nothing, 0.)
