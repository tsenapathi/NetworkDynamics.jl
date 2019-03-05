begin
    include("src/NetworkDynamics_exp.jl")

    using LightGraphs
    using LinearAlgebra
    using DifferentialEquations
    using Plots

    G =barabasi_albert(30, 4)

    edges! = [NetworkDynamics_exp.diffusion_edge! for e in edges(G)]
    vertices! = [NetworkDynamics_exp.diffusion_vertex! for e in vertices(G)]

    dim_v = 2 * ones(Int32, nv(G))
    dim_e = 2 * ones(Int32, ne(G))

    rhs = NetworkDynamics_exp.multi_static(vertices!, edges!, G, dim_v, dim_e)

    prob = ODEProblem(rhs,rand(sum(dim_v)),(0.,5.))
end

sol = solve(prob)

plot(sol,legend=false)
