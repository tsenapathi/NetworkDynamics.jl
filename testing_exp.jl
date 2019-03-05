using LightGraphs

function test_diffusive()
    G =barabasi_albert(30, 4)

    edges! = [NetworkDynamics_exp.diffusion_edge! for e in edges(G)]
    vertices! = [NetworkDynamics_exp.diffusion_vertex! for e in vertices(G)]

    dim_v = ones(Int32, nv(G))
    dim_e = ones(Int32, ne(G))

    L = laplacian_matrix(G)

    rhs = NetworkDynamics_exp.multi_static(vertices!, edges!, G, dim_v, dim_e)
    rhs2 = (dx, x, p, t) -> dx .= - L * x

    for _ in 1:10

        x0 = rand(sum(dim_v))
        dx0 = rand(sum(dim_v))
        dx1 = rand(sum(dim_v))

        p = nothing
        t = 0.

        rhs(dx0, x0, p, t)
        rhs2(dx1, x0, p, t)

        @assert isapprox(dx0, dx1)
    end
    true
end

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
