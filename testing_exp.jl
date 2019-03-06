begin
    using LightGraphs
    using LinearAlgebra
    using DifferentialEquations
    using Plots
end

include("src/NetworkDynamics_exp.jl")

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
    G =barabasi_albert(30, 4)

    edges! = [NetworkDynamics_exp.diffusion_edge! for e in edges(G)]
    vertices! = [NetworkDynamics_exp.diffusion_vertex! for e in vertices(G)]

    dim_v = 2 * ones(Int32, nv(G))
    dim_e = 2 * ones(Int32, ne(G))

    rhs = NetworkDynamics_exp.multi_static(vertices!, edges!, G, dim_v, dim_e)

    ic = rand(sum(dim_v))

    prob = ODEProblem(rhs,ic,(0.,50.))
end

sol = solve(prob)

begin
    x_f = sol[:, end]
    x_f1 = sum(ic[1:2:end])/30
    x_f2 = sum(ic[2:2:end])/30

    @assert isapprox(sum(ic[1:2:end]), sum(x_f[1:2:end]))
    @assert isapprox(sum(ic[2:2:end]), sum(x_f[2:2:end]))

    @assert isapprox(x_f1, x_f[1])
    @assert isapprox(x_f2, x_f[2])
end

plot(sol,legend=false)
