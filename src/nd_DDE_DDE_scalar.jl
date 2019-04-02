module nd_DDE_DDE_scalar_mod

using Parameters
using LightGraphs
using LinearAlgebra

export nd_DDE_DDE_scalar

@with_kw struct nd_DDE_DDE_scalar
    vertices!
    edges!
    s_e
    t_e
    e
    e_int
    e_s
    e_d
    num_e
    num_v
    tau_s
    tau_d
end

struct indexed_h
    h
    idxs
end

function (ih::indexed_h)(args...)
    ih.h(args...; idxs = ih.idxs)
end

function e_s_delayed(h, p, t, tau_s, s_e, i, num_v, num_e)
    [h(p,t-tau_s[i])[num_v + j] for j in 1:num_e if s_e[j] == i]
end

function e_d_delayed(h, p, t, tau_d, t_e, i, num_v, num_e)
    [h(p,t-tau_d[i])[num_v + j] for j in 1:num_e if t_e[j] == i]
end

function (d::nd_DDE_DDE_scalar)(dx, x, h, p, t)
    for i in 1:d.num_e
        d.e[i] .= x[d.num_v+i]
        d.edges![i](view(dx,d.num_v+i), x[d.num_v+i], indexed_h(h,d.num_v +i), x[d.s_e[i]], x[d.t_e[i]], indexed_h(h,d.s_e[i]), indexed_h(h,d.t_e[i]), p[d.num_v + i], t)
    end
    for i in 1:d.num_v
        d.vertices![i](view(dx,i), x[i], indexed_h(h,i), d.e_s[i], d.e_d[i], e_s_delayed(h, p, t, d.tau_s, d.s_e, i, d.num_v, d.num_e), e_d_delayed(h, p, t, d.tau_d, d.t_e, i, d.num_v, d.num_e), p[i], t)
    end
    nothing
end

function (d::nd_DDE_DDE_scalar)(dx, x, h, p::Nothing, t)
    for i in 1:d.num_e
        d.e[i] .= x[d.num_v+i]
        d.edges![i](view(dx,d.num_v+i),x[d.num_v+i], indexed_h(h,d.num_v +i), x[d.s_e[i]], x[d.t_e[i]], indexed_h(h,d.s_e[i]), indexed_h(h,d.t_e[i]), p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](view(dx,i), x[i], indexed_h(h,i), d.e_s[i], d.e_d[i], e_s_delayed(h, p, t, d.tau_s, d.s_e, i, d.num_v, d.num_e), e_d_delayed(h, p, t, d.tau_d, d.t_e, i, d.num_v, d.num_e), p, t)
    end
    nothing
end

function nd_DDE_DDE_scalar(vertices!, edges!, s_e::Array{Int64}, t_e::Array{Int64}, tau_s, tau_d)
    num_e = length(edges!)
    num_v = length(vertices!)

    # Will get longer once we have more variables per edge
    e_int = zeros(num_e)

    # This will be views to more than one variable eventually
    e = [
        view(e_int, i)
        for i in 1:num_e ]

    # Create an array of views into the edges that selects only the sources
    # for each vertex
    e_s = [
         [e[j] for j in 1:num_e if s_e[j] == i]
        for i in 1:num_v ]
    # Create an array of views into the edges that selects only the targets
    # for each vertex
    e_d = [
        [e[j] for j in 1:num_e if t_e[j] == i]
        for i in 1:num_v ]
    nd_DDE_DDE_scalar(vertices!, edges!, s_e, t_e, e, e_int, e_s, e_d, num_e, num_v, tau_s, tau_d)
end

"""
    When called with a graph, we construct the source and target vectors."""
function nd_DDE_DDE_scalar(vertices!, edges!, g::AbstractGraph, tau_s, tau_d)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    nd_DDE_DDE_scalar(vertices!, edges!, s_e, t_e, tau_s, tau_d)
end

end # module
