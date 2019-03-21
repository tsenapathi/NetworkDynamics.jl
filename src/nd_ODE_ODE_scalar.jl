module nd_ODE_ODE_scalar_mod

using Parameters
using LightGraphs
using LinearAlgebra

export nd_ODE_ODE_scalar

@with_kw struct nd_ODE_ODE_scalar
    vertices!
    edges!
    s_e
    t_e
    e
    e_int
    e_s
    e_t
    num_e
    num_v
end


function (d::nd_ODE_ODE_scalar)(dx, x, p, t)
    for i in 1:d.num_e
        d.e[i] .= x[d.num_v+i]
        d.edges![i](view(dx,d.num_v+i),x[d.num_v+i], x[d.s_e[i]], x[d.t_e[i]], p[d.num_v + i], t)
    end
    for i in 1:d.num_v
        d.vertices![i](view(dx,i), x[i], sum.(d.e_s[i]), sum.(d.e_t[i]), p[i], t)
    end
    nothing
end

function (d::nd_ODE_ODE_scalar)(dx, x, p::Nothing, t)
    for i in 1:d.num_e
        d.e[i] .= x[d.num_v+i]
        d.edges![i](view(dx,d.num_v+i),x[d.num_v+i], x[d.s_e[i]], x[d.t_e[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](view(dx,i), x[i], sum.(d.e_s[i]), sum.(d.e_t[i]), p, t)
    end
    nothing
end

function nd_ODE_ODE_scalar(vertices!, edges!, s_e, t_e)
    num_e = length(edges!)
    num_v = length(vertices!)

    # Will get longer once we have more variables per edge
    e_int = rand(num_e)

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
    e_t = [
        [e[j] for j in 1:num_e if t_e[j] == i]
        for i in 1:num_v ]
    nd_ODE_ODE_scalar(vertices!, edges!, s_e, t_e, e, e_int, e_s, e_t, num_e, num_v)
end

"""
    When called with a graph, we construct the source and target vectors."""

function nd_ODE_ODE_scalar(vertices!, edges!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    nd_ODE_ODE_scalar(vertices!, edges!, s_e, t_e)
end

end # module
