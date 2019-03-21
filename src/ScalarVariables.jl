module ScalarVariables

using Parameters
using LightGraphs
using LinearAlgebra

export diffusive_network_dynamics

"""
    diffusive_network_dynamics(L, nodes)
    L: Matrix
    nodes: scalar function ``x \\rightarrow nodes(x)``"""
@with_kw struct diffusive_network_dynamics{T}
    L::AbstractArray{T,2}
    nodes::Function
end

"""
    (dnd::diffusive_network_dynamics)(dx, x, p, t)

    Calling a struct of type diffusive network dynamics implements the ODE:

        ``\\frac{dx_i}{dt} = nodes(x_i) + \\sum_j L_{ij} x_j``

    dnd = diffusive_network_dynamics(L, nodes)
    dnd(dx, x, p, t)"""

function (dnd::diffusive_network_dynamics)(dx, x, p, t)
    mul!(dx, dnd.L, x) # dx .= L * x
    dx_temp = 0
    for i in 1:length(dx)
        dnd.nodes(dx_temp, x[i], p, t)
        dx[i] = dx_temp - x[i]
    end
    nothing
end

"""
    When called with a graph, the dynamics defaults to using the laplacian.
"""
function diffusive_network_dynamics(g::AbstractGraph, nodes)
    diffusive_network_dynamics(laplacian_matrix(g), nodes)
end

# The main, fully flexible network dynamics implementation.

"""
    The key functions or function arrays are:

    nodes: ``nodes!(dx, x, [l]_s,[l]_t, p, t)``

    lines: ``lines!(dl, l, x_s, x_t, p, t)``

    Given edges ``e``, ans nodes ``n``, as well as an orientation encoded by
    the source function ``s(e)`` and the target function ``t(e)``
    this implements the system of ODEs:

    ``\\frac{dx_n}{dt} = dx_n``

    ``\\frac{dl_e}{dt} = dl_e``

    with ``dx`` and ``dl`` calculated by

    ``[l]_s = [l_e \\text{ if } s(e) = n]``

    ``[l]_t = [l_e \\text{ if } t(e) = n]``

    ``nodes![n](dx_n, x_n, [l]_s, [l]_t, p_n, t)``

    ``lines![e](dl_e, l_e, x_{s(e)}, x_{t(e)}, p_e, t)``

    Alternative design:

    Something that relaxes to a diffusive network would for example be
    implemented by

        lines = (dl, l, x_1, x_2) -> dl .= 1000. * ((x_1 - x_2) - l)
        agg = (list_1, list_2) -> sum(list_1) - sum(list_2)"""


export scalar_static_lines

"""
Documentation!!
"""
@with_kw struct scalar_static_lines
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
    no_parameters
end


function (d::scalar_static_lines)(dx, x, p, t)
    if d.no_parameters == true
        for i in 1:d.num_e
            d.edges![i](d.e[i], x[d.s_e[i]], x[d.t_e[i]], p, t)
        end
        for i in 1:d.num_v
            d.vertices![i](view(dx,i), x[i], sum.(d.e_s[i]), sum.(d.e_t[i]), p, t)
        end
    else
        for i in 1:d.num_e
            d.edges![i](d.e[i], x[d.s_e[i]], x[d.t_e[i]], p[d.num_v + i], t)
        end
        for i in 1:d.num_v
            d.vertices![i](view(dx,i), x[i], sum.(d.e_s[i]), sum.(d.e_t[i]), p[i], t)
        end
    end
    nothing
end

function scalar_static_lines(vertices!, edges!, s_e, t_e, no_parameters)
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
    # Create an array of views into the lines that selects only the targets
    # for each node
    e_t = [
        [e[j] for j in 1:num_e if t_e[j] == i]
        for i in 1:num_v ]

    scalar_static_lines(vertices!, edges!, s_e, t_e, e, e_int, e_s, e_t, num_e, num_v,no_parameters)
end

"""
When called with a graph, we construct the source and target vectors.
"""
function scalar_static_lines(vertices!, edges!, g::AbstractGraph; no_parameters=true)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    scalar_static_lines(vertices!, edges!, s_e, t_e, no_parameters)
end


export scalar_dynamic_lines

@with_kw struct scalar_dynamic_lines
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
    no_parameters
end


function (d::scalar_dynamic_lines)(dx, x, p, t)
    if d.no_parameters == false
        for i in 1:d.num_e
            d.e[i] .= x[d.num_v+i]
            d.edges![i](view(dx,d.num_v+i),x[d.num_v+i], x[d.s_e[i]], x[d.t_e[i]], p[d.num_v + i], t)
        end
        for i in 1:d.num_v
            d.vertices![i](view(dx,i), x[i], sum.(d.e_s[i]), sum.(d.e_t[i]), p[i], t)
        end
    else
        for i in 1:d.num_e
            d.e[i] .= x[d.num_v+i]
            d.edges![i](view(dx,d.num_v+i),x[d.num_v+i], x[d.s_e[i]], x[d.t_e[i]], p, t)
        end
        for i in 1:d.num_v
            d.vertices![i](view(dx,i), x[i], sum.(d.e_s[i]), sum.(d.e_t[i]), p, t)
        end
    end
    nothing
end

function scalar_dynamic_lines(vertices!, edges!, s_e, t_e,no_parameters)
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
    scalar_dynamic_lines(vertices!, edges!, s_e, t_e, e, e_int, e_s, e_t, num_e, num_v,no_parameters)
end

"""
    When called with a graph, we construct the source and target vectors."""

function scalar_dynamic_lines(vertices!, edges!, g::AbstractGraph; no_parameters = true)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    scalar_dynamic_lines(vertices!, edges!, s_e, t_e,no_parameters)
end

end # module
