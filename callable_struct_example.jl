using Parameters
using LightGraphs
using DifferentialEquations
using Plots




function diffusion_edge!(e, v_s, v_d, p, t)
    e .= v_s .- v_d
    nothing
end

function diffusion_vertex!(dv, v, e_ss, e_ds, p, t)
    # Note that e_ss and e_ds might be empty, the code needs to be able to deal
    # with this situation.
    dv .= 0

    for e_s in e_ss
        dv .-= e_s
    end

    for e_d in e_ds
        dv .+= e_d
    end

    nothing
end

@with_kw struct multi_static
    edges!
    vertices!
    num_v # Number of vertices
    num_e # Number of edges
    e_int # Variables living on edges
    e_idx # Array of Array of indices of variables in e_int belonging to edges
    s_idx # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx # Array of Array of indices of variables in x belonging to vertex
    e_s # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function (d::multi_static)(dx, x, p, t)
    @views begin

    for i in 1:d.num_e
        # d.l_e[i] = d.l_e_int[l_idx[i]] as view
        d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p, t)
    end

    end
    nothing
end

"""
dim_v is an array of the number of variables per vertex
dim_e is an array of the number of variables per edge
"""
function multi_static(vertices!, edges!, s_e, d_e, dim_v, dim_e)
    num_v = length(dim_v)
    num_e = length(dim_e)

    e_int = zeros(sum(dim_e))

    # x has length sum(dim_v)
    # v_idx is an Array of Array of indices of variables in x belonging to vertex

    counter = 1
    v_idx = [zeros(Int32, dim) for dim in dim_v]
    for i in 1:num_v
        v_idx[i] .= collect(counter:counter + dim_v[i] - 1)
        counter += dim_v[i]
    end

    counter = 1
    e_idx = [zeros(Int32, dim) for dim in dim_e]
    for i in 1:num_e
        e_idx[i] .= collect(counter:counter + dim_e[i] - 1)
        counter += dim_e[i]
    end

    # For every vertex, and for every edge, if the source of the edge is that vertex, take the view on the variables of that edge.
    # Thus e_s[i] is an array of views onto the variables of the edges for which i is the source.
    e_s = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == d_e[i_e]] for i_v in 1:num_v]

    s_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]


    multi_static(
    edges!,
    vertices!,
    num_v, # Number of vertices
    num_e, # Number of edges
    e_int, # Variables living on edges
    e_idx, # Array of Array of indices of variables in e_int belonging to edges
    s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d) # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function multi_static(vertices!, edges!, g::AbstractGraph, dim_v, dim_e)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]
    multi_static(vertices!, edges!, s_e, d_e, dim_v, dim_e)
end

@with_kw struct multi_dyn
    edges!
    vertices!
    num_v # Number of vertices
    num_e # Number of edges
    e_int # Variables living on edges
    e_idx # Array of Array of indices of variables in e_int belonging to edges
    e_x_idx
    s_idx # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx # Array of Array of indices of variables in x belonging to vertex
    e_s # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function (d::multi_dyn)(dx, x, p, t)
    @views begin

    for i in 1:d.num_e
        # d.l_e[i] = d.l_e_int[l_idx[i]] as view
        d.edges![i](d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p, t)
    end
    for i in 1:d.num_v
        d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p, t)
    end

    end
    nothing
end

"""
dim_v is an array of the number of variables per vertex
dim_e is an array of the number of variables per edge
"""
function multi_dyn(vertices!, edges!, s_e, d_e, dim_v, dim_e)
    num_v = length(dim_v)
    num_e = length(dim_e)

    e_int = zeros(sum(dim_e))

    # x has length sum(dim_v)
    # v_idx is an Array of Array of indices of variables in x belonging to vertex

    counter = 1
    v_idx = [zeros(Int32, dim) for dim in dim_v]
    for i in 1:num_v
        v_idx[i] .= collect(counter:counter + dim_v[i] - 1)
        counter += dim_v[i]
    end

    counter = 1
    e_idx = [zeros(Int32, dim) for dim in dim_e]
    for i in 1:num_e
        e_idx[i] .= collect(counter:counter + dim_e[i] - 1)
        counter += dim_e[i]
    end

    e_x_idx = [e_idx[i] .+= sum(dim_v) for i in 1:num_e]
    # For every vertex, and for every edge, if the source of the edge is that vertex, take the view on the variables of that edge.
    # Thus e_s[i] is an array of views onto the variables of the edges for which i is the source.
    e_s = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == d_e[i_e]] for i_v in 1:num_v]

    s_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]


    multi_dyn(
    edges!,
    vertices!,
    num_v, # Number of vertices
    num_e, # Number of edges
    e_int, # Variables living on edges
    e_idx, # Array of Array of indices of variables in e_int belonging to edges
    e_x_idx,
    s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d) # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
end

function multi_dyn(vertices!, edges!, g::AbstractGraph, dim_v, dim_e)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]
    multi_dyn(vertices!, edges!, s_e, d_e, dim_v, dim_e)
end

g= barabasi_albert(10,4)
edges! = [diffusion_edge! for e in edges(g)]
vertices! = [diffusion_vertex! for v in vertices(g)]
x0 = rand(10)
dim_v = ones(Int32,10)
dim_e = ones(Int32,24)

a=multi_static(vertices!,edges!,g,dim_v,dim_e)

aprob = ODEProblem(a,x0,(0.,2.))
sol = solve(aprob)
plot(sol, legend = false)

b = multi_dyn(vertices!,edges!,g,dim_v,dim_e)

bprob = ODEProblem(b,x0,(0.,2.))
sol2 = solve(bprob)
plot(sol2, legend = false)

e_idx = [zeros(Int32, dim) for dim in dim_e]
counter = 1
e_idx = [zeros(Int32, dim) for dim in dim_e]
for i in 1:a.num_e
    e_idx[i] .= collect(counter:counter + dim_e[i] - 1)
    counter += dim_e[i]
end
[e_idx[i] .+= sum(dim_v) for i in 1:a.num_e]
