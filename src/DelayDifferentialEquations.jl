begin
    using LightGraphs
    using LinearAlgebra
    using DifferentialEquations
    using Plots
    using Parameters
end

h(p, t) = ones(3)

 p0 = 0.2;  q0 = 0.3;  v0 = 1;  d0 = 5
 p1 = 0.2;  q1 = 0.3;  v1 = 1;  d1 = 1
 d2 = 1;  beta0 = 1;  beta1 = 1;  tau = 1
function bc_model(du,u,h,p,t)
  du[1] = (v0/(1+beta0*(h(p, t-tau)[3]^2))) * (p0 - q0)*u[1] - d0*u[1]
  du[2] = (v0/(1+beta0*(h(p, t-tau)[3]^2))) * (1 - p0 + q0)*u[1] +
          (v1/(1+beta1*(h(p, t-tau)[3]^2))) * (p1 - q1)*u[2] - d1*u[2]
  du[3] = (v1/(1+beta1*(h(p, t-tau;idxs=3)^2))) * (1 - p1 + q1)*u[2] - d2*u[3]
end

lags=[tau]
tspan = (0.0,10.0)
u0 = [1.0,1.0,1.0]
prob = DDEProblem(bc_model,u0,h,tspan)

alg = MethodOfSteps(Tsit5())
sol=solve(prob,alg)

plot(sol)


@with_kw struct static_line_delayed
    nodes!
    lines!
    s_e
    t_e
    l_e
    l_e_int
    l_s
    l_t
    len_l
    len_n
end

function (d::static_line_delayed)(dx, x, h, p, t)
    for i in 1:d.len_l
        d.lines![i](d.l_e[i],x[d.s_e[i]], x[d.t_e[i]], indexed_h(h,d.s_e[i]), indexed_h(h,d.t_e[i]), p, t)
    end
    for i in 1:d.len_n
        d.nodes![i](view(dx,i),x[i],indexed_h(h,i), sum.(d.l_s[i]), sum.(d.l_t[i]), p, t)
    end
    nothing
end

struct indexed_h
    h
    idxs
end
function (ih::indexed_h)(args...)
    ih.h(args...; idxs = ih.idxs)
end

# h_ss[i] would be the data of the ith source/destination edge...


function static_line_delayed(nodes!, lines!, s_e, t_e)
    len_l = length(lines!)
    len_n = length(nodes!)

    # Will get longer once we have more variables per line
    l_e_int = rand(len_l)

    # This will be views to more than one variable eventually
    l_e = [
        view(l_e_int, i)
        for i in 1:len_l ]

    # Create an array of views into the lines that selects only the sources
    # for each node
    l_s = [
         [l_e[j] for j in 1:len_l if s_e[j] == i]
        for i in 1:len_n ]
    # Create an array of views into the lines that selects only the targets
    # for each node
    l_t = [
        [l_e[j] for j in 1:len_l if t_e[j] == i]
        for i in 1:len_n ]

    static_line_delayed(nodes!, lines!, s_e, t_e, l_e, l_e_int, l_s, l_t, len_l, len_n)
end

function static_line_delayed(nodes!, lines!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    static_line_delayed(nodes!, lines!, s_e, t_e)
end

g=barabasi_albert(4,2)
tau = 1
taus = rand(4)
nodes! = [(dx,x,h,l_s,l_t, p,t) -> dx .=   sum(l_t) .- sum(l_s)  for v in vertices(g)]
lines! = [(l,x_s,x_t,h_s,h_t,p,t) -> l .= x_s .- h_t(p,1.01t) for e in edges(g)]
test= static_line_delayed(nodes!,lines!,g)

x0 = [1.0,2.0,3.0,4.0]
h(p,t;kwargs...) = [1.0,2.0,3.0,4.0]

testprob = DDEProblem(test,x0,h,(0.,2.))
sol =solve(testprob)

plot(sol)

function delaytest(du,u,h,p,t)
    @. du = -u
#    du[1] += h(p,t - 2.)[2]
    du[[1,2]] += h(p,t - 2.;idxs=[1,2])[1:2]
#    du[1] += indexed_h(h,2)(p,t - 2.)[1]
    nothing
end

testprob2 = DDEProblem(delaytest,x0,h,(0.,10.),constant_lags=[2.])
sol = solve(testprob2)
plot(sol)

@with_kw struct dynamic_lines_delayed
    edges!
    vertices!
    s_e
    d_e
    num_v # Number of vertices
    num_e # Number of edges
    e_int # Variables living on edges
    e_idx # Array of Array of indices of variables in e_int belonging to edges
    e_x_idx # Array of Array of indices of variables in x belonging to edges
    s_idx # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx # Array of Array of indices of variables in x belonging to vertex
    e_s # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
    tau_s
    tau_d
    no_parameters
end


function e_s_delayed(h,i,s_e,num_v,num_e,e_x_idx)
    [indexed_h(h,e_x_idx[i_e])  for i_e in 1:num_e  if i == s_e[i_e]]
end

function e_d_delayed(h,i,d_e,num_v,num_e,e_x_idx)
    [indexed_h(h,e_x_idx[i_e]) for i_e in 1:num_e if i == d_e[i_e]]
end


function (d::dynamic_lines_delayed)(dx, x, h, p, t)
    @views begin
    if d.no_parameters == true
        for i in 1:d.num_e
            d.e_int[d.e_idx[i]] .= x[d.e_x_idx[i]]
            d.edges![i](dx[d.e_x_idx[i]], x[d.e_x_idx[i]], indexed_h(h,d.e_idx[i]), x[d.s_idx[i]], x[d.d_idx[i]],indexed_h(h,d.s_idx[i]),indexed_h(h,d.d_idx[i]), p, t)
        end
        for i in 1:d.num_v
            d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], indexed_h(h,d.v_idx[i][1]), d.e_s[i], d.e_d[i], e_s_delayed(h,i,d.s_e,d.num_v,d.num_e,d.e_x_idx), e_d_delayed(h,i,d.d_e,d.num_v,d.num_e,d.e_x_idx) , p, t, d.tau_s[i],d.tau_d[i])
        end
    else
        for i in 1:d.num_e
            d.e_int[d.e_idx[i]] .= x[d.e_x_idx[i]]
            d.edges![i](dx[d.e_x_idx[i]],d.e_int[d.e_idx[i]], x[d.s_idx[i]], x[d.d_idx[i]], p[d.num_v +i], t)
        end
        for i in 1:d.num_v
            d.vertices![i](dx[d.v_idx[i]], x[d.v_idx[i]], d.e_s[i], d.e_d[i], p[i], t)
        end
    end
    end
    nothing
end

"""
dim_v is an array of the number of variables per vertex
dim_e is an array of the number of variables per edge
"""
function dynamic_lines_delayed(vertices!, edges!, s_e, d_e, dim_v, dim_e, tau_s, tau_d, no_parameters)
    num_v = length(dim_v)
    num_e = length(dim_e)

    e_int = zeros(sum(dim_e))
    e_int_2 = zeros(sum(dim_e))

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

    e_x_idx = [e_idx[i] .+ sum(dim_v) for i in 1:num_e]
    # For every vertex, and for every edge, if the source of the edge is that vertex, take the view on the variables of that edge.
    # Thus e_s[i] is an array of views onto the variables of the edges for which i is the source.
    e_s = [[view(e_int, e_idx[i_e]) for i_e in 1:num_e if i_v == s_e[i_e]] for i_v in 1:num_v]
    e_d = [[view(e_int_2, e_idx[i_e]) for i_e in 1:num_e if i_v == d_e[i_e]] for i_v in 1:num_v]

    s_idx = [v_idx[s_e[i_e]] for i_e in 1:num_e]
    d_idx = [v_idx[d_e[i_e]] for i_e in 1:num_e]

    dynamic_lines_delayed(
    edges!,
    vertices!,
    s_e,
    d_e,
    num_v, # Number of vertices
    num_e, # Number of edges
    e_int, # Variables living on edges
    e_idx, # Array of Array of indices of variables in e_int belonging to edges
    e_x_idx, #Array of Array of indices of variables in x belonging to edges
    s_idx, # Array of Array of indices of variables in x belonging to source vertex of edge
    d_idx, # Array of Array of indices of variables in x belonging to destination vertex of edge
    v_idx, # Array of Array of indices of variables in x belonging to vertex
    e_s, # Array of Array of views on the variables in e_int of the edges that are source of a vertex
    e_d, # Array of Array of views on the variables in e_int of the edges that are destination of a vertex
    tau_s,
    tau_d,
    no_parameters)
end

function dynamic_lines_delayed(vertices!, edges!, g::AbstractGraph, dim_v, dim_e, tau_s, tau_d; no_parameters = true)
    s_e = [src(e) for e in edges(g)]
    d_e = [dst(e) for e in edges(g)]
    dynamic_lines_delayed(vertices!, edges!, s_e, d_e, dim_v, dim_e, tau_s, tau_d, no_parameters)
end

dim_v = ones(Int32, nv(g))
dim_e = ones(Int32,ne(g))

s = [src(e) for e in edges(g)]
d = [dst(e) for e in edges(g)]

tau_s = [[0.01*rand(1)[1] for i_e in 1:ne(g) if i_v == s[i_e]] for i_v in 1:nv(g)]
tau_d = [[0.01*rand(1)[1] for i_e in 1:ne(g) if i_v == d[i_e]] for i_v in 1:nv(g)]

lines! = [(dl,l,h_l,x_s,x_t,h_s,h_t,p,t) -> dl .= x_s .- x_t .- l for e in edges(g)]

function vertex!(dv,v,h,e_s,e_d,h_ss,h_ds,p,t,tau_s,tau_d)
    dv .= -h(p,t)[1]
    i = 1
    j = 1
    for h_s in h_ss
        #dv .+= h_s(p,t-tau_s[i])
        i += 1
    end
    for h_d in h_ds
        #dv .-= h_d(p,t-tau_d[i])
        j += 1
    end
    for e in e_s
        #dv .-= e
    end
    for e in e_d
        #dv .+= e
    end
end

vertices! = [vertex! for vertex in vertices(g)]

dld = dynamic_lines_delayed(vertices!,lines!,g,dim_v,dim_e,tau_s,tau_d)

x0 = rand(8)

h0 = x0

dld(x0,x0,h0,nothing,0.)

dld_prob = DDEProblem(dld,x0,h0,(0.,2.))

sol= solve(dld_prob)

plot(sol,legend = false, vars = 1:4)

A=[e_s_delayed(h0,i,dld.e_s,dld.num_v,dld.num_e,dld.e_x_idx) for i in 1:8]
