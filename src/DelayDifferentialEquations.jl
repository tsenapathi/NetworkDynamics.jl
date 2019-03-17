begin
    using LightGraphs
    using LinearAlgebra
    using DifferentialEquations
    using Plots
    using Parameters
end

h(p, t) = ones(3)

const p0 = 0.2; const q0 = 0.3; const v0 = 1; const d0 = 5
const p1 = 0.2; const q1 = 0.3; const v1 = 1; const d1 = 1
const d2 = 1; const beta0 = 1; const beta1 = 1; const tau = 1
function bc_model(du,u,h,p,t)
  du[1] = (v0/(1+beta0*(h(p, t-tau)[3]^2))) * (p0 - q0)*u[1] - d0*u[1]
  du[2] = (v0/(1+beta0*(h(p, t-tau)[3]^2))) * (1 - p0 + q0)*u[1] +
          (v1/(1+beta1*(h(p, t-tau)[3]^2))) * (p1 - q1)*u[2] - d1*u[2]
  du[3] = (v1/(1+beta1*(h(p, t-tau)[3]^2))) * (1 - p1 + q1)*u[2] - d2*u[3]
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
        d.lines![i](d.l_e[i],x[d.s_e[i]], x[d.t_e[i]], h(p,t-1)[d.s_e[i]], h(p,t-1)[d.t_e[i]], p, t)
    end
    for i in 1:d.len_n
        d.nodes![i](view(dx,i),x[i], h(p,t-1)[i], sum.(d.l_s[i]), sum.(d.l_t[i]), p, t)
    end
    nothing
end

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

"""
When called with a graph, we construct the source and target vectors.
"""
function static_line_delayed(nodes!, lines!, g::AbstractGraph)
    s_e = [src(e) for e in edges(g)]
    t_e = [dst(e) for e in edges(g)]
    static_line_delayed(nodes!, lines!, s_e, t_e)
end

g=barabasi_albert(4,2)
tau = 1
nodes! = [(dx,x,h,l_s,l_t,p,t) -> dx .= h + sum(sum.(l_t)) - sum(sum.(l_s)) for v in vertices(g)]
lines! = [(l,x_s,x_t,h_s,h_t,p,t) -> l .= x_s .- x_t for e in edges(g)]

test= static_line_delayed(nodes!,lines!,g)

x0 = [1.0,2.0,3.0,4.0]
h(p,t) = [1.0,2.0,3.0,4.0]

testprob = DDEProblem(test,x0,h,(0.,2.))
sol =solve(testprob)

plot(sol)
