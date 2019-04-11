using NetworkDynamics
using LightGraphs
using LinearAlgebra
using DifferentialEquations

g = barabasi_albert(10,5)

function vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0
    for e in e_s
        dv .-= e
    end
    for e in e_d
        dv .+= e
    end
    nothing
end

function ddevertex!(dv,v,h_v,e_s,e_d,h_s,h_d,p,t)
    dv .= 0
    for h in h_s
        dv .-= h
    end
    for h in h_d
        dv .+= h
    end
end

odeedge! = (dl,l,v_s,v_d,p,t) -> dl .= 1000*(v_s - v_d - l)
staticedge! = (l,v_s,v_d,p,t) -> l .= v_s - v_d
ddeedge! = (de,e,h_e,v_s,v_d,h_s,h_d,p,t) -> de .= 1000*(v_s - v_d - e)

#testing the scalar constructors

odevertices = [ODEVertex(vertex!,1) for v in vertices(g)]
odeedges = [ODEEdge(odeedge! ,1) for e in edges(g)]
staticedges = [StaticEdge(staticedge!, 1) for e in edges(g)]
ddevertices = [DDEVertex(ddevertex!, 1, 0, 0) for v in vertices(g)]
ddeedges = [DDEEdge(ddeedge!, 1) for e in edges(g)]


test1 = network_dynamics(odevertices,staticedges,g)

x0 = rand(35)
h0(p,t; idxs = 1:35) = x0

test_prob1 = ODEProblem(test1,x0,(0.,2.))

test_sol1 = solve(test_prob1)

using Plots

plot(test_sol1, legend = false, vars = 1:10)

test2 = network_dynamics(odevertices,odeedges,g)

test_prob2 = ODEProblem(test2,x0,(0.,2.))

test_sol2 = solve(test_prob2)

plot(test_sol2,legend = false, vars= 1:10)


test3 = network_dynamics(ddevertices,ddeedges,g)

test_prob3 = DDEProblem(test3,x0,h0,(0.,2.))

test_sol3 = solve(test_prob3)

plot(test_sol3,legend = false, vars = 1:10)

for i in 1:10
    @test isapprox(test_sol1[end][i], test_sol2[end][i], atol=0.001)
end

for i in 1:10
    @test isapprox(test_sol1[end][i], test_sol3[end][i], atol=0.001)
end
