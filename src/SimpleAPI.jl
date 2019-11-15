export homogeneous_scalar_network_with_sum

# """
# Implement dx_i/dt = f(x_i, p, t) + \sum_j A_{ij} h(x_i, x_j, p, t) with scalar x
# """
function homogeneous_scalar_network_with_sum(f, h, graph)

    function e!(e, x_s, x_t, p, t)
        e[1] = h(x_s[1], x_t[1], p, t)
        e[2] = h(x_t[1], x_s[1], p, t)
        nothing
    end

    function v!(dv, v, e_s, e_t, p, t)
        dv[1] = f(v[1], p, t) + sum( [e[1] for e in e_s] ) + sum( [e[2] for e in e_t] )
        nothing
    end

    odevertex = ODEVertex(f! = v!, dim = 1)
    staticedge = StaticEdge(f! = e!, dim = 2)

    vertex_list = [odevertex for v in vertices(graph)]
    edge_list = [staticedge for e in edges(graph)]

    network_dynamics(vertex_list, edge_list, graph)
end
