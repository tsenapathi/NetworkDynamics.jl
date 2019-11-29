export scalar_network_with_sum

# """
# Implement dx_i/dt = f(x_i, p, t) + \sum_j A_{ij} h(x_i, x_j, p, t) with scalar x
# """

"""
This function builds the homogeneous network of the form:

`dx_i/dt = f(x_i, p, t) + \\sum_j A_{ij} h(x_i, x_j, p, t)` with scalar variable `x`
"""
function scalar_network_with_sum(f, h, graph, p)

    @inline function hsws_e!(e, x_s, x_t, p, t)
        e[1] = h(x_s[1], x_t[1], p, t)
        e[2] = h(x_t[1], x_s[1], p, t)
        nothing
    end

    @inline function hsws_v!(dv, v, e_s, e_t, p, t)
        dv[1] = f(v[1], p, t) + sum( [e[1] for e in e_s] ) + sum( [e[2] for e in e_t] )
        nothing
    end

    odevertex = ODEVertex(f! = hsws_v!, dim = 1)
    staticedge = StaticEdge(f! = hsws_e!, dim = 2)

    network_dynamics(odevertex, staticedge, graph)
end
