var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#NetworkDynamics.jl-1",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "NetworkDynamics.jl",
    "title": "Overview",
    "category": "section",
    "text": "This package implements functions for defining and studying dynamics on networks. The key construction is a callable struct compatible with the DifferentialEquations.jl calling syntax.nd = network_dynamics(nodes!, lines!, s_e, t_e)\nnd(dx, x, p, t)The key functions, or function arrays from which a network dynamics is constructed are:nodesn(dx x l_sl_t p t) linese(dl l x_s x_t p t)Given edges e, and nodes n, as well as an orientation encoded by the source function s(e) and the target function t(e) this implements the system of ODEs:fracdx_ndt = dx_nfracdl_edt = dl_ewith dx and dl calculated byl_s = l_e text if  s(e) = nl_t = l_e text if  t(e) = nnodesn(dx_n x_n l_n p_n t)linese(dl_e l_e x_s(e) x_t(e) p_e t)Something that relaxes to a diffusive network would for example be implemented bylines = (dl, l, x_1, x_2) -> dl .= 1000. * ((x_1 - x_2) - l)\nagg = (list_1, list_2) -> sum(list_1) - sum(list_2)"
},

{
    "location": "#NetworkDynamics.diffusive_network_dynamics",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.diffusive_network_dynamics",
    "category": "type",
    "text": "diffusive_network_dynamics(L, nodes)\n\nL: Matrix nodes: scalar function x rightarrow nodes(x)\n\n\n\n\n\n"
},

{
    "location": "#NetworkDynamics.diffusive_network_dynamics-NTuple{4,Any}",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.diffusive_network_dynamics",
    "category": "method",
    "text": "(dnd::diffusive_network_dynamics)(dx, x, p, t)\n\nCalling a struct of type diffusive network dynamics implements the ODE:\n\nfracdx_idt = nodes(x_i) + sum_j L_ij x_j\n\ndnd = diffusive_network_dynamics(L, nodes)\ndnd(dx, x, p, t)\n\n\n\n\n\n"
},

{
    "location": "#NetworkDynamics.diffusive_network_dynamics-Tuple{LightGraphs.AbstractGraph,Any}",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.diffusive_network_dynamics",
    "category": "method",
    "text": "When called with a graph, the dynamics defaults to using the laplacian.\n\n\n\n\n\n"
},

{
    "location": "#NetworkDynamics.network_dynamics",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.network_dynamics",
    "category": "type",
    "text": "The three key functions are:\n\nnodes: nodes(dx x l p t)\n\naggregate: agg(l_sl_t)\n\nlines: lines(dl l x_s x_t p t)\n\nGiven edges e, ans nodes n, as well as an orientation encoded by the source function s(e) and the target function t(e) this implements the system of ODEs:\n\nfracdx_ndt = dx_n\n\nfracdl_edt = dl_e\n\nwith dx and dl calculated by\n\nl_s = l_e text if  s(e) = n\n\nl_t = l_e text if  t(e) = n\n\nl_n = aggn(l_s l_t)\n\nnodesn(dx_n x_n l_n p_n t)\n\nlinese(dl_e l_e x_s(e) x_t(e) p_e t)\n\nSomething that relaxes to a diffusive network would for example be implemented by\n\nlines = (dl, l, x_1, x_2) -> dl .= 1000. * ((x_1 - x_2) - l)\nagg = (list_1, list_2) -> sum(list_1) - sum(list_2)\n\n\n\n\n\n"
},

{
    "location": "#NetworkDynamics.network_dynamics_on_the_line",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.network_dynamics_on_the_line",
    "category": "type",
    "text": "edge_list: list of edges nodes: scalar function x rightarrow nodes(x) lines: scalar function x rightarrow lines(x)\n\nImplements the ODE fracdx_idt = nodes(x_i) + sum_j A_ij lines(x_i - x_j)\n\nwhere A is the adjacency matrix of B\n\n\n\n\n\n"
},

{
    "location": "#NetworkDynamics.network_dynamics_on_lines-Tuple{LightGraphs.AbstractGraph,Any}",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.network_dynamics_on_lines",
    "category": "method",
    "text": "When called with a graph, the dynamics defaults to using the laplacian.\n\n\n\n\n\n"
},

{
    "location": "#API-1",
    "page": "NetworkDynamics.jl",
    "title": "API",
    "category": "section",
    "text": "Modules = [NetworkDynamics]"
},

]}
