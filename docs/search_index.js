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
    "text": "This package implements functions for defining and studying dynamics on networks. The key construction is a callable struct compatible with the DifferentialEquations.jl calling syntax.nd = network_dynamics(nodes!, lines!, s_e, t_e, dim_n, dim_l; symbols_n=nothing, symbols_l=nothing)\nnd(dx, x, p, t)The first two parameters are the functions, or function arrays from which a network dynamics is constructed:nodesn(dx x l_s l_t p t) linese(dl l x_s x_t p t)The arrays dimn and diml encode the dimensionality of x and l variables. The arrays se and te encode the network structure by giving the source and target of each edge.Optionally we can also specify an array of symbols per edge and node that allow convenience access to the different nodal and line dimensions.Given edges e, and nodes n, as well as an orientation encoded by the source function s(e) and the target function t(e) this implements the system of ODEs:fracdx_ndt = dx_nfracdl_edt = dl_ewith dx and dl calculated byl_s = l_e text if  s(e) = nl_t = l_e text if  t(e) = nnodesn(dx_n x_n l_s l_t p_n t)linese(dl_e l_e x_s(e) x_t(e) p_e t)Something that relaxes to a diffusive network would for example be implemented bylines = (dl, l, x_1, x_2) -> dl .= 1000. * ((x_1 - x_2) - l)\nnodes = (dx_n, x_n, l_s, l_t, p_n, t) -> dx_n .= f(x_n) - (sum(l_s) - sum(l_t))The package also supplies a node and a line type that combine the the node function, the dimensionality and the symbols, as well as a constructor for the network dynamics from these types. Further Constructors are provided for LightGraphs:nd=network_dynamics(nodes::Array{ODE_Node}, lines::Array{ODE_Line}, G::AbstractGraph)"
},

{
    "location": "#Static-lines-1",
    "page": "NetworkDynamics.jl",
    "title": "Static lines",
    "category": "section",
    "text": "For static line relations we similarly have:sl_nd = static_lines_network_dynamics(nodes!, lines!, s_e, t_e, ...)\nsl_nd(dx, x, p, t)With the convention for lines given by:linese(l x_s x_t p t)and otherwise as above. This implements the system of ODEs:fracdx_ndt = dx_nwith dx calculated bylinese(l_e x_s(e) x_t(e) p_e t)l_s = l_e text if  s(e) = nl_t = l_e text if  t(e) = nnodesn(dx_n x_n l_s l_t p_n t)A diffusive network would be implemented bylines = (l, x_1, x_2) -> l .= x_1 - x_2\nnodes = (dx_n, x_n, l_s, l_t, p_n, t) -> dx_n .= f(x_n) - (sum(l_s) - sum(l_t))The alternative constructor is given by:sl_nd=static_lines_network_dynamics(nodes::Array{ODE_Node}, lines::Array{Static_Line}, G::AbstractGraph)"
},

{
    "location": "#Network-DAEs-1",
    "page": "NetworkDynamics.jl",
    "title": "Network DAEs",
    "category": "section",
    "text": "Design question: Don\'t do implicit DAEs, support mass matrices everywhere by default. This adds two more (optional) arrays of vectors to the signature of the constructor and the ODENode and the StaticLine types.The advantage is that we can then deal with only three types of dynamics, ODE, SDE and DDE. Leaning towards yes on this. The important special case of all ODE can then be done via a performance optimization.This would suggest splitting this into three sub packages, ODE, SDE and DDE. Each sub package should define its own XDENode and XDEline as well as promotion rules from other XDENode and XDEline types."
},

{
    "location": "#Network-SDEs-1",
    "page": "NetworkDynamics.jl",
    "title": "Network SDEs",
    "category": "section",
    "text": ""
},

{
    "location": "#Network-DDEs-1",
    "page": "NetworkDynamics.jl",
    "title": "Network DDEs",
    "category": "section",
    "text": ""
},

{
    "location": "#Convenience-functions-for-symbolic-access-to-node-variables-1",
    "page": "NetworkDynamics.jl",
    "title": "Convenience functions for symbolic access to node variables",
    "category": "section",
    "text": ""
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
    "location": "#NetworkDynamics.diffusive_network_dynamics-Tuple{AbstractGraph,Any}",
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
    "text": "The key functions or function arrays are:\n\nnodes: nodes(dx x l_sl_t p t)\n\nlines: lines(dl l x_s x_t p t)\n\nGiven edges e, ans nodes n, as well as an orientation encoded by the source function s(e) and the target function t(e) this implements the system of ODEs:\n\nfracdx_ndt = dx_n\n\nfracdl_edt = dl_e\n\nwith dx and dl calculated by\n\nl_s = l_e text if  s(e) = n\n\nl_t = l_e text if  t(e) = n\n\nnodesn(dx_n x_n l_s l_t p_n t)\n\nlinese(dl_e l_e x_s(e) x_t(e) p_e t)\n\nAlternative design:\n\nSomething that relaxes to a diffusive network would for example be implemented by\n\nlines = (dl, l, x_1, x_2) -> dl .= 1000. * ((x_1 - x_2) - l)\nagg = (list_1, list_2) -> sum(list_1) - sum(list_2)\n\n\n\n\n\n"
},

{
    "location": "#NetworkDynamics.network_dynamics_on_the_line",
    "page": "NetworkDynamics.jl",
    "title": "NetworkDynamics.network_dynamics_on_the_line",
    "category": "type",
    "text": "edge_list: list of edges nodes: scalar function x rightarrow nodes(x) lines: scalar function x rightarrow lines(x)\n\nImplements the ODE fracdx_idt = nodes(x_i) + sum_j A_ij lines(x_i - x_j)\n\nwhere A is the adjacency matrix of B\n\n\n\n\n\n"
},

{
    "location": "#NetworkDynamics.network_dynamics_on_lines-Tuple{AbstractGraph,Any}",
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
