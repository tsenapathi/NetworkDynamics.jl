var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "General",
    "title": "General",
    "category": "page",
    "text": ""
},

{
    "location": "#NetworkDynamics.jl-1",
    "page": "General",
    "title": "NetworkDynamics.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "General",
    "title": "Overview",
    "category": "section",
    "text": "This package implements functions for defining and studying dynamics on networks. The key construction is a callable struct compatible with the DifferentialEquations.jl calling syntax.nd = network_dynamics(vertices!, edges!, g)\nnd(dx, x, p, t)The first two parameters are the functions, or function arrays from which a network dynamics is built. The last parameter g is a graph encoding the network constructed with the LightGraphs.jl package. Note that the type network_dynamics is only a placeholder for the dynamics types we will now specify."
},

{
    "location": "#Static-lines-1",
    "page": "General",
    "title": "Static lines",
    "category": "section",
    "text": ""
},

{
    "location": "#Scalar-variables-1",
    "page": "General",
    "title": "Scalar variables",
    "category": "section",
    "text": "We will first look into the case where there is only a scalar variable on each vertex. A dynamical network with static lines (static meaning that the current on an edge depends solely on the values on the nodes it connects) is created via the scalarstaticlines function:ssl = scalar_static_lines(vertices!, edges!, g)The functions vertices! and edges! are of the form:vertices![n](dv[n],v[n],e_s[n],e_t[n],p,t)\nedges![m](e[m],v_s,v_t,p,t)  Specifically, the given variables are:e_s[n] = [e[m] if s[m] == n for m in 1:length(lines!)]\ne_t[n] = [e[m] if t[m] == n for m in 1:length(lines!)]\nv_s= v[s[m]]\nv_t= v[t[m]]The vectors s and t contain the information about the source and target of each edge, i.e. s[1] == 2 -> The source of edge 1 is vertex 2. The function creates these vectors from the given graph, they can be accessed via the calling syntax ssl.se or ssl.te. The vectors es[n] and et[n] are containing the in- and outgoing edge values (or currents) of vertex n in the form of an array. Thus, one would classically sum over these in vertices!, but one is not restricted on doing this.For example, a system of equations describing a simple diffusive network would be:using LightGraphs\ng= barabasi_albert(10,5)\nvertices! = [(dv,v,l_s,l_t,p,t) -> dv .= sum(e_s) .- sum(e_t) for vertex in vertices(g)]\nedges! = [(e,v_s,v_t,p,t) -> e .= v_s .- v_t for edge in edges(g)]Here, the diffusiveness lies within the lines! function. It states that there is only a current between two nodes if these nodes have a non-equal value. This current then ultimatively leads to an equilibrium in which the value on any connected node is equal.Note that one should (for performance reasons) and actually needs to put a dot before the mathematical operators. This is due to the use of views in the internals of the scalarstaticlines function.We finally want to solve the defined diffusive system. This we do by using the well-known package DifferentialEquations.jl (see here). We also need to specify a set of initial values x0 as well as a time interval t for which we are solving the problem:using DifferentialEquations\nusing Plots\nx0 = rand(10)\nt = (0.,2.)\nssl = scalar_static_lines(vertices!,edges!,g)\nssl_prob = ODEProblem(ssl,x0,t)\nsol = solve(ssl_prob)\nplot(sol, legend = false)(Image: )As one would expect in a diffusive network, the values on the vertices converge."
},

{
    "location": "#Vector-variables-1",
    "page": "General",
    "title": "Vector variables",
    "category": "section",
    "text": "In most cases, one is interested in problems with multiple variables on the vertices, one can deal with these problems by using the, compared to the previous more general, static_lines function:sl = static_lines(vertices!, edges!, g, dim_v, dim_e)In comparison to the scalarstaticlines function, here one also has to specify the dimension of the variables on each vertex and edge. The vectors dimv and dime have entries for every vertex and edge in the graph, the value of the entry fixes the number of variables. As a simple example, let\'s again look at a problem with scalar variables. The functions vertices! and edges! are as before:dim_v = ones(Int32, length(vertices!))\ndim_e = ones(Int32, length(edges!))\nsl = static_lines(vertices!, edges!, g, dim_v, dim_e)Every entry of the dimension vectors is one, so the variables on each vertex and edge are 1-dimensional. Note that one has to specify the ones type as Integer.The purpose of this function though is that we want treat higher dimensional problems. So let\'s look at an 2-dimensional problem.dim_v = 2 * ones(Int32, length(vertices!))\ndim_e = 2 * ones(Int32, length(edges!))Now we fixed the dimension on every vertex and edge to 2, following from that the vertices! and edges! functions are also 2-dimensional. We again want to look at a diffusion problem, for that we have to change the vertices! function  due to a problem occuring with the sum function, namely that it is not able to deal with an Any-type empty set:function vertex!(dv, v, e_s, e_d, p, t)\n    dv .= 0\n    for e in e_s\n        dv .-= e\n    end\n    for e in e_d\n        dv .+= e\n    end\n    nothing\nend\n\nvertices! = [vertex! for vertex in vertices(g)]The edge! function stays the same. The arguments appearing in the functions are now all 2-dimensional, so in principle one could put for example non-diagonal matrices into play to establish an interaction between the different variables. As we did not do that here, what we expect is that we again get the same solution as before, but twice. To solve it, we now need to give the solver twice as many initial values x0, the pattern in the x vector goes [vertex1variable1,vertex1variable2,vertex2_variable1,...]:x0 = rand(20)\nt = (0.,2.)\nsl = static_lines(vertices!,edges!,g,dim_v,dim_e)\nsl_prob = ODEProblem(ssl,x0,t)\nsol = solve(ssl_prob)\nplot(sol, legend = false)#figureAs we see, we get the solution of two independent diffusive networks. (Here one could put a more complex exmaple to showcase the functionality.)"
},

{
    "location": "#Dynamic-lines-1",
    "page": "General",
    "title": "Dynamic lines",
    "category": "section",
    "text": ""
},

{
    "location": "#Scalar-variables-2",
    "page": "General",
    "title": "Scalar variables",
    "category": "section",
    "text": "In general, currents do not solely depend on the vertex values they are connecting, but rather depend on its own value in some sort. For the case of scalar variables, we may use the function scalardynamiclines:sdl = scalar_dynamic_lines(vertices!,edges!,g)The function arguments are now of the following form:vertices![n](dv[n],v[n],e_s[n],e_t[n],p,t)\nedges![m](de[m],e[m],v_s,v_t,p,t)Compared to the static lines case with scalar variables, the vertices! function keeps its structure whereas the edges! function gets the new argument de[m]. This de[m] is the derivative of the edge value of edge m. Let\'s look at a simple example: A system with dynamic lines which decay to the usual diffusive system:vertices! = [(dv,v,l_s,l_t,p,t) -> dv .= sum(e_s) .- sum(e_t) for vertex in vertices(g)]\nedges! = [(de,e,v_s,v_t,p,t) -> de .= 1000*(v_s .- v_t .- e) for edge in edges(g)]The change compared to the example for the static case should be clear; the factor of 1000 is just accelerating the decay. Again, we can quite simply solve this system. One has to be aware though that now one needs initial values for the vertices and the edges! These are given in the order x0 = [vertex1,vertex2,...,edge1,edge2,...]:g = barabasi_albert(10,5) #generates a graph with 10 vertices and 25 edges\nx0 = rand(10 + 25)\nt = (0.,2.)\nsdl = scalar_dynamic_lines(vertices!,edges!,g)\nsdl_prob = ODEProblem(sdl,x0,t)\nsol = solve(sdl_prob)\nplot(sol, legend = false , vars = 1:10)(Hier sollte ein Bild sein)We see that the plot looks pretty much the same as for the static lines case. That is, because we included the factor of 1000 in the edges! function. Note that we added the argument vars to the plot function, this gives us solely the first 10 arguments of x which are the vertices. One could also get just the edge values by writing vars = 11:35 if one wishes."
},

{
    "location": "#Vector-variables-2",
    "page": "General",
    "title": "Vector variables",
    "category": "section",
    "text": "The step here is not a hard one, if one read through the previous Vector variables section. We can treat a system of vector variables with dynamic lines with the function dynamic_lines:dl = dynamic_lines(vertices!,edges!,g,dim_v,dim_e)One has to apply the same change to the vertices! function as for the static_lines function. Otherwise, everything should be clear. For the example, we take the decaying dynamic lines and just make two independent networks as for the Static lines:dim_v = 2 * ones(Int32, length(vertices!))\ndim_e = 2 * ones(Int32, length(edges!))\ng = barabasi_albert(10,5)\n\nfunction vertex!(dv, v, e_s, e_d, p, t)\n    dv .= 0\n    for e in e_s\n        dv .-= e\n    end\n    for e in e_d\n        dv .+= e\n    end\n    nothing\nend\n\nvertices! = [vertex! for vertex in vertices(g)]\nedges! = [(de,e,v_s,v_t,p,t) -> de .= 1000*(v_s .- v_t .- e) for edge in edges(g)]\n\ndl = dynamic_lines(vertices!,edges!,g,dim_v,dim_e)\n\nx0 = rand(10 + 10 + 25 + 25)\nt= (0.,2.)\ndl_prob = ODEProblem(dl,x0,t)\nsol= solve(dl_prob)\nplot(sol, legend = false, vars = 1:20)(Bild)We get the same pattern as for the scalar case, just twice."
},

{
    "location": "Functions_and_Constructors/#",
    "page": "-",
    "title": "-",
    "category": "page",
    "text": ""
},

{
    "location": "Functions_and_Constructors/#Functions-1",
    "page": "-",
    "title": "Functions",
    "category": "section",
    "text": "The Dynamics for the whole Network is constructed from functions for the single vertices and edges. There are several types:ODEVertex(vertexfunction!, dimension, massmatrix, sym)\nStaticEdge(edgefunction!, dimension)\nODEEdge(edgefunction!, dimension, massmatrix, sym)\n````\n\n\n### ODEVertex\n\nThe arguments mean the following: vertexfunction is catching the dynamics of a single vertex depending on the vertex value itself as well as in- and outgoing currents (or edges). An example for such a function would be:\njulia function vertexfunction!(dv, v, es, ed, p, t)   dv .= 0   for e in es     dv .-= e   end   for e in ed     dv .+= e   end end\nThe e_s and e_d are arrays containing the edges that have the decribed vertex as source and destination. Other arguments coincide with the usual ODE function arguments. The vertexfunction given to ODEVertex always needs to have the shown argument structure. Note the importance of the broadcast structure of the equations (the dot before every operator), this is necessary due to the use of views in the internal functions.\n\ndimension is the number of Variables on the Vertex.\n\nmassmatrix is the mass matrix M, i.e.\n@docs M*dv = vertexfunction!\nsym are the symbols of the Vertex. If one had for example a vertex with a frequency and some angle, one would construct sym via:\n@docs sym = [:omega, :phi]\nThis makes it easier to later fish out the interesting variables one wants to look at.\n\nOne may also call ODEVertex via:\n@docs ODEVertex(vertexfunction!, dimension)\nThe function then defaults to using the identity as mass matrix and [:v for i in 1:dimension] as symbols.\n\n\n### StaticEdge\n\nStatic here means, that the edge value described by edgefunction! solely depends on the vertex values the edge connects. One very simple and natural example is a diffusive system:\n@julia edgefunction! = (e, vs, vd, p, t) -> e .= vs .- vd\nv_s and v_d are the vertex values of the edges source and destination. There is no derivative of the edge value involved, hence we call these problems static.\n\ndimension: see ODEVertex\n\n### ODEEdge\n\nFor Problems where the derivative of the edge value indeed depends on the edge value itself, we use the ODEEdge function. Another simple and natural example is a system that quickly diffuses to the static case:\n@julia edgefunction! = (de, e, vs, vd, p, t) -> de .= 1000 * (vs .- vd .- e)\ndimension: see ODEVertex\n\nmassmatrix: see ODEVertex\n\nsym: see ODEVertex\n\nAlso, one can construct an ODEEdge by only giving the first two arguments:\n@docs ODEEdge(edgefunction!, dimension)\nThen the function defaults to using the identity as mass matrix as well as using [:e for in 1:dimension] as sym.\n\n\n\n\n## Constructor\n\nThe central constructor of this package is network_dynamics, this function demands an array of VertexFunction and EdgeFunction as well as a graph (see LightGraphs), and returns an ODEFunction which one can easily solve via the tools given in DifferentialEquations.jl. One calls it via:\n@docs network_dynamics(Array{VertexFunction}, Array{EdgeFunction}, graph) ```"
},

{
    "location": "StaticEdges/#",
    "page": "-",
    "title": "-",
    "category": "page",
    "text": ""
},

{
    "location": "StaticEdges/#Static-lines-1",
    "page": "-",
    "title": "Static lines",
    "category": "section",
    "text": ""
},

{
    "location": "StaticEdges/#Scalar-variables-1",
    "page": "-",
    "title": "Scalar variables",
    "category": "section",
    "text": "We will first look into the case where there is only a scalar variable on each vertex. A dynamical network with static lines (static meaning that the current on an edge depends solely on the values on the nodes it connects) is created via the scalarstaticlines function:ssl = scalar_static_lines(vertices!, edges!, g)The functions vertices! and edges! are of the form:vertices![n](dv[n],v[n],e_s[n],e_t[n],p,t)\nedges![m](e[m],v_s,v_t,p,t)  Specifically, the given variables are:e_s[n] = [e[m] if s[m] == n for m in 1:length(lines!)]\ne_t[n] = [e[m] if t[m] == n for m in 1:length(lines!)]\nv_s= v[s[m]]\nv_t= v[t[m]]The vectors s and t contain the information about the source and target of each edge, i.e. s[1] == 2 -> The source of edge 1 is vertex 2. The function creates these vectors from the given graph, they can be accessed via the calling syntax ssl.se or ssl.te. The vectors es[n] and et[n] are containing the in- and outgoing edge values (or currents) of vertex n in the form of an array. Thus, one would classically sum over these in vertices!, but one is not restricted on doing this.For example, a system of equations describing a simple diffusive network would be:using LightGraphs\ng= barabasi_albert(10,5)\nvertices! = [(dv,v,l_s,l_t,p,t) -> dv .= sum(e_s) .- sum(e_t) for vertex in vertices(g)]\nedges! = [(e,v_s,v_t,p,t) -> e .= v_s .- v_t for edge in edges(g)]Here, the diffusiveness lies within the lines! function. It states that there is only a current between two nodes if these nodes have a non-equal value. This current then ultimatively leads to an equilibrium in which the value on any connected node is equal.Note that one should (for performance reasons) and actually needs to put a dot before the mathematical operators. This is due to the use of views in the internals of the scalarstaticlines function.We finally want to solve the defined diffusive system. This we do by using the well-known package DifferentialEquations.jl (see here). We also need to specify a set of initial values x0 as well as a time interval t for which we are solving the problem:using DifferentialEquations\nusing Plots\nx0 = rand(10)\nt = (0.,2.)\nssl = scalar_static_lines(vertices!,edges!,g)\nssl_prob = ODEProblem(ssl,x0,t)\nsol = solve(ssl_prob)\nplot(sol, legend = false)(Image: )As one would expect in a diffusive network, the values on the vertices converge."
},

{
    "location": "DynamicEdges/#",
    "page": "-",
    "title": "-",
    "category": "page",
    "text": ""
},

{
    "location": "DynamicEdges/#Dynamic-lines-1",
    "page": "-",
    "title": "Dynamic lines",
    "category": "section",
    "text": ""
},

{
    "location": "DynamicEdges/#Scalar-variables-1",
    "page": "-",
    "title": "Scalar variables",
    "category": "section",
    "text": "In general, currents do not solely depend on the vertex values they are connecting, but rather depend on its own value in some sort. For the case of scalar variables, we may use the function scalardynamiclines:sdl = scalar_dynamic_lines(vertices!,edges!,g)The function arguments are now of the following form:vertices![n](dv[n],v[n],e_s[n],e_t[n],p,t)\nedges![m](de[m],e[m],v_s,v_t,p,t)Compared to the static lines case with scalar variables, the vertices! function keeps its structure whereas the edges! function gets the new argument de[m]. This de[m] is the derivative of the edge value of edge m. Let\'s look at a simple example: A system with dynamic lines which decay to the usual diffusive system:vertices! = [(dv,v,l_s,l_t,p,t) -> dv .= sum(e_s) .- sum(e_t) for vertex in vertices(g)]\nedges! = [(de,e,v_s,v_t,p,t) -> de .= 1000*(v_s .- v_t .- e) for edge in edges(g)]The change compared to the example for the static case should be clear; the factor of 1000 is just accelerating the decay. Again, we can quite simply solve this system. One has to be aware though that now one needs initial values for the vertices and the edges! These are given in the order x0 = [vertex1,vertex2,...,edge1,edge2,...]:g = barabasi_albert(10,5) #generates a graph with 10 vertices and 25 edges\nx0 = rand(10 + 25)\nt = (0.,2.)\nsdl = scalar_dynamic_lines(vertices!,edges!,g)\nsdl_prob = ODEProblem(sdl,x0,t)\nsol = solve(sdl_prob)\nplot(sol, legend = false , vars = 1:10)(Hier sollte ein Bild sein)We see that the plot looks pretty much the same as for the static lines case. That is, because we included the factor of 1000 in the edges! function. Note that we added the argument vars to the plot function, this gives us solely the first 10 arguments of x which are the vertices. One could also get just the edge values by writing vars = 11:35 if one wishes."
},

{
    "location": "DynamicEdges/#Vector-variables-1",
    "page": "-",
    "title": "Vector variables",
    "category": "section",
    "text": "The step here is not a hard one, if one read through the previous Vector variables section. We can treat a system of vector variables with dynamic lines with the function dynamic_lines:dl = dynamic_lines(vertices!,edges!,g,dim_v,dim_e)One has to apply the same change to the vertices! function as for the static_lines function. Otherwise, everything should be clear. For the example, we take the decaying dynamic lines and just make two independent networks as for the Static lines:dim_v = 2 * ones(Int32, length(vertices!))\ndim_e = 2 * ones(Int32, length(edges!))\ng = barabasi_albert(10,5)\n\nfunction vertex!(dv, v, e_s, e_d, p, t)\n    dv .= 0\n    for e in e_s\n        dv .-= e\n    end\n    for e in e_d\n        dv .+= e\n    end\n    nothing\nend\n\nvertices! = [vertex! for vertex in vertices(g)]\nedges! = [(de,e,v_s,v_t,p,t) -> de .= 1000*(v_s .- v_t .- e) for edge in edges(g)]\n\ndl = dynamic_lines(vertices!,edges!,g,dim_v,dim_e)\n\nx0 = rand(10 + 10 + 25 + 25)\nt= (0.,2.)\ndl_prob = ODEProblem(dl,x0,t)\nsol= solve(dl_prob)\nplot(sol, legend = false, vars = 1:20)(Bild)We get the same pattern as for the scalar case, just twice."
},

]}
