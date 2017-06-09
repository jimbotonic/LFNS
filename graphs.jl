using Graphs, Logging, Distances

@Logging.configure(level=DEBUG)

# vertex type
type Bus
	id::Int64
	name::AbstractString
	# 0:PQ, 1:Q \theta, 2: PV, 3: V\theta
	bus_type::Int
	init_voltage::Float64
	final_voltage::Float64
	# base voltage in KV
	base_voltage::Float64
	angle::Float64
	# ENTSOE convention: Re(load) > 0, Re(generation) < 0
	load::Complex{Float64}
	generation::Complex{Float64}
	Q_min::Float64
	Q_max::Float64
	P_min::Float64
	P_max::Float64
	sh_conductance::Float64
	sh_susceptance::Float64
	# normalized longitude and latitude (between 0 and 1)
	lng::Float64
	lat::Float64
	
	# default constructor
	function Bus(id::Int64, name::AbstractString, bus_type::Int, init_voltage::Float64, final_voltage::Float64, base_voltage::Float64, angle::Float64, load::Complex{Float64}, generation::Complex{Float64}, Q_min::Float64, Q_max::Float64, P_min::Float64, P_max::Float64, sh_conductance::Float64, sh_susceptance::Float64, lng::Float64, lat::Float64)
		new(id, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, P_min, P_max, sh_conductance, sh_susceptance, lng, lat)
	end

	# simple constructors
	function Bus(id::Int64, t::Float64, p::Float64)
		if p < 0
			# load
			new(id, "$id", 2, 1., 1., 1., t, -p, 0., 0., 0., 0., 0., 0., 0., 0., 0.)
		else
			# generation
			new(id, "$id", 2, 1., 1., 1., t, 0., p, 0., 0., 0., 0., 0., 0., 0., 0.)
		end
	end 
	
	function Bus(id::Int64, name::AbstractString)
		new(id, name, 2, 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
	end 
end

# edge type
type Line
	id::Int64
	# source and target bus
	source::Bus
	target::Bus
	# 0: normal, 1: transformer 
	line_type::Int
	# default on 
	# 0,1,2-> ON
	# 7,8,9 -> OFF  
	line_status::Bool
	# Z = R + iX, Y = 1/Z
	admittance::Complex{Float64}
	sh_susceptance::Float64
	# 1 for a standard power line
	# transformer on the source side (IEEE)
	s_ratio::Complex{Float64}
	# transformer on the target side (ENTSOE)
	t_ratio::Complex{Float64}

	# default constructor
	function Line(id::Int64, source::Bus, target::Bus, line_type::Int, line_status::Bool, admittance::Complex{Float64}, sh_susceptance::Float64, s_ratio::Complex{Float64}, t_ratio::Complex{Float64})
		new(id, source, target, line_type, line_status, admittance, sh_susceptance, s_ratio, t_ratio)
	end

	# simple constructor
	function Line(id::Int64, source::Bus, target::Bus, y::Complex{Float64})
		new(id, source, target, 0, true, y, 0., 1., 1.)
	end
end

# redefine some of the base functions
Graphs.source(e::Line, g::Graphs.AbstractGraph{Bus,Line}) = e.source
Graphs.target(e::Line, g::Graphs.AbstractGraph{Bus,Line}) = e.target
Graphs.revedge(e::Line) = Line(e.id, e.target, e.source, e.line_type, e.line_status, e.admittance, e.sh_susceptance, e.s_ratio, e.t_ratio)
Graphs.vertex_index(v::Bus) = v.id
Graphs.edge_index(e::Line) = e.id

# get the connected component vertex ids containing the slack bus
function get_slack_component_ids(g::Graphs.AbstractGraph{Bus,Line})
	cs = connected_components(g)
	for c in cs
		if in(3,Int[v.bus_type for v in c])
			return Int[v.id for v in c]
		end
	end
end


# get the vector of active powers
function get_active_powers(g::Graphs.AbstractGraph{Bus,Line})
	n = length(vertices(g))
	P = zeros(Float64,n)
	for v in vertices(g)
		P[v.id] = real(v.generation)-real(v.load)
	end
	
	return P
end

# get the vector of reactive powers
function get_reactive_powers(g::Graphs.AbstractGraph{Bus,Line})
	n = length(vertices(g))
	Q = zeros(Float64,n)
	for v in vertices(g)
		Q[v.id] = imag(v.generation)-imag(v.load)
	end
	
	return Q
end

# get the admittance matrix
function get_admittance_matrix(g::Graphs.AbstractGraph{Bus,Line})
	n = length(vertices(g))
	Y = spzeros(Complex{Float64},n,n)

##	YY = SparseMatrixCSC{Complex{Float64},Int64}(Y)
##	Using this type convertion from an Array of Complex Floats, we obtained some non-zero values at non-deterministic locations in matrix YY.

	for e in edges(g)
		s = e.source.id
		t = e.target.id
		Y[s,t] = -e.admittance
		Y[t,s] = -e.admittance
	end
	
	return Y - spdiagm(vec(sum(Y,1)))
end	


# get the vector of angles
function get_angles(g::Graphs.AbstractGraph{Bus,Line})
	T = Array{Float64,1}()
	for v in vertices(g)
		push!(T,v.angle)
	end
	return mod(T+pi,2*pi)-pi
end
	
# get the vector of voltage amplitudes
function get_voltages(g::Graphs.AbstractGraph{Bus,Line})
	V = Array{Float64,1}()
	for v in vertices(g)
		push!(V,v.base_voltage)
	end
	return V
end
	
# get the adjacency matrix
function get_adjacency_matrix(g::Graphs.AbstractGraph{Bus,Line})
	n = length(vertices(g))
	A = zeros(Int64,n,n)
	for edge in edges(g)
		s = edge.source.id
		t = edge.target.id
		A[s,t] = 1
		A[t,s] = 1
	end
	return A
end

# get the incidence matrix
function get_incidence_matrix(g::Graphs.AbstractGraph{Bus,Line})
	n = length(vertices(g))
	m = length(edges(g))
	I = zeros(Int64,n,m)
	for edge in edges(g)
		s = edge.source.id
		t = edge.target.id
		I[s,edge.id] = 1
		I[t,edge.id] = -1
	end
	return I
end

# get the principal component as an array of vertices
function get_principal_component(g::Graphs.AbstractGraph{Bus,Line})
	cs = connected_components(g)
	return cs[indmax([length(c) for c in cs])]
end

# get subgraph
# 
# if keep_vertices is true, get the sugraph induced by the specified list of vertices
# otherwise, get the subgraph induced by the set of vertices remaining after the specified list is removed
function get_subgraph(g::Graphs.AbstractGraph{Bus,Line}, vids::Array{Int64,1}, keep_vids::Bool=true)
	if keep_vids
		# sorted array of subgraph vertex ids
		svids = sort(vids)
	else
		gvids = Int64[v.id for v in vertices(g)]
		svids = sort(setdiff(gvids,vids)) 
	end

	nvs = Bus[]
	nes = Line[]

	# old to new vertex id mappings
	o2n_ids = Dict{Int64, Int64}()

	vcounter = 1
	for v in vertices(g)
		if v.id in svids
			o2n_ids[v.id] = vcounter
			nv = Bus(vcounter, v.name, v.bus_type, v.init_voltage, v.final_voltage, v.base_voltage, v.angle, v.load, v.generation, v.Q_min, v.Q_max, v.P_min, v.P_max, v.sh_conductance, v.sh_susceptance, v.lng, v.lat)
			push!(nvs, nv)
			vcounter += 1
		end
	end

	ecounter = 1
	for edge in edges(g)
		sid = edge.source.id
		tid = edge.target.id
		# if both endpoints of the edge are in the subgraph
		if length(searchsorted(svids, sid)) > 0 && length(searchsorted(svids, tid)) > 0
			source = nvs[o2n_ids[sid]]
			target = nvs[o2n_ids[tid]]
			ne = Line(ecounter, source, target, edge.line_type, edge.line_status, edge.admittance, edge.sh_susceptance, edge.s_ratio, edge.t_ratio)
			push!(nes, ne)
			ecounter += 1
		end
	end
	
	return graph(nvs, nes, is_directed=false)
end

# prune the graph by removing the specified list of  edge ids
#
# NB: there is no need to renumber the vertex ids
function get_pruned_graph(g::Graphs.AbstractGraph{Bus,Line}, reids::Array{Int64,1})
	# sort array of subgraph edge ids to be removed
	sort!(reids)
	nvs = Bus[]
	nes = Line[]

	for v in vertices(g)
		nv = Bus(v.id, v.name, v.bus_type, v.init_voltage, v.final_voltage, v.base_voltage, v.angle, v.load, v.generation, v.Q_min, v.Q_max, v.P_min, v.P_max, v.sh_conductance, v.sh_susceptance, v.lng, v.lat)
		push!(nvs, nv)
	end

	ecounter = 1
	for edge in edges(g)
		# if the current edge is not in the sorted list
		if length(searchsorted(reids, edge.id)) == 0
			ne = Line(ecounter, edge.source, edge.target, edge.line_type, edge.line_status, edge.admittance, edge.sh_susceptance, edge.s_ratio, edge.t_ratio)
			push!(nes, ne)
			ecounter += 1
		end
	end
	
	return graph(nvs, nes, is_directed=false)
end

# prune the graph by removing the specified list of edges specified as pairs (source_id,target_id
#
# NB: there is no need to renumber the vertex ids
function get_pruned_graph(g::Graphs.AbstractGraph{Bus,Line}, redges::Array{Pair{Int,Int},1})
	nvs = Bus[]
	nes = Line[]

	for v in vertices(g)
		# id, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, P_min, P_max, sh_conductance, sh_susceptance, lat, lng)
		#nv = copy(v)
		nv = Bus(v.id, v.name, v.bus_type, v.init_voltage, v.final_voltage, v.base_voltage, v.angle, v.load, v.generation, v.Q_min, v.Q_max, v.P_min, v.P_max, v.sh_conductance, v.sh_susceptance, v.lng, v.lat)
		push!(nvs, nv)
	end

	ecounter = 1
	for edge in edges(g)
		# if the current edge is not in the sorted list
		if !(Pair{Int,Int}(edge.source.id,edge.target.id) in redges)
			#ne = copy(edge)
			#ne.id = ecounter
			ne = Line(ecounter, edge.source, edge.target, edge.line_type, edge.line_status, edge.admittance, edge.sh_susceptance, edge.s_ratio, edge.t_ratio)
			push!(nes, ne)
			ecounter += 1
		end
	end
	
	return graph(nvs, nes, is_directed=false)
end

# get avg / max degree of the specified graph
function get_avg_min_max_degree(g::Graphs.AbstractGraph{Bus,Line})
	degrees = Int64[out_degree(v,g) for v in vertices(g)]
	return mean(degrees), minimum(degrees), maximum(degrees)
end

# get the cycle base of the specified graph
# 
# see http://graphsjl-docs.readthedocs.org/en/latest/algorithms.html
#
# Paton algorithm: 
# 	http://www.cs.kent.edu/~dragan/GraphAn/CycleBasis/p514-paton.pdf
# 	https://code.google.com/p/niographs/source/browse/src/main/java/net/ognyanov/niographs/PatonCycleBase.java
function get_cycle_base(g::Graphs.AbstractGraph{Bus,Line})
	# get minimum spanning tree edges
	ew = ones(length(edges(g)))
	ke, kw = prim_minimum_spantree(g, ew, vertices(g)[1])
	
	# get edges not belonging to the tree (edges(g) \ edges(sp))
	reids = setdiff(Int64[edge.id for edge in edges(g)],Int64[edge.id for edge in ke])

	# extract spanning tree
	gt = get_pruned_graph(g, reids)

	#@debug("# vertices (g): ", length(vertices(g)))
	#@debug("# edges (g): ", length(edges(g)))
	#@debug("---")
	#@debug("# vertices (mst): ", length(vertices(gt)))
	#@debug("# edges (mst): ", length(edges(gt)))
	#@debug("---")
	#@debug("# shortcut edges: ", length(reids))
	
	# set of fundamental cycles
	cycles = Array{Array{Int64,1},1}()

	for eid in reids
		edge = edges(g)[eid]
		# select the least edge vertex
		s,t = sort([edge.source.id,edge.target.id])
		# pergorm Djikstra shortest path algorithm from the source s
		dsp = dijkstra_shortest_paths(gt, vertices(gt)[s])
		c = enumerate_indices(dsp.parent_indices, t)
		push!(cycles, c)
	end

	return cycles
end

# naïve implementation of Pagerank algorithm for an undirected graph
function my_pagerank(g::Graphs.AbstractGraph{Bus,Line},pr::Array{Float64,1},damping::Float64=0.85,epsilon::Float64=1e-4)
	vs = vertices(g)
	n = length(vs)
	pr2 = copy(pr)
	while true
		for v in vs
			nv = 0.
			# get v children in the reverse graph
			nei = out_neighbors(v,g)
			if length(nei) > 0
				for p in nei
					nv +=  pr[p.id]/out_degree(p,g)
				end
			end
			pr2[v.id] = (1-damping)/n+damping*nv
		end
		d = euclidean(pr,pr2) 
		d <= epsilon && break
		pr = copy(pr2)
	end
	return pr2
end


# creates a copy of a graph
function copy_graph(g::Graphs.AbstractGraph{Bus,Line})
	vs = Bus[]
	es = Line[]
	
	for v in vertices(g)
		push!(vs,Bus(copy(v.id),copy(v.name),copy(v.bus_type),copy(v.init_voltage),copy(v.final_voltage),copy(v.base_voltage),copy(v.angle),copy(v.load),copy(v.generation),copy(v.Q_min),copy(v.Q_max),copy(v.P_min),copy(v.P_max),copy(v.sh_conductance),copy(v.sh_susceptance),copy(v.lng),copy(v.lat)))
	end
	
	for e in edges(g)
		push!(es,Line(copy(e.id),vs[e.source.id],vs[e.target.id],copy(e.line_type),copy(e.line_status),copy(e.admittance),copy(e.sh_susceptance),copy(e.s_ratio),copy(e.t_ratio)))
	end
	
	return graph(vs,es,is_directed=is_directed(g))
end

	
