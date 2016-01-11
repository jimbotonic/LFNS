using Graphs, Logging

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
	
	# default constructor
	function Bus(id::Int64, name::AbstractString, bus_type::Int, init_voltage::Float64, final_voltage::Float64, base_voltage::Float64, angle::Float64, load::Complex{Float64}, generation::Complex{Float64}, Q_min::Float64, Q_max::Float64, P_min::Float64, P_max::Float64, sh_conductance::Float64, sh_susceptance::Float64)
		new(id, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, P_min, P_max, sh_conductance, sh_susceptance)
	end

	# simple constructors
	function Bus(id::Int64, t::Float64, p::Float64)
		if p < 0
			# load
			new(id, "$id", 2, 1., 1., 1., t, -p, 0., 0., 0., 0., 0., 0., 0.)
		else
			# generation
			new(id, "$id", 2, 1., 1., 1., t, 0., p, 0., 0., 0., 0., 0., 0.)
		end
	end 
	
	function Bus(id::Int64, name::AbstractString)
		new(id, name, 2, 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
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

# get the principal component vertex ids
function get_principal_component(g::Graphs.AbstractGraph{Bus,Line})
	cs = connected_components(g)
	return cs[indmax([length(c) for c in cs])]
end

# get the subgraph induced by the specified vertices
function get_subgraph(g::Graphs.AbstractGraph{Bus,Line}, vs::Array{Bus,1})
	# sorted array of subgraph vertex ids
	svids = sort(Int64[v.id for v in vs])
	nvs = Bus[]
	nes = Line[]

	# old to new vertex id mappings
	o2n_ids = Dict{Int64, Int64}()

	vcounter = 1
	for v in vs
		o2n_ids[v.id] = vcounter
		nv = Bus(vcounter, v.name, v.bus_type, v.init_voltage, v.final_voltage, v.base_voltage, v.angle, v.load, v.generation, v.Q_min, v.Q_max, v.P_min, v.P_max, v.sh_conductance, v.sh_susceptance)
		push!(nvs, nv)
		vcounter += 1
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
function get_pruned_graph(g::Graphs.AbstractGraph{Bus,Line}, reids::Array{Int64,1})
	# sort array of subgraph edge ids to be removed
	sort!(reids)
	nvs = Bus[]
	nes = Line[]

	for v in vertices(g)
		nv = Bus(v.id, v.name, v.bus_type, v.init_voltage, v.final_voltage, v.base_voltage, v.angle, v.load, v.generation, v.Q_min, v.Q_max, v.P_min, v.P_max, v.sh_conductance, v.sh_susceptance)
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

	@debug("# vertices (g): ", length(vertices(g)))
	@debug("# edges (g): ", length(edges(g)))
	@debug("---")
	@debug("# vertices (gt): ", length(vertices(gt)))
	@debug("# edges (gt): ", length(edges(gt)))
	@debug("---")
	@debug("# shortcut edges: ", length(reids))
	
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

function direct_cycles(cycles::Array{Array{Int64,1},1})
	# find least cycle (cycle containing the least vertex whose sum is minimum)
	cycles1 = filter(x -> 1 in x, cycles)
	lc = filter(x -> sum(x) == minimum(sum(cycles)))

	
end

# generate a double cycle with a bus where p is injected
function generate_double_cycle(l::Int,c::Int,r::Int,p::Float64)
	vs = Bus[]
	es = Line[]	
	ecounter = 1

	# left branch 1:l
	for i in 1:l
		# t=0., p=0.
		push!(vs,Bus(i,0.,0.))
	end
	for i in 2:l
		push!(es, Line(ecounter,vs[i-1],vs[i],1im))
		ecounter += 1
	end
	# central branch (l+1):(c+l)
	for i in (l+1):(c+l)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+2):(c+l)
		push!(es, Line(ecounter,vs[i-1],vs[i],1im))
		ecounter += 1
	end
	# right branch (l+c+1):(c+l+r)
	for i in (l+c+1):(c+l+r)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+c+2):(c+l+r)
		push!(es, Line(ecounter,vs[i-1],vs[i],1im))
		ecounter += 1
	end

	# 2 last remaining vertices
	b_in = Bus(l+c+r+1,0.,p)
	b_out = Bus(l+c+r+2,0.,-p)

	push!(es, Line(ecounter,vs[l+c+r+1],vs[1],1im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+1],vs[l+1],1im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+1],vs[l+c+1],1im))
	ecounter += 1
	
	push!(es, Line(ecounter,vs[l+c+r+2],vs[l],1im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+2],vs[l+c],1im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+2],vs[l+c+r],1im))
	
	return graph(vs, es, is_directed=false)
end

