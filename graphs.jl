using Graphs

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
	load::Complex{Float64}
	generation::Complex{Float64}
	Q_min::Float64
	Q_max::Float64
	P_min::Float64
	P_max::Float64
	sh_conductance::Float64
	sh_susceptance::Float64
end

# edge type
type Line
	source::Bus
	target::Bus
	# 0: normal, 1: transformer 
	line_type::Int
	# default on 
	# 0,1,2-> ON
	# 7,8,9 -> OFF  
	line_status::Bool
	resistance::Float64
	reactance::Float64
	sh_susceptance::Float64
	# 1 for a standard power line
	# transformer on the source side (IEEE)
	s_ratio::Complex{Float64}
	# transformer on the target side (ENTSOE)
	t_ratio::Complex{Float64}
end

# redefine some of the base functions
Graphs.source(e::Line, g::Graphs.AbstractGraph{Bus,Line}) = e.source
Graphs.target(e::Line, g::Graphs.AbstractGraph{Bus,Line}) = e.target
Graphs.revedge(e::Line) = Line(e.target, e.source, e.line_type, e.line_status, e.resistance, e.reactance, e.sh_susceptance, e.s_ratio, e.t_ratio)

# get the connected component ids containing the slack bus
function get_slack_component_ids(g::Graphs.AbstractGraph{Bus,Line})
	cs = connected_components(g)
	for c in cs
		if in(3,Int[v.bus_type for v in c])
			return Int[v.id for v in c]
		end
	end
end

#function find_connected_graph(vs::Vector{Bus}, es::Vector{Line})
#	n = length(vs)
#	A = zeros(Int64, n,n)
#    	for edge in es
#			if edge.line_status
#				A[edge.source.id, edge.target.id] = 1
#				A[edge.target.id, edge.source.id] = 1
#			end
#    	end
#	queue = filter(v -> v.bus_type == 3, vs)[1].id
#	connected = []
#	while(length(queue) != 0)
#        connected  = [connected; queue[1]]
#	x = findin(A[:,queue[1]],1)
#        queue = [queue; setdiff(x, connected)]
#        queue = unique(queue[2:end])
#    end
#	#is_connected = falses(n,1)
#	#is_connected[connected] = true
#	return sort(connected)
#end


