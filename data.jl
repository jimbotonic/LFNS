# node type
type Node
	id::Int64
	name::AbstractString
	bus_type::Int
	final_voltage::Float64
	angle::Float64
	load::Complex{Float64}
	generation::Complex{Float64}
	init_voltage::Float64
	Q_min::Float64
	Q_max::Float64
	sh_conductance::Float64
	sh_susceptance::Float64
end

# edge type
type Edge
	source_id::Int64
	target_id::Int64
	line_type::Int
	resistance::Float64
	reactance::Float64
	sh_susceptance::Float64
	ratio::Float64
end

# load IEEE Solved Load Flow Data file
function load_IEEE_SLFD(filename::AbstractString)
	# load file
	f = open(filename, "r")
	lines = readlines(f)
	close(f)

	nodes = Node[]
	edges = Edge[]

	node_section = "BUS DATA FOLLOWS"
	edge_section = "BRANCH DATA FOLLOWS"
	end_section = "-999"

	in_node_section = false
	in_edge_section = false

	for l in lines
		if startswith(l, end_section)
			in_node_section = false
			in_edge_section = false
		end

		if in_node_section
			id = parse(Int, strip(l[1:5]))
			name = strip(l[6:15])
			bus_type = parse(Int, strip(l[26:26]))
			final_voltage = float(strip(replace(l[28:33],',','.')))
			angle = float(strip(replace(l[34:40],',','.')))
			
			active_load = float(strip(replace(l[41:49],',','.')))
			reactive_load = float(strip(replace(l[50:58],',','.')))
			load = complex(active_load, reactive_load)

			active_generation = float(strip(replace(l[59:67],',','.')))
			reactive_generation = float(strip(replace(l[68:75],',','.')))
			generation = complex(active_generation, reactive_generation)

			init_voltage = float(strip(replace(l[85:90],',','.')))
			Q_min = float(strip(replace(l[91:98],',','.')))
			Q_max = float(strip(replace(l[99:106],',','.')))
			sh_conductance = float(strip(replace(l[107:114],',','.')))
			sh_susceptance = float(strip(replace(l[115:122],',','.')))
			
			n = Node(id, name, bus_type, final_voltage, angle, load, generation, init_voltage, Q_min, Q_max, sh_conductance, sh_susceptance)
			push!(nodes, n)
			println(n)
		elseif in_edge_section
			source_id = parse(Int, strip(l[1:5]))
			target_id = parse(Int, strip(l[6:11]))
			line_type = parse(Int, strip(l[19:19]))
			resistance = float(strip(replace(l[22:29],',','.')))
			reactance = float(strip(replace(l[32:39],',','.')))
			sh_susceptance = float(strip(replace(l[43:50],',','.')))
			ratio = float(strip(replace(l[77:83],',','.')))

			edge = Edge(source_id, target_id, line_type, resistance, reactance, sh_susceptance, ratio)
			push!(edges, edge)
			println(edge)
		end

		if startswith(l, node_section)
			in_node_section = true
		elseif startswith(l, edge_section)
			in_edge_section = true
		end
	end	

	return nodes,edges
end

# export data in graphml
function export_graphml(filename::AbstractString, nodes::Array{Node,1}, edges::Array{Edge,1})
	graphmlFile = open(filename, "w")

	write(graphmlFile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
	write(graphmlFile, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns/graphml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"  xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/graphml\">\n")
	write(graphmlFile, "<graph edgedefault=\"undirected\">\n")

	for n in nodes
		write(graphmlFile, "<node id=\"" * string(n.id) * "\">\n")
		write(graphmlFile, "	<data key=\"name\">" * n.name * "</data>\n")
		write(graphmlFile, "	<data key=\"bus_type\">" * string(n.bus_type) * "</data>\n")
		write(graphmlFile, "	<data key=\"final_voltage\">" * string(n.final_voltage) * "</data>\n")
		write(graphmlFile, "	<data key=\"angle\">" * string(n.angle) * "</data>\n")
		write(graphmlFile, "	<data key=\"load\">" * string(n.load) * "</data>\n")
		write(graphmlFile, "	<data key=\"generation\">" * string(n.generation) * "</data>\n")
		write(graphmlFile, "	<data key=\"init_voltage\">" * string(n.init_voltage) * "</data>\n")
		write(graphmlFile, "	<data key=\"Q_min\">" * string(n.Q_min) * "</data>\n")
		write(graphmlFile, "	<data key=\"Q_max\">" * string(n.Q_max) * "</data>\n")
		write(graphmlFile, "	<data key=\"sh_conductance\">" * string(n.sh_conductance) * "</data>\n")
		write(graphmlFile, "	<data key=\"sh_susceptance\">" * string(n.sh_susceptance) * "</data>\n")
		write(graphmlFile, "</node>\n")
	end

	for edge in edges
		write(graphmlFile, "<edge id=\"" * string(edge.source_id) *"|" * string(edge.target_id) *"\" source=\"" * string(edge.source_id) * "\" target=\"" * string(edge.target_id) * "\">\n")
		write(graphmlFile, "	<data key=\"line_type\">" * string(edge.line_type) * "</data>\n")
		write(graphmlFile, "	<data key=\"resistance\">" * string(edge.resistance) * "</data>\n")
		write(graphmlFile, "	<data key=\"reactance\">" * string(edge.reactance) * "</data>\n")
		write(graphmlFile, "	<data key=\"sh_susceptance\">" * string(edge.sh_susceptance) * "</data>\n")
		write(graphmlFile, "	<data key=\"ratio\">" * string(edge.ratio) * "</data>\n")
		write(graphmlFile, "</edge>\n")
	end

	write(graphmlFile, "</graph>\n")
	write(graphmlFile, "</graphml>\n")
	close(graphmlFile)									
end

nodes,edges = load_IEEE_SLFD(ARGS[1])
export_graphml("test.graphml", nodes, edges)

