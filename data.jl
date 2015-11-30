using Logging, DataFrames

@Logging.configure(level=CRITICAL)

# ENTSOE standard voltage levels
ENTSOE_VOLTAGE_LEVELS = Dict(0=> 750., 1=> 380., 2=> 220., 3=> 150., 4=> 120., 5=> 110., 6=> 70., 7=> 27., 8=> 330., 9=> 500.)

# usual base power value
S_BASE = 100.

# node type
type Node
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
type Edge
	source::Node
	target::Node
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


# load ENTSOE file
function load_ENTSOE(filename::AbstractString)
	# load file
	f = open(filename, "r")
	lines = readlines(f)
	close(f)

	nodes = Node[]
	edges = Edge[]

	node_section = "##ZCH"
	edge_section1 = "##L"
	edge_section2 = "##T"
	edge_section3 = "##R"
	end_section = "##C"
	
	in_node_section = false
	in_edge_section1 = false
	in_edge_section2 = false
	in_edge_section3 = false

	# auto-incremented node id
	nid = 1
	# number of standard power line
	nline = 0
	# transformer id
	tid = 1
	# node name -> node
	name_node = Dict{AbstractString,Node}()
	# transformer name -> auto incremented id
	trans_name_id = Dict{AbstractString,Int64}()
	# node id -> base voltage (used to compute the per-unit values)
	id_bv = Dict{Int64, Float64}()
	
	for l in lines
		if startswith(l, end_section)
			in_node_section = false
			in_edge_section1 = false
			in_edge_section2 = false
			in_edge_section3 = false
		end

		if in_node_section
			if !startswith(l,'#')
				name = strip(l[1:8])
				voltage_level = parse(Int, strip(l[7:7]))
				base_voltage = ENTSOE_VOLTAGE_LEVELS[voltage_level]
				bus_type = parse(Int, strip(l[25:25]))
				# when disconnected final_voltage might empty
				if length(strip(l[27:32])) != 0
					final_voltage = float(strip(replace(l[27:32],',','.')))
				else
					final_voltage = 0.
				end
				angle = 0.

				active_load = float(strip(replace(l[34:40],',','.')))
				reactive_load = float(strip(replace(l[42:48],',','.')))
				load = complex(active_load, reactive_load)

				active_generation = float(strip(replace(l[50:56],',','.')))
				reactive_generation = float(strip(replace(l[58:64],',','.')))
				generation = complex(active_generation, reactive_generation)

				init_voltage = final_voltage/base_voltage

				# sometimes the line is ended before the P/Q limits
				if length(l) > 85
					Q_min = float(strip(replace(l[82:88],',','.')))
					Q_max = float(strip(replace(l[90:96],',','.')))
					P_min = float(strip(replace(l[66:72],',','.')))
					P_max = float(strip(replace(l[74:80],',','.')))
				else
					Q_min = 0.
					Q_max = 0.
					P_min = 0.
					P_max = 0.
				end					
				sh_conductance = 0.
				sh_susceptance = 0.
				
				n = Node(nid, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, P_min, P_max, sh_conductance, sh_susceptance)
				push!(nodes, n)
				@info("adding node: ", n)

				# complete utils dictionaries and increment node id
				name_node[name] = n
				id_bv[nid] = base_voltage
				nid += 1
			end
		elseif in_edge_section1
			name1 = strip(l[1:8])
			name2 = strip(l[10:17])
			source = name_node[name1]
			target = name_node[name2]
			# bloc 1 -> 0
			line_type = 0
			# line status 0,1,2 -> ON 7,8,9 -> OFF
			lsc = parse(Int, strip(l[21:21]))
			if lsc == 0 || lsc == 1 || lsc == 2
				line_status = true 
			elseif lsc == 7 || lsc == 8 || lsc == 9
				line_status = false 
			end
			# compute per_unit value
			per_unit = S_BASE / id_bv[source.id]^2 
			resistance = float(strip(replace(l[23:28],',','.'))) * per_unit
			reactance = float(strip(replace(l[30:35],',','.'))) * per_unit
			# value is given in micro Siemens
			sh_susceptance = 1e-6*float(strip(replace(l[37:44],',','.'))) / per_unit
			# use default value
			s_ratio = 1.
			t_ratio = 1.

			edge = Edge(source, target, line_type, line_status, resistance, reactance, sh_susceptance, s_ratio, t_ratio)
			push!(edges, edge)
			@info("adding edge: ", edge)

			# count the number of standard power line
			nline += 1
		elseif in_edge_section2
			name1 = strip(l[1:8])
			name2 = strip(l[10:17])
			source = name_node[name1]
			target = name_node[name2]
			# bloc 2 -> 1
			line_type = 1
			# line status 0,1,2 -> ON 7,8,9 -> OFF
			lsc = parse(Int, strip(l[21:21]))
			if lsc == 0 || lsc == 1 || lsc == 2
				line_status = true 
			elseif lsc == 7 || lsc == 8 || lsc == 9
				line_status = false 
			end
			# compute per_unit value
			per_unit = S_BASE / id_bv[source.id]^2 
			resistance = float(strip(replace(l[41:46],',','.'))) * per_unit
			reactance = float(strip(replace(l[48:53],',','.'))) * per_unit
			sh_susceptance = 1e-6*float(strip(replace(l[55:62],',','.'))) / per_unit
			ratio1 = float(strip(replace(l[23:27],',','.')))
			ratio2 = float(strip(replace(l[29:33],',','.')))
			# t_ratio in pu
			t_ratio = (ratio1/ratio2)*(nodes[target.id].base_voltage/nodes[source.id].base_voltage)
			s_ratio = 1.

			edge = Edge(source, target, line_type, line_status, resistance, reactance, sh_susceptance, s_ratio, t_ratio)
			push!(edges, edge)
			@info("adding edge: ", edge)

			# transformer name
			t_name = l[1:19]
			# complete t_name -> transformer id dictionary
			trans_name_id[t_name] = tid
			# increment transformer id
			tid += 1
		elseif in_edge_section3
			# for more information about parameters used in the following,
			# see Quality of Datasetes and Calculation published by ENSTO-E.
			t_name = l[1:19]
			# see if there is something to read between cols 21 and 25
			# if yes -> phase regulation
			# if not -> angle regulation
			if length(strip(l[21:50])) != 0
				selection_criterion = length(strip(l[21:25]))
				if selection_criterion != 0
					# n_p is the current tap position
					n_p = float(strip(l[30:32]))
					# du is the voltage change per tap in percent
					du = float(strip(replace(l[21:25],',','.')))
					# Theta is the regulation angle
					Theta = 0.
				else
					n_p = float(strip(replace(l[55:57],',','.')))
					du = float(strip(replace(l[40:44],',','.')))
					Theta = (pi/180.0)*float(strip(replace(l[46:50],',','.')))		
				end
			else
				n_p = 0
				du = 1
				Theta = 0
			end
			# rho is the amplitude of the regulation factor
			rho = 1/sqrt((0.01*n_p*du*sin(Theta))^2+(1+0.01*n_p*du*cos(Theta))^2)
			# alpha is the phase of the regulation factor
			alpha = atan(0.01*n_p*du*sin(Theta)/(1+0.01*n_p*du*cos(Theta)))
			# change the transformation ratio of the transformers with regulation
			edges[trans_name_id[t_name]+nline].t_ratio *= rho*exp(-alpha*im)
		end

		if startswith(l, node_section)
			in_node_section = true
		elseif startswith(l, edge_section1)
			in_edge_section1 = true
		elseif startswith(l, edge_section2)
			in_edge_section2 = true
		elseif startswith(l, edge_section3)
			in_edge_section3 = true
		end
	end

	return nodes,edges
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
	
	# node id -> node
	id_node = Dict{Int64,Node}()

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
			# use default value
			base_voltage = S_BASE^2
			
			active_load = float(strip(replace(l[41:49],',','.')))
			reactive_load = float(strip(replace(l[50:58],',','.')))
			load = complex(active_load, reactive_load)

			active_generation = -float(strip(replace(l[59:67],',','.')))
			reactive_generation = -float(strip(replace(l[68:75],',','.')))
			generation = complex(active_generation, reactive_generation)

			init_voltage = float(strip(replace(l[85:90],',','.')))
			Q_min = float(strip(replace(l[91:98],',','.')))
			Q_max = float(strip(replace(l[99:106],',','.')))
			sh_conductance = float(strip(replace(l[107:114],',','.')))
			sh_susceptance = float(strip(replace(l[115:122],',','.')))
			
			n = Node(id, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, 0., 0., sh_conductance, sh_susceptance)
			push!(nodes, n)
			@info("adding node: ", n)
			
			# complete id -> node dictionary
			id_node[id] = n
		elseif in_edge_section
			source_id = parse(Int, strip(l[1:5]))
			target_id = parse(Int, strip(l[6:11]))
			line_type = parse(Int, strip(l[19:19]))
			# NB: per_unit = 1 here
			resistance = float(strip(replace(l[22:29],',','.')))
			reactance = float(strip(replace(l[32:39],',','.')))
			sh_susceptance = float(strip(replace(l[43:50],',','.')))
			rtemp = float(strip(replace(l[77:83],',','.')))
			if rtemp > 0 
				s_ratio = 1./rtemp
			else
				s_ratio = 1.
			end
			t_ratio = 1.

			edge = Edge(id_node[source_id], id_node[target_id], line_type, true, resistance, reactance, sh_susceptance, s_ratio, t_ratio)
			push!(edges, edge)
			@info("adding edge: ", edge)
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
		write(graphmlFile, "	<data key=\"init_voltage\">" * string(n.init_voltage) * "</data>\n")
		write(graphmlFile, "	<data key=\"final_voltage\">" * string(n.final_voltage) * "</data>\n")
		write(graphmlFile, "	<data key=\"base_voltage\">" * string(n.base_voltage) * "</data>\n")
		write(graphmlFile, "	<data key=\"angle\">" * string(n.angle) * "</data>\n")
		write(graphmlFile, "	<data key=\"load\">" * string(n.load) * "</data>\n")
		write(graphmlFile, "	<data key=\"generation\">" * string(n.generation) * "</data>\n")
		write(graphmlFile, "	<data key=\"Q_min\">" * string(n.Q_min) * "</data>\n")
		write(graphmlFile, "	<data key=\"Q_max\">" * string(n.Q_max) * "</data>\n")
		write(graphmlFile, "	<data key=\"P_min\">" * string(n.P_min) * "</data>\n")
		write(graphmlFile, "	<data key=\"P_max\">" * string(n.P_max) * "</data>\n")
		write(graphmlFile, "	<data key=\"sh_conductance\">" * string(n.sh_conductance) * "</data>\n")
		write(graphmlFile, "	<data key=\"sh_susceptance\">" * string(n.sh_susceptance) * "</data>\n")
		write(graphmlFile, "</node>\n")
	end

	for edge in edges
		write(graphmlFile, "<edge id=\"" * string(edge.source.id) *"|" * string(edge.target.id) *"\" source=\"" * string(edge.source.id) * "\" target=\"" * string(edge.target.id) * "\">\n")
		write(graphmlFile, "	<data key=\"line_type\">" * string(edge.line_type) * "</data>\n")
		if edge.line_status
			write(graphmlFile, "	<data key=\"line_status\">1</data>\n")
		else
			write(graphmlFile, "	<data key=\"line_status\">0</data>\n")
		end
		write(graphmlFile, "	<data key=\"resistance\">" * string(edge.resistance) * "</data>\n")
		write(graphmlFile, "	<data key=\"reactance\">" * string(edge.reactance) * "</data>\n")
		write(graphmlFile, "	<data key=\"sh_susceptance\">" * string(edge.sh_susceptance) * "</data>\n")
		write(graphmlFile, "	<data key=\"s_ratio\">" * string(edge.s_ratio) * "</data>\n")
		write(graphmlFile, "	<data key=\"t_ratio\">" * string(edge.t_ratio) * "</data>\n")
		write(graphmlFile, "</edge>\n")
	end

	write(graphmlFile, "</graph>\n")
	write(graphmlFile, "</graphml>\n")
	close(graphmlFile)									
end

# write array to CSV file (with no header)
function export_csv_data(M, fn::AbstractString)
	# M is a vector, i.e., Array{Float64,1}
	if length(size(M)) == 1
		writetable(fn, DataFrame(a = M), header=false)
	# M is a matrix
	elseif length(size(M)) == 2
		writetable(fn, DataFrame(M), header=false)
	end
end

# load CSV data from file
# 
# NB: file is assumed to be in CSV format with no header
function load_csv_data(fn::AbstractString)
	return readtable(fn, header = false)
end

