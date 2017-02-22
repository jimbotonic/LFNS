using Logging, DataFrames, JLD

include("graphs.jl")

@Logging.configure(level=ERROR)

# ENTSOE standard voltage levels
ENTSOE_VOLTAGE_LEVELS = Dict(0=> 750., 1=> 380., 2=> 220., 3=> 150., 4=> 120., 5=> 110., 6=> 70., 7=> 27., 8=> 330., 9=> 500.)

# usual base power value
S_BASE = 100.

# load ENTSOE file
function load_ENTSOE(filename::AbstractString)
	# load file
	f = open(filename, "r")
	lines = readlines(f)
	close(f)

	vs = Bus[]
	es = Line[]

	node_section = "##ZCH"
	edge_section1 = "##L"
	edge_section2 = "##T"
	edge_section3 = "##R"
	end_section = "##C"
	
	in_node_section = false
	in_edge_section1 = false
	in_edge_section2 = false
	in_edge_section3 = false

	# incremented node id
	nid = 1
	# incremented edge id
	eid = 1
	# number of standard power line
	nline = 0
	# transformer id
	tid = 1
	# node name -> node
	name_node = Dict{AbstractString,Bus}()
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
				
				v = Bus(nid, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, P_min, P_max, sh_conductance, sh_susceptance)
				push!(vs, v)
				@info("adding vertex: ", v)

				# complete utils dictionaries and increment node id
				name_node[name] = v
				id_bv[nid] = base_voltage
				nid += 1
			end
		elseif in_edge_section1
			name1 = strip(l[1:8])
			name2 = strip(l[10:17])
			line_source = name_node[name1]
			line_target = name_node[name2]
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
			per_unit = S_BASE / id_bv[line_source.id]^2 
			resistance = float(strip(replace(l[23:28],',','.'))) * per_unit
			reactance = float(strip(replace(l[30:35],',','.'))) * per_unit
			admittance = 1./(resistance+reactance*im)
			# value is given in micro Siemens
			sh_susceptance = 1e-6*float(strip(replace(l[37:44],',','.'))) / per_unit
			# use default value
			s_ratio = 1.
			t_ratio = 1.

			edge = Line(eid, line_source, line_target, line_type, line_status, admittance, sh_susceptance, complex(s_ratio), complex(t_ratio))
			push!(es, edge)
			@info("adding edge: ", edge)

			# increment edge id
			eid += 1
			# count the number of standard power line
			nline += 1
		elseif in_edge_section2
			name1 = strip(l[1:8])
			name2 = strip(l[10:17])
			line_source = name_node[name1]
			line_target = name_node[name2]
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
			per_unit = S_BASE / id_bv[line_source.id]^2 
			resistance = float(strip(replace(l[41:46],',','.'))) * per_unit
			reactance = float(strip(replace(l[48:53],',','.'))) * per_unit
			admittance = 1./(resistance+reactance*im)
			sh_susceptance = 1e-6*float(strip(replace(l[55:62],',','.'))) / per_unit
			ratio1 = float(strip(replace(l[23:27],',','.')))
			ratio2 = float(strip(replace(l[29:33],',','.')))
			# t_ratio in pu
			t_ratio = (ratio1/ratio2)*(vs[line_target.id].base_voltage/vs[line_source.id].base_voltage)
			s_ratio = 1.

			edge = Line(eid, line_source, line_target, line_type, line_status, admittance, sh_susceptance, complex(s_ratio), complex(t_ratio))
			push!(es, edge)
			@info("adding edge: ", edge)

			# transformer name
			t_name = l[1:19]
			# complete t_name -> transformer id dictionary
			trans_name_id[t_name] = tid
			# increment edge and transformer id
			eid += 1
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
					# theta is the regulation angle
					theta = 0.
				else
					n_p = float(strip(replace(l[55:57],',','.')))
					du = float(strip(replace(l[40:44],',','.')))
					theta = (pi/180.0)*float(strip(replace(l[46:50],',','.')))		
				end
			else
				n_p = 0
				du = 1
				theta = 0
			end
			# rho is the amplitude of the regulation factor
			rho = 1/sqrt((0.01*n_p*du*sin(theta))^2+(1+0.01*n_p*du*cos(theta))^2)
			# alpha is the phase of the regulation factor
			alpha = atan(0.01*n_p*du*sin(theta)/(1+0.01*n_p*du*cos(theta)))
			# change the transformation ratio of the transformers with regulation
			es[trans_name_id[t_name]+nline].t_ratio *= rho*exp(-alpha*im)
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

	return graph(vs, es, is_directed=false)
end

# load IEEE Solved Load Flow Data file
function load_IEEE_SLFD(filename::AbstractString)
	# load file
	f = open(filename, "r")
	lines = readlines(f)
	close(f)

	vs = Bus[]
	es = Line[]
	
	# incremented edge id
	eid = 1

	node_section = "BUS DATA FOLLOWS"
	edge_section = "BRANCH DATA FOLLOWS"
	end_section = "-999"

	in_node_section = false
	in_edge_section = false
	
	# node id -> node
	id_node = Dict{Int64,Bus}()

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
			# convert angles defined in degrees into radians
			angle = float(strip(replace(l[34:40],',','.')))/(180*pi)
			# use default value
			base_voltage = S_BASE^2
			
			active_load = float(strip(replace(l[41:49],',','.')))
			reactive_load = float(strip(replace(l[50:58],',','.')))
			load = complex(active_load, reactive_load)

			active_generation = float(strip(replace(l[59:67],',','.')))
			reactive_generation = float(strip(replace(l[68:75],',','.')))
			generation = complex(active_generation, reactive_generation)
			
			if bus_type == 0
				# for PQ bus, bus voltage is set to its base voltage
				init_voltage = 1.
			else
				init_voltage = float(strip(replace(l[85:90],',','.')))
			end
			Q_min = float(strip(replace(l[91:98],',','.')))
			Q_max = float(strip(replace(l[99:106],',','.')))
			sh_conductance = float(strip(replace(l[107:114],',','.')))
			sh_susceptance = float(strip(replace(l[115:122],',','.')))
			
			v = Bus(id, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, 0., 0., sh_conductance, sh_susceptance,0.,0.)
			push!(vs, v)
			@info("adding vertex: ", v)
			
			# complete id -> node dictionary
			id_node[id] = v
		elseif in_edge_section
			source_id = parse(Int, strip(l[1:5]))
			target_id = parse(Int, strip(l[6:11]))
			line_type = parse(Int, strip(l[19:19]))
			# NB: per_unit = 1 here
			resistance = float(strip(replace(l[22:29],',','.')))
			reactance = float(strip(replace(l[32:39],',','.')))
			admittance = 1./(resistance+reactance*im)
			sh_susceptance = float(strip(replace(l[43:50],',','.')))
			rtemp = float(strip(replace(l[77:83],',','.')))
			phase_shift=float(strip(replace(l[84:90],',','.')))*(pi/180)
			if rtemp > 0 
				s_ratio = (1./rtemp)*exp(im*phase_shift)
			else
				s_ratio = 1.*exp(im*phase_shift)
			end
			t_ratio = 1.

			edge = Line(eid, id_node[source_id], id_node[target_id], line_type, true, admittance, sh_susceptance, complex(s_ratio), complex(t_ratio))
			push!(es, edge)
			@info("adding edge: ", edge)

			# increment edge id
			eid += 1
		end

		if startswith(l, node_section)
			in_node_section = true
		elseif startswith(l, edge_section)
			in_edge_section = true
		end
	end	

	return graph(vs, es, is_directed=false)
end


# load Matpower Solved Load Flow Data file
function load_MATPOWER_SLFD(filename::AbstractString)
	# load file
	f = open(filename, "r")
	lines = readlines(f)
	close(f)

	vs = Bus[]
	es = Line[]
	
	# incremented edge id
	eid = 1

	# incremented node id
	nid = 1	

	node_section = "mpc.bus ="
	gen_section = "mpc.gen ="
	edge_section = "mpc.branch ="
	end_section = "];"

	in_node_section = false
	in_gen_section = false
	in_edge_section = false
	
	# node id -> node
	id_node = Dict{AbstractString,Bus}()

	for l in lines	
		if startswith(l, end_section)
			in_node_section = false
			in_gen_section = false
			in_edge_section = false
		end
		
		if in_node_section
			data = split(l)
			name = data[1]

			bus_type = parse(Int64, data[2])
			final_voltage = parse(Float64, data[8])
			# convert angles defined in degrees into radians
			angle = parse(Float64,data[9])*pi/180
			# base voltage value in kV
			base_voltage = parse(Float64,data[10])
			
			active_load = parse(Float64,data[3])
			reactive_load = parse(Float64,data[4])
			load = complex(active_load, reactive_load)

			# Generation is taken into account in the generation_section, 
			active_generation = 0.0
			reactive_generation = 0.0
			generation = complex(active_generation, reactive_generation)
			Q_max=0.0
			Q_min=0.0
			P_max=0.0
			P_min=0.0

			if bus_type > 1
				init_voltage = 1. # The initial voltage is set in the generator section
			else
				# for PQ bus, bus voltage is set to its base voltage
				bus_type=0				
				init_voltage = 1.
			end
			sh_conductance = parse(Float64,data[5])/S_BASE
			sh_susceptance = parse(Float64,data[6])/S_BASE
			
			v = Bus(nid, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, P_min, P_max, sh_conductance, sh_susceptance,0.,0.)
			push!(vs, v)
			@info("adding vertex: ", v)
			
			# increment node id
			nid += 1

			# complete id -> node dictionary			
			id_node[name] = v

		elseif in_gen_section
			data = split(l)
			name = data[1]
			
			v=id_node[name]
			active_generation = parse(Float64,data[2])
			reactive_generation = parse(Float64,data[3])
			generation = complex(active_generation, reactive_generation)

			v.generation=generation					
				
			v.init_voltage=parse(Float64,data[6])
			v.Q_max=parse(Float64,data[4])
			v.Q_min=parse(Float64,data[5])
			v.P_max=parse(Float64,data[9])
			v.P_min=parse(Float64,data[10])

		elseif in_edge_section
			data = split(l)
			source_id = data[1]
			target_id = data[2]
			resistance = parse(Float64, data[3]) 
			reactance =  parse(Float64, data[4]) 
			admittance = 1./(resistance+reactance*im)
			sh_susceptance = parse(Float64, data[5]) 
			rtemp = parse(Float64, data[9])
			phase_shift= parse(Float64, data[10])*(pi/180)
			line_status = convert(Bool,parse(Int,data[11]))

			#normal lines, transformers and phase shifters
			if rtemp == 0. && phase_shift == 0.
				line_type=0
				s_ratio=1.
			elseif rtemp > 0.
				line_type=1
				s_ratio = (1./rtemp)*exp(-im*phase_shift)
			else 
				line_type=1
				s_ratio = 1.*exp(-im*phase_shift)
			end
			t_ratio = 1.

			edge = Line(eid, id_node[source_id], id_node[target_id], line_type, line_status, admittance, sh_susceptance, complex(s_ratio), complex(t_ratio))
			push!(es, edge)
			@info("adding edge: ", edge)

			# increment edge id
			eid += 1
		end

		if startswith(l, node_section)
			in_node_section = true
		elseif startswith(l, gen_section)
			in_gen_section = true
		elseif startswith(l, edge_section)
			in_edge_section = true
		end
	end	

	return graph(vs, es, is_directed=false)
end


# export graph to graphml
function export_graphml(g::Graphs.AbstractGraph{Bus,Line}, filename::AbstractString)
	vs = vertices(g)
	es = edges(g)
	graphmlFile = open(filename, "w")

	write(graphmlFile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
	write(graphmlFile, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n")
	write(graphmlFile, "	xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n")  
	write(graphmlFile, "	xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n")
        write(graphmlFile, "	http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n")
	write(graphmlFile, "<graph id=\"G\" edgedefault=\"undirected\">\n")

	# strip non-standard characters
	# stripc0{T<:AbstractString}(a::T) = replace(a, r"[\x00-\x1f\x7f]", "")
	stripc0x{T<:AbstractString}(a::T) = replace(a, r"[^\x20-\x7e]", "")

	for n in vs
		write(graphmlFile, "	<node id=\"" * string(n.id) * "\">\n")
		write(graphmlFile, "		<data key=\"name\">" * stripc0x(n.name) * "</data>\n")
		write(graphmlFile, "		<data key=\"bus_type\">" * string(n.bus_type) * "</data>\n")
		write(graphmlFile, "		<data key=\"init_voltage\">" * string(n.init_voltage) * "</data>\n")
		write(graphmlFile, "		<data key=\"final_voltage\">" * string(n.final_voltage) * "</data>\n")
		write(graphmlFile, "		<data key=\"base_voltage\">" * string(n.base_voltage) * "</data>\n")
		write(graphmlFile, "		<data key=\"angle\">" * string(n.angle) * "</data>\n")
		write(graphmlFile, "		<data key=\"load\">" * string(n.load) * "</data>\n")
		write(graphmlFile, "		<data key=\"generation\">" * string(n.generation) * "</data>\n")
		write(graphmlFile, "		<data key=\"Q_min\">" * string(n.Q_min) * "</data>\n")
		write(graphmlFile, "		<data key=\"Q_max\">" * string(n.Q_max) * "</data>\n")
		write(graphmlFile, "		<data key=\"P_min\">" * string(n.P_min) * "</data>\n")
		write(graphmlFile, "		<data key=\"P_max\">" * string(n.P_max) * "</data>\n")
		write(graphmlFile, "		<data key=\"sh_conductance\">" * string(n.sh_conductance) * "</data>\n")
		write(graphmlFile, "		<data key=\"sh_susceptance\">" * string(n.sh_susceptance) * "</data>\n")
		write(graphmlFile, "		<data key=\"lat\">" * string(n.lat) * "</data>\n")
		write(graphmlFile, "		<data key=\"lng\">" * string(n.lng) * "</data>\n")
		write(graphmlFile, "	</node>\n")
	end

	for edge in es
		write(graphmlFile, "	<edge source=\"" * string(edge.source.id) * "\" target=\"" * string(edge.target.id) * "\">\n")
		write(graphmlFile, "		<data key=\"line_type\">" * string(edge.line_type) * "</data>\n")
		if edge.line_status
			write(graphmlFile, "		<data key=\"line_status\">1</data>\n")
		else
			write(graphmlFile, "		<data key=\"line_status\">0</data>\n")
		end
		write(graphmlFile, "		<data key=\"admittance\">" * string(edge.admittance) * "</data>\n")
		write(graphmlFile, "		<data key=\"sh_susceptance\">" * string(edge.sh_susceptance) * "</data>\n")
		write(graphmlFile, "		<data key=\"s_ratio\">" * string(edge.s_ratio) * "</data>\n")
		write(graphmlFile, "		<data key=\"t_ratio\">" * string(edge.t_ratio) * "</data>\n")
		write(graphmlFile, "	</edge>\n")
	end

	write(graphmlFile, "</graph>\n")
	write(graphmlFile, "</graphml>\n")
	close(graphmlFile)									
end

# export graph topology to graphml
function export_topology_graphml(g::Graphs.AbstractGraph{Bus,Line}, filename::AbstractString)
	vs = vertices(g)
	es = edges(g)
	graphmlFile = open(filename, "w")

	write(graphmlFile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
	write(graphmlFile, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n")
	write(graphmlFile, "	xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n")  
	write(graphmlFile, "	xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n")
        write(graphmlFile, "	http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n")
	write(graphmlFile, "<graph id=\"G\" edgedefault=\"undirected\">\n")

	for n in vs
		write(graphmlFile, "	<node id=\"" * string(n.id) * "\"/>\n")
	end

	for edge in es
		write(graphmlFile, "	<edge source=\"" * string(edge.source.id) * "\" target=\"" * string(edge.target.id) * "\"/>\n")
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
# NB: dimensions are (# rows, # columns)
function load_csv_data(fn::AbstractString)
	return readtable(fn, header = false)
end

# load serialized adjacency list
#
# @deprecated
function load_serialized(filename::AbstractString)
	x = open(filename, "r") do file
		deserialize(file)
	end
	return x
end

# serialize graph
#
# @deprecate
function serialize_to_file(x::Any, filename::AbstractString)
	open(filename, "w") do file
		serialize(file, x)
	end
end			

# load serialized JLS data
function load_jls_serialized(filename::AbstractString)
	x = open(filename, "r") do file
		deserialize(file)
	end
	return x
end

# serialize data to JLS format
function serialize_to_jls(x::Any, filename::AbstractString)
	open("$filename.jls", "w") do file
		serialize(file, x)
	end
end

# load serialized JLD data
#
# NB: better long-term backwards compatibility than jls files 
function load_jld_serialized(name::AbstractString, filename::AbstractString)
	x = jldopen(filename, "r") do file
    		read(file, name)
  	end
	return x
end

# serialize data to JLD format
#
# NB: better long-term backwards compatibility than jls files 
function serialize_to_jld(x::Any, name::AbstractString, filename::AbstractString)
	jldopen("$filename.jld", "w") do file
		write(file, name, x)
	end
end	

# load graph from the specified Y and P files
#
# CSV files with no-header and comma-separated are expected
# P_fn contains one float per line
# Y_fn: node_id1, node_id2, G_value (0 in the non-dissipative case), B value
function load_graph(P_fn::AbstractString, Y_fn::AbstractString)
	# load the data from Y and P file
	P_df = load_csv_data(P_fn)
	Y_df = load_csv_data(Y_fn)

	# the size of the network is assumed to be the # of rows in P0 
	P = collect(P_df[1])
	n = length(P)

	vertices = Bus[]
	for i in 1:n
		push!(vertices, Bus(i, 0., Float64(P[i])))
	end
	
	edges = Line[]
	for i in 1:size(Y_df,1)
		source = vertices[Y_df[i,1]]
		target = vertices[Y_df[i,2]]
		y = Float64(Y_df[i,3])+Float64(Y_df[i,4])*im
		push!(edges, Line(i, source, target, y))
	end
	
	return graph(vertices, edges, is_directed=false)
end


# load graph from the specified Y, T and P files
#
# CSV files with no-header and comma-separated are expected
# P_fn contains one float per line
# T_fn contains on float per line
# Y_fn: node_id1, node_id2, G_value (0 in the non-dissipative case), B value
function load_graph2(P_fn::AbstractString, T_fn::AbstractString, Y_fn::AbstractString)
	# load the data from Y and P file
	P_df = load_csv_data(P_fn)
	T_df = load_csv_data(t_fn)
	Y_df = load_csv_data(Y_fn)

	# the size of the network is assumed to be the # of rows in P0 
	P = collect(P_df[1])
	T = collect(T_df[1])
	n = length(P)

	vertices = Bus[]
	for i in 1:n
		push!(vertices, Bus(i, Float64(T[i]), Float64(P[i])))
	end
	
	edges = Line[]
	for i in 1:size(Y_df,1)
		source = vertices[Y_df[i,1]]
		target = vertices[Y_df[i,2]]
		y = Float64(Y_df[i,3])+Float64(Y_df[i,4])*im
		push!(edges, Line(i, source, target, y))
	end
	
	return graph(vertices, edges, is_directed=false)
end

