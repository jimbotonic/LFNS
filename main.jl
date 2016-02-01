using ArgParse

include("init.jl")
include("data.jl")
include("solvers.jl")
include("graphs.jl")

# parse program arguments
function parse_cl()
	s = ArgParseSettings()
	@add_arg_table s begin
		"--solver"
			help = "solver to be used (e.g., NR, RK, ...)"
		"--ft"
			help = "file type to be loaded (e.g., IEEE, ENTSOE)"
			required = false
		"--fn"
			help = "file name"
			required = false
		"--p_fn"
			help = "file name for P vector"
			required = false
		"--y_fn"
			help = "file name for admittance matrix Y"
			required = false
		"--g_fn"
			help = "file name for serialized graph"
			required = false
	end
	return parse_args(s)
end

pargs = parse_cl()
solver = pargs["solver"]

if solver == "NR"
	ftype = pargs["ft"]
	fn = pargs["fn"]
	if ftype == "IEEE"
		# load file
		g = load_IEEE_SLFD(fn)
	elseif ftype == "ENTSOE"
		g = load_ENTSOE(fn)
	end
	
	vs = vertices(g)
	n = length(vs) 
	
	@info("# vertices: ", length(vs))
	@info("# edges: ", length(edges(g)))
	
	# export graph to graphml
	#export_graphml(g, "my_export.graphml")

	o_args = Dict{Symbol,Any}()
	o_args[:g] = g
	o_args[:bootstrap_iter] = 3
	s = Simulator(g,NR_solver,o_args,100.,1e-8,15)
	
	# launch the simulation
	simulation(s)
	state = s.states[1]
	
	export_csv_data(state.V, "V_out.csv")
	state.T = state.T*180/pi
	export_csv_data(state.T, "T_out.csv")
elseif solver == "RK"
	p_fn = pargs["p_fn"] # vector of initial powers 
	y_fn = pargs["y_fn"] # initial admittance matrix
	g = load_graph(p_fn,y_fn) # load Admittance matrix and injected/consumed powers

	o_args = Dict{Symbol,Any}()
	o_args[:h] = 1e-2
	s = Simulator(g,RK_solver1,o_args,1.,1e-11,round(Int64,1e5))
	
	# launch the simulation
	simulation(s)
	state = s.states[1]

	@info("# iteration: ", state.n_iter)
	export_csv_data(state.T, "T_out.csv")
	export_csv_data(state.Tdot, "Tdot_out.csv")
elseif solver == "SD"
	p_fn = pargs["p_fn"]  
	y_fn = pargs["y_fn"] 
	g = load_graph(p_fn,y_fn) 

	o_args = Dict{Symbol,Any}()
	o_args[:d] = 1e-2
	s = Simulator(g,SD_solver,o_args,1.,1e-6,round(Int64,1e5))
	export_csv_data(imag(s.Y), "B.csv")
	
	# launch the simulation
	simulation(s)
	state = s.states[1]

	export_csv_data(state.T, "T_out.csv")
elseif solver == "KR"
	g = load_serialized(pargs["g_fn"])
	n = length(vertices(g))
	sb = 1.
	max_iter = round(Int64,1e5)

	method = "SD"

	if method == "SD"
		o_args = Dict{Symbol,Any}()
		o_args[:d] = 1
		epsilon = 1e-8
		s = Simulator(g,SD_solver,o_args,sb,epsilon,max_iter)
	elseif method == "RK"
		o_args = Dict{Symbol,Any}()
		o_args[:h] = 3e-2
		epsilon = 1e-8
		s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)
	end
	
	max_degree = 17.
	alpha = 0.16
	A = init_unif_dist(n)
	counter = 1
	step = 1e-3
	for t in 0:step:alpha
#		t = 3e-3
		P = init_P1(A,t,max_degree)
#		writecsv("P.csv",P)
		change_P(s.g,P)
		simulation(s)
		@info("step $counter/", (alpha/step))
		counter += 1
#		state = s.states[1]
#		export_csv_data(state.T, "T_out.csv")
	end
end
