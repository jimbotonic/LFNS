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
		"--ssolver"
			help = "sub solver to be used (e.g., NR, RK, ...)"
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
		"--u_fn"
			help = "file name for uniform random distribution in [-1,1]"
			required = false
		"--bt_fn"
			help = "file name of the boostrap T vector"
			required = false
		"--iter"
			help = "iteration number"
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
	export_graphml(g, "my_export.graphml")

	o_args = Dict{Symbol,Any}()
	o_args[:g] = g
	o_args[:bootstrap_iter] = 0
	s = Simulator(g,NR_solver,o_args,100.,1e-8,50)
	
	# launch the simulation
	state = simulation(s)
	
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
	state = simulation(s)

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
	state = simulation(s)

	export_csv_data(state.T, "T_out.csv")
elseif solver == "KR"
	function get_state(s::Simulator,A::Array{Float64,1},alpha::Float64,max_value::Float64)
		P = init_P1(A,alpha,max_value)
		change_P(s.g,P)
		return simulation(s)
	end

	g = load_serialized(pargs["g_fn"])
	n = length(vertices(g))
	
	# initialize simulation parameters
	sb = 1.
	max_iter = round(Int64,1e5)
	
	epsilon = 1e-6
	# max degree of Eurogrid is 17
	max_value = 1.

	ssolver = pargs["ssolver"]
	if ssolver == "SD"
		o_args = Dict{Symbol,Any}()
		o_args[:d] = 1
		s = Simulator(g,SD_solver,o_args,sb,epsilon,max_iter)
	elseif ssolver == "RK"
		o_args = Dict{Symbol,Any}()
		o_args[:h] = 3e-2
		s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)
	end
	
	#U = init_unif_dist(n)
	#export_csv_data(U, "U.csv")
	u_fn = pargs["u_fn"]
	U = collect(load_csv_data(u_fn)[1])
	iter = parse(Int,pargs["iter"])
	u_name = basename(u_fn)[1:end-4]
	
	# if we have an initial bootstraping vector T
	bt_fn = pargs["bt_fn"]
	if typeof(bt_fn) != Void
		bT = collect(load_csv_data(bt_fn)[1])
		change_T(s.g,bT)
	end

	states = Array{Array{Float64,1},1}()

	for j in 1:iter
		# start at 0 (ie, "flat start")
		t = (j-1)/1000
		state = get_state(s,U,t,max_value)
		
		@debug("# simulation $j (alpha=$t)")
		@debug("----------")

		push!(states,state.T)
		#export_fn = "T_" * "$iter" * "_" * "$ssolver" * "_" * "$epsilon" * "_" * "$u_name" * ".csv"
    		#export_csv_data(state.T, export_fn)
	end

	serialize_to_file(states, "states_$max_value.jld")
end
