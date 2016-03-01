using ArgParse, Logging, ConfParser

include("init.jl")
include("data.jl")
include("solvers.jl")
include("graphs.jl")

# @Logging.configure(level=INFO)

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
			help = "iteration number (if parallel=yes : iteration number within the subinterval of alpha"
			required = false
		"--g_fn"
			help = "file name for serialized graph"
			required = false
		"--parallel"
			default = "1"
			help = "1: no parallelization | 2: subinterval | 3: parallelization"
			required = false
		"--nprocs"
			help = "# of processes to allocate for parallelization"
			required = false
		"--alpha_i"
			help = "first alpha value (parallel=2)"
			required = false
		"--alpha_interval_length"
			help = "length of alpha interval (parallel=2)"
			required = false
	end
	return parse_args(s)
end

# parse arguments
pargs = parse_cl()
solver = pargs["solver"]

# loading config file
conf = ConfParse("config.ini")
parse_conf!(conf)

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
	function get_state(s::Simulator,P_ref::Array{Float64,1},t::Float64)
		P = P_ref*t
		@debug("norm P (2/Inf): ", norm(P,2), "/", norm(P,Inf))
		change_P(s.g,P)
		return simulation(s)
	end
	
	g = load_serialized(pargs["g_fn"])
	n = length(vertices(g))
	
	# initialize simulation parameters
	sb = parse(Float64,retrieve(conf,"solvers","base_voltage"))
	max_iter = round(Int,parse(Float64,retrieve(conf,"solvers","max_iter")))
	epsilon = parse(Float64,retrieve(conf,"solvers","epsilon"))
	# NB: Eurogrid PC has 6021 nodes and a max degree of 17
	# critical value is between 0.24 and 0.255
	max_value = parse(Float64,retrieve(conf,"eurogrid","max_value"))

	ssolver = pargs["ssolver"]
	if ssolver == "SD"
		o_args = Dict{Symbol,Any}()
		o_args[:d] = parse(Int,retrieve(conf,"sd","d"))
		s = Simulator(g,SD_solver,o_args,sb,epsilon,max_iter)
	elseif ssolver == "RK"
		o_args = Dict{Symbol,Any}()
		o_args[:h] = parse(Float64,retrieve(conf,"rk","h"))
		s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)
	end

	# generate a uniform distribution
	#U = init_unif_dist(n)
	#export_csv_data(U, "U.csv")
	
	u_fn = pargs["u_fn"]
	U = collect(load_csv_data(u_fn)[1])
	iter = parse(Int,pargs["iter"])
	# remove file extension from base name
	u_name = basename(u_fn)[1:end-4]

	P_ref = init_P3(U)
	
	# if we have an initial bootstraping vector T
	bt_fn = pargs["bt_fn"]
	if typeof(bt_fn) != Void
		bT = collect(load_csv_data(bt_fn)[1])
		change_T(s.g,bT)
	end

	# array of computed states 
	states = Dict{Float64,State}()

	# parallel mode	
	par = parse(Int,pargs["parallel"])
	
	# no parallelization (simulation on the whole interval of alpha values)
	if par == 1
		tic()
		for j in 1:iter
			# start at 0 (ie, "flat start")
			alpha = (j-1)/iter	
			if j > 1
				change_T(s.g,state.T)
			end	
			@info("simulation # $j (alpha=$alpha max=$max_value)")
			t = alpha*max_value
			state = get_state(s,P_ref,t)
			@info("----------")
			
			states[t] = state
		end
		toc()
		serialize_to_file(states, "states_$u_name-$iter-$max_value.jld")
	# simulation on a subinterval of alpha : [alpha_i , alpha_i + alpha_interval_length]	
	elseif par == 2
		tic()
		alpha_i = parse(Float64,pargs["alpha_i"])
		alpha_interval_length = parse(Float64,pargs["alpha_interval_length"])
		for j in 1:iter
			alpha = alpha_i + (j-1)*alpha_interval_length/iter
			if j > 1
				change_T(s.g,state.T)
			end
			@info("simulation # $j (alpha=$alpha max=$max_value)")
			t = alpha*max_value
			state = get_state(s,P_ref,t)
			@info("----------")
			
			states[t] = state
		end
		toc()
		serialize_to_file(states, "states_$u_name-$max_value-$alpha_i.jld")
	elseif par == 3
		# number of processes to be used
		npcs = parse(Int,pargs["nprocs"])
		addprocs(npcs)
		# include the required code
		@everywhere include("init.jl")
		@everywhere include("data.jl")
		@everywhere include("solvers.jl")
		@everywhere include("graphs.jl")
		@everywhere function get_state(s::Simulator,P_ref::Array{Float64,1},alpha::Float64)
			P = P_ref*alpha
			@debug("norm P (2/Inf): ", norm(P,2), "/", norm(P,Inf))
			change_P(s.g,P)
			return alpha,simulation(s)
		end
		tic()
		@sync results = pmap(get_state,Simulator[s for j in 1:iter],Array{Float64,1}[P_ref*max_value for j in 1:iter],Float64[((j-1)/iter) for j in 1:iter])	
		toc()
		
		for r in results
			states[r[1]] = r[2]
		end		
		serialize_to_file(states, "states_$u_name-$iter-$max_value.jld")
	end
end
