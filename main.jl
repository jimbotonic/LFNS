using ArgParse, Logging

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
			help = "iteration number, if parallel=yes : iteration number within the subinterval of alpha"
			required = false
		"--g_fn"
			help = "file name for serialized graph"
			required = false
		"--parallel"
			default = "no"
			help = "yes if parallel computing"
			required = false
		"--alpha_i"
			help = "first alpha value when parallelizing"
			required = false
		"--step_num"
			help = "number of steps for each instance when parallelizing"
			required = false
		"--alpha_interval_length"
			help = "length of alpha interval when parallelizing"
			required = false
	end
	return parse_args(s)
end

# parse arguments
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
	@everywhere function get_state(s::Simulator,P_ref::Array{Float64,1},alpha::Float64,max_value::Float64)
		P = P_ref*alpha*max_value
		@debug("norm P (2/Inf): ", norm(P,2), "/", norm(P,Inf))
		change_P(s.g,P)
		return t,simulation(s)
	end
	
	g = load_serialized(pargs["g_fn"])
	n = length(vertices(g))
	
	# initialize simulation parameters
	sb = 1.
	max_iter = round(Int64,1e6)
	
	epsilon = 1e-3
	# NB: Eurogrid PC has 6021 nodes and a max degree of 17
	# critical value is between 800 and 850
	max_value = 600.

	ssolver = pargs["ssolver"]
	if ssolver == "SD"
		o_args = Dict{Symbol,Any}()
		o_args[:d] = 1
		s = Simulator(g,SD_solver,o_args,sb,epsilon,max_iter)
	elseif ssolver == "RK"
		o_args = Dict{Symbol,Any}()
		o_args[:h] = .12
		s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)
	elseif ssolver == "RK2"
		o_args = Dict{Symbol,Any}()
		o_args[:h] = .12
		s = Simulator(g,RK_solver2,o_args,sb,epsilon,max_iter)
	end

	# generate a uniform distribution
	#U = init_unif_dist(n)
	#export_csv_data(U, "U.csv")
	
	u_fn = pargs["u_fn"]
	U = collect(load_csv_data(u_fn)[1])
	iter = parse(Int,pargs["iter"])
	u_name = basename(u_fn)[1:end-4]

	P_ref = init_P3(U)
	
	# if we have an initial bootstraping vector T
	bt_fn = pargs["bt_fn"]
	if typeof(bt_fn) != Void
		bT = collect(load_csv_data(bt_fn)[1])
		change_T(s.g,bT)
	end

	# array of computed states 
	states = Dict{Int,Array{Float64,1}}()

	# parallel mode	
	par = parse(Int,pargs["parallel"])
	
	# no parallelization (simulation on the whole interval of alpha values)
	if par == 1
		tic()
		for j in 1:iter
			# start at 0 (ie, "flat start")
			t = (j-1)/iter	
			if j > 1
				change_T(s.g,state.T)
			end	
			@info("simulation # $j (alpha=$t max=$max_value)")
			i,state = get_state(s,P_ref,t,max_value)
			@info("----------")
			
			states[i] = state.T
		end
		toc()
		serialize_to_file(states, "states_$u_name-$iter-$max_value.jld")
	# simulation on a subinterval of alpha : [alpha_i , alpha_i + alpha_interval_length]	
	elseif par == 2
		tic()
		# step_num = parse(Int,pargs["step_num"])
		alpha_i = parse(Float64,pargs["alpha_i"])
		alpha_interval_length = parse(Float64,pargs["alpha_interval_length"])
		for j in 1:iter
			t = alpha_i + (j-1)*alpha_interval_length/iter
			if j > 1
				change_T(s.g,state.T)
			end
			@info("simulation # $j (alpha=$t max=$max_value)")
			i,state = get_state(s,P_ref,t,max_value)
			@info("----------")
			
			states[i] = state.T
		end
		toc()
		serialize_to_file(states, "states_$u_name-$max_value-$alpha_i.jld")
	elseif par == 3
		addprocs(4)
		tic()
		results = pmap(get_state,Simulator[s for j in 1:iter],Array{Float64,1}[P_ref for j in 1:iter],Float64[(j-1)/iter for j in 1:iter],Float64[max_value for j in 1:iter])	
		toc()
	end
end
