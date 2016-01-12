using ArgParse

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
	#export_graphml("my_export.graphml", g)

	o_args = Dict{Symbol,Any}()
	o_args[:g] = g
	o_args[:bootstrap_iter] = 3
	s = Simulator(g,NR_solver,o_args,100.,1e-8,15)
	
	# launch the simulation
	simulation(s)
	state = s.states[1]
	
	export_csv_data(state.V, "v.csv")
	state.T = state.T*180/pi
	export_csv_data(state.T, "t.csv")
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

	println(state.n_iter)
	export_csv_data(state.T, "t.csv")
	export_csv_data(state.Tdot, "tdot.csv")
elseif solver == "SD"
	p_fn = pargs["p_fn"]  
	y_fn = pargs["y_fn"] 
	g = load_graph(p_fn,y_fn) 

	o_args = Dict{Symbol,Any}()
	o_args[:d] = 1e-2
	s = Simulator(g,SD_solver,o_args,1.,1e-6,round(Int64,1e5))

	# launch the simulation
	simulation(s)
	state = s.states[1]

	export_csv_data(state.T, "t.csv")
end
