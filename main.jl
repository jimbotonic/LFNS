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
			help = "file name for P0 vector"
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

	# initialize simulation data
	T = zeros(Float64, n)
	# n-dimensional vector of 0s
	V = Float64[v.init_voltage for v in vs]
	# node ids whose bus type is 0
	PQ_ids = Int64[v.id for v in filter(v -> v.bus_type == 0, vs)]
	# set PQ bus voltages to 1 pu
	V[PQ_ids] = 1.
	# node id whose bus type is 3
	slack_id = filter(v -> v.bus_type == 3, vs)[1].id
	Y,P0,Q0 = init_NR_data(g)

	# set of node ids which are part of the connected component containing the slack bus 
	sc_ids = get_slack_component_ids(g)
	Y = Y[sc_ids, sc_ids]
	V = V[sc_ids]
	T = T[sc_ids]
	P0 = P0[sc_ids]
	Q0 = Q0[sc_ids]
	# bus types of the slack component 	
	bus_type = Int64[v.bus_type for v in vs]
	# PQ bus positions in the slack component
	PQ_pos = findin(bus_type[sc_ids], 0)
	# slack position in the slack component
	slack_pos = findin(bus_type[sc_ids],3)[1]
	# bootstrap simulation
	V,T = GS_solver(V, T, Y, P0, Q0, PQ_pos, slack_pos, 3)
	# Newton-Raphson solver
	V,T,n_iter = NR_solver(V, T, Y, P0, Q0, PQ_pos, slack_pos, 1e-8, 15)
	export_csv_data(V, "v.csv")
	T = T*180/pi
	export_csv_data(T, "t.csv")
elseif solver == "RK"
	p_fn = pargs["p_fn"] # vector of initial powers 
	y_fn = pargs["y_fn"] # initial admittance matrix
	Y,P0 = load_RK_data(p_fn,y_fn) # load Admittance matrix and injected/consumed powers
	V = ones(length(P0)) # set all voltages to 1
	T = zeros(length(P0)) # flat start all angles set to zero
	h, epsilon, step_max = 0.01, 1e-11, round(Int64,1e5)

	T,Tdot,n_iter=RK_solver1(T, h, V, Y, P0, epsilon, step_max)

	#println(n_iter)
	export_csv_data(T, "t.csv")
	export_csv_data(Tdot, "tdot.csv")
elseif solver == "SD"

end
