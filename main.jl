include("data.jl")
include("solver.jl")

# load file
nodes,edges = load_IEEE_SLFD(ARGS[1])
#nodes,edges = load_ENTSOE(ARGS[1])

# export graph to graphml
# export_graphml(ARGS[2], nodes, edges)

# initialize simulation data
n = length(nodes) 
T = zeros(Float64, n)
# n-dimensional vector of 0s
V = Float64[n.init_voltage for n in nodes]
# node ids whose bus type is 0
PQ_ids = Int64[n.id for n in filter(n -> n.bus_type == 0, nodes)]
# node id whose bus type is 3
slack_id = filter(n -> n.bus_type == 3, nodes)[1].id

Y,P0,Q0 = init_simulation_data(nodes,edges)

	for i in 1:n
		for j in i:n
			@printf("[%d,%d]  %10.2f - %10.2f\n", i, j, real(Y[i,j]), imag(Y[i,j]))
		end
end
#NR_solver(Y, V, T, P0, Q0, PQ_ids, slack_id, epsilon::Float=1e-4, iter_max::Int=50)
