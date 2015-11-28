include("init.jl")
include("data.jl")
include("solvers.jl")
include("graphs.jl")

# load file
#nodes,edges = load_IEEE_SLFD(ARGS[1])
nodes,edges = load_ENTSOE(ARGS[1])
t_ratio = Complex{Float64}[edge.t_ratio for edge in edges]
export_csv_data(t_ratio, "ratio.csv")
# export graph to graphml
#export_graphml(ARGS[2], nodes, edges)

# initialize simulation data
n = length(nodes) 
T = zeros(Float64, n)
# n-dimensional vector of 0s
V = Float64[n.init_voltage for n in nodes]
# node ids whose bus type is 0
PQ_ids = Int64[n.id for n in filter(n -> n.bus_type == 0, nodes)]
# set PQ bus voltages to 1 pu
V[PQ_ids] = 1.
# node id whose bus type is 3
slack_id = filter(n -> n.bus_type == 3, nodes)[1].id

Y,P0,Q0 = init_NR_data(nodes,edges)

# for ENTSOE
id_c = find_connected_graph(nodes, edges)
Y = Y[id_c,id_c]
V = V[id_c]
T = T[id_c]
P0 = P0[id_c]
Q0 = Q0[id_c]
bus_type = Int64[n.bus_type for n in nodes]
PQ_ids = findin(bus_type[id_c],0)
slack_id = findin(bus_type[id_c],3)[1]
#	for i in 1:n
#		for j in i:n
#			@printf("[%d,%d]  %10.2f - %10.2f\n", i, j, real(Y[i,j]), imag(Y[i,j]))
#		end
#end
#B = readcsv("B.csv")
#G = readcsv("G.csv")
#Y = G +B*im
#println(Y)
V,T = GS_solver(V, T, Y, P0, Q0, PQ_ids, slack_id, 3)
#export_csv_data(V, "v.csv")
#T = T*180/pi
#export_csv_data(T, "t.csv")
export_csv_data(imag(Y), "B.csv")
V,T,n_iter = NR_solver(V, T, Y, P0, Q0, PQ_ids, slack_id, 1e-8, 15)
export_csv_data(PQ_ids, "PQ.csv")
export_csv_data(V, "v.csv")
T = T*180/pi
export_csv_data(T, "t.csv")
