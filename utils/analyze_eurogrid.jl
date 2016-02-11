include("../data.jl")
include("../solvers.jl")
include("../graphs.jl")

# save generated graphs?
SAVE = true

# load graph
g = load_serialized("../data/eurogrid/eurogrid.jld")

println("-> INITIAL GRAPH")
println("# vertices: ", length(vertices(g)))
println("# edges: ", length(edges(g)))
println("---")

cs = connected_components(g)

println("-> COMPONENTS")
println("# components: ", length(cs))
println("components size: ", Int64[length(c) for c in cs])
println("---")

mc = get_principal_component(g)

gp = get_subgraph(g, mc)
SAVE && serialize_to_file(gp, "../data/eurogrid/eurogrid_pc.jld")
SAVE && export_graphml(gp, "../data/eurogrid/eurogrid_pc.graphml")
SAVE && export_topology_graphml(gp, "../data/eurogrid/eurogrid_pc_topology.graphml")

avg_deg,min_deg,max_deg = get_avg_min_max_degree(gp)

println("-> PRINCIPAL COMPONENT")
println("# vertices: ", length(vertices(gp)))
println("# edges: ", length(edges(gp)))
println("avg/min/max degree: $avg_deg/$min_deg/$max_deg")
println("---")

n = length(vertices(gp))
# get adjacency matrix
A = get_adjacency_matrix(gp)
# set null non-diagonal entries to Inf
for i in 1:n
	for j in 1:n
		if i != j
			if A[i,j] == 0
				A[i,j] = 1e2
			end
		end
	end
end

dists = floyd_warshall(A)
serialize_to_file(dists, "../data/eurogrid/eurogrid_pc_fw.jld")

ew = ones(length(edges(gp)))
ke, kw = prim_minimum_spantree(gp, ew, vertices(gp)[1])
#println("# tree edges: ", length(ke))

reids = setdiff(Int64[edge.id for edge in edges(gp)],Int64[edge.id for edge in ke])
#println("# out-tree edges: ", length(reids))

gpt = get_pruned_graph(gp, reids)
SAVE && serialize_to_file(gpt, "../data/eurogrid/eurogrid_pc_st.jld")

println("---")
println("-> SPANNING TREE OF PRINCIPAL COMPONENT")
println("# vertices: ", length(vertices(gpt)))
println("# edges: ", length(edges(gpt)))
println("---")
 
cycles = get_cycle_base(gp)
SAVE && serialize_to_file(cycles, "../data/eurogrid/eurogrid_pc_cb.jld")

println("-> FUNDAMENTAL CYCLE BASIS")
println("# basis cycles: ", length(cycles))
