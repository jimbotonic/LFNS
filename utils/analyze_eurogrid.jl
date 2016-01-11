include("../data.jl")
include("../solvers.jl")
include("../graphs.jl")


# load graph
g = load_serialized("../data/eurogrid.jld")

println("# vertices: ", length(vertices(g)))
println("# edges: ", length(edges(g)))
println("---")

cs = connected_components(g)
println("# components: ", length(cs))
println("components size: ", [length(c) for c in cs])
println("---")

# display small components node names
#for i in 2:length(cs)
#	println([b.name for b in cs[i]])
#	println("---")
#end

mc = get_principal_component(g)
println("# vertices of principal component: ", length(mc))
gp = get_subgraph(g, mc)
serialize_to_file(gp, "../data/eurogrid_pc.jld")

println("# vertices (pc): ", length(vertices(gp)))
println("# edges (pc): ", length(edges(gp)))
println("---")

ew = ones(length(edges(gp)))
ke, kw = prim_minimum_spantree(gp, ew, vertices(gp)[1])

println("# tree edges: ", length(ke))

reids = setdiff(Int64[edge.id for edge in edges(gp)],Int64[edge.id for edge in ke])
println("# out-tree edges: ", length(reids))

gpt = get_pruned_graph(gp, reids)
serialize_to_file(gpt, "../data/eurogrid_pc_st.jld")

n = length(vertices(gpt))
m = length(edges(gpt))
println("# vertices (pc st): $n")
println("# edges (pc st): $m")
println("---")
 
cycles = get_cycle_base(gp)
println("# basis cycles: ", length(cycles))
println("basis cycle 1: ", cycles[1])
serialize_to_file(cycles, "../data/eurogrid_pc_cb.jld")
