include("data.jl")

# load file
#nodes,edges = load_IEEE_SLFD(ARGS[1])
nodes,edges = load_ENTSOE(ARGS[1])

# export graph to graphml
export_graphml(ARGS[2], nodes, edges)
