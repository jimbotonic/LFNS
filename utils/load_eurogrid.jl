using Logging, DataFrames, Graphs

include("graphs.jl")
include("data.jl")

@Logging.configure(level=ERROR)

# load graph from the specified Y and P files
#
# CSV files with no-header and comma-separated are expected
# E_fn: node_id1, node_id2
function load_graph(E_fn::AbstractString)
	# load the data from Y and P file
	E_df = load_csv_data(E_fn)

	names = AbstractString[]
	
	for i in 1:size(E_df,1)
		push!(names,strip(E_df[i,1][1:12]))
		push!(names,strip(E_df[i,2][1:12]))
	end
	
	# the size of the network is assumed to be the # of rows in P0 
	names = unique(names)
	n = length(names)

	vertices = Bus[]
	name_id = Dict{AbstractString,Int64}()
	counter = 1
	for na in names
		name_id[na] = counter
		push!(vertices, Bus(counter,na))
		counter += 1
	end
	
	edges = Line[]
	for i in 1:size(E_df,1)
		sn = strip(E_df[i,1][1:12])
		tn = strip(E_df[i,2][1:12])
		source = vertices[name_id[sn]]
		target = vertices[name_id[tn]]
		push!(edges, Line(source, target, 1.im))
	end
	
	return graph(vertices, edges, is_directed=false)
end

g = load_graph(ARGS[1])

println(g)

serialize_to_file(g,"eurogrid.jld")
export_graphml(g,"eurogrid.graphml")
