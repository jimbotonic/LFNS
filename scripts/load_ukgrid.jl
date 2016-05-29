using Logging, DataFrames, Graphs

include("../graphs.jl")
include("../data.jl")

@Logging.configure(level=DEBUG)

# load graph
#
# CSV files with no-header and space-separated are expected
# E_fn: node_id1, node_id2, B_value
# pos_fn: lat,lng
function load_graph(E_fn::AbstractString, pos_fn::AbstractString)
	# load the data
	E_df = load_csv_data(E_fn)
	pos_df = load_csv_data(pos_fn)

	# adding vertices
	vertices = Bus[]
	counter = 1
	for i in 1:size(pos_df,1)
		bus = Bus(counter,string(counter))
		bus.lng = pos_df[i,1]
		bus.lat = pos_df[i,2]
		push!(vertices, bus)
		counter += 1
	end
	
	# adding edges
	s = Set{UInt64}()
	ecounter = 1
	edges = Line[]
	for i in 1:size(E_df,1)
		sn = string(E_df[i,1]+1)
		tn = string(E_df[i,2]+1)
		si = parse(Int,sn)
		ti = parse(Int,tn)
		if si != ti
			if si < ti
				i1 = si
				i2 = ti
				n1 = sn
				n2 = tn
			else
				i1 = ti
				i2 = si
				n1 = tn
				n2 = sn
			end
			h = hash(n1*n2)
			if !(h in s)
				source = vertices[i1]
				target = vertices[i2]
				@debug("adding edge ($i1,$i2)")
				push!(edges, Line(ecounter, source, target, -1.im))
				push!(s,h)
				ecounter += 1
			end
		end
	end
	
	return graph(vertices, edges, is_directed=false)
end

g = load_graph(ARGS[1],ARGS[2])
println(g)

serialize_to_file(g,"../data/ukgrid/ukgrid.jld")
export_graphml(g,"../data/ukgrid/ukgrid.graphml")
