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

#serialize_to_file(g,"../data/ukgrid/ukgrid.jld")
#export_graphml(g,"../data/ukgrid/ukgrid.graphml")

###
# making graph planar
###

vs = vertices(g)
# correct lng of vertex 108 and 58
vs[108].lng = vs[107].lng
vs[58].lng = vs[59].lng

# remove edge 7-42
redges = Array{Pair{Int,Int},1}()
push!(redges,Pair{Int,Int}(7,42))
ng = get_pruned_graph(g, redges)
println(ng)

# remove cut component vertices (110,111,112,113,114,115,116)
n2g = get_subgraph(ng,Int[110,111,112,113,114,115,116],false)
println(n2g)

# remove sinks and branches
n3g = get_subgraph(n2g,Int[18,23,30,47,53,55,84,85,88,98,99,103,105],false)
println(n3g)

serialize_to_file(g,"../data/ukgrid/ukgrid_planar_no-sink-no-cut.jld")
export_graphml(g,"../data/ukgrid/ukgrid_planar_no-sink-no-cut.graphml")
