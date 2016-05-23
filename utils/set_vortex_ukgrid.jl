using Logging

include("../init.jl")
include("../data.jl")
include("../graphs.jl")

@Logging.configure(level=DEBUG)

g = load_serialized(ARGS[1])

X = Float64[]
Y = Float64[]
X2 = Float64[]
Y2 = Float64[]
nvids = Set{Int}()
on_vid = Dict{Int,Int}()
oid_v = Dict{Int,Bus}()
T = Float64[]

counter = 1
for v in vertices(g)
	push!(X,v.lat)
	push!(Y,v.lng)
	# ignore sinks
	if out_degree(v,g) > 1
		push!(nvids,v.id)
		on_vid[v.id] = counter
		oid_v[v.id] = v
		push!(X2,v.lat)
		push!(Y2,v.lng)
		counter += 1
	end
end

#xc,yc = mean(X2),mean(Y2)
xc,yc = mean(X),mean(Y)
println("Center: $xc:$yc")

# get contour cycle
bcycle = get_geolocalized_graph_contour_cycle(g)
@debug(bcycle)

# 
for i in 1:length(bcycle)
	t = compute_edge_angle(xc, yc, X[bcycle[i]], Y[bcycle[i]])
	println(rad2deg(t))
end
