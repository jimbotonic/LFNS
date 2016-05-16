using Logging

include("../data.jl")
include("../graphs.jl")

@Logging.configure(level=DEBUG)

g = load_serialized(ARGS[1])

X = Float64[]
Y = Float64[]
D = Int[]

for v in vertices(g)
	push!(X,v.lng)
	push!(Y,v.lat)
	if out_degree(v,g) > 1
		push!(D,1)
	else
		push!(D,0)
	end 
end

xc,yc = mean(X),mean(Y)
println("Center: $xc:$yc")
