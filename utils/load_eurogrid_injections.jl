using Logging, DataFrames, Graphs

include("graphs.jl")
include("data.jl")

@Logging.configure(level=DEBUG)
       
V_df = load_csv_data(ARGS[1])
g = load_serialized(ARGS[2])
n = length(vertices(g))

P = Dict{AbstractString, Float64}()

for i in 1:size(V_df,1)
	name = strip(V_df[i,1][1:12])
	power = V_df[i,2]
	if haskey(P,name)
		P[name] += power
	else
		P[name] = power
	end
end

vnames = Dict{AbstractString,Int}()
for v in vertices(g)
	vnames[v.name] = v.id
end

A = zeros(Float64,length(vertices(g)))

for n in keys(P)
	if haskey(vnames,n)
		A[vnames[n]] = P[n] 
	end
end

## LEIP
#println(A[931])
## PENLY
#println(A[2944])
#println(A[2945])
#println(A[2946])
## CHAM
#println(A[1196])
#println(A[1197])


println(A[1:20])
println("Sum: ", sum(A))
println("Min: ", minimum(A))
println("Max: ", maximum(A))

println("Max consumer: ", vertices(g)[indmax(A)].name)
println("Max producer: ", vertices(g)[indmin(A)].name)
	
export_csv_data(A, "Euro_P.csv")
