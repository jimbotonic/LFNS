include("../graphs.jl")
include("../data.jl")

g = load_serialized(ARGS[1])
n = length(vertices(g))

pr = ones(Float64,n)
pr = pr/n
# PR reference value
# pr = load_serialized(ARGS[2])

pr = my_pagerank(g,pr,0.85,1e-11)

println(pr[1:10])
println(sum(pr))
println(maximum(pr))
println(minimum(pr))

serialize_to_file(pr,"pr.jld")
export_csv_data(pr,"pr.csv")
