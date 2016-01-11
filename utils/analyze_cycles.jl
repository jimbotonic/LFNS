include("../data.jl")
include("../solvers.jl")
include("../graphs.jl")

# load CB
cb = load_serialized("../data/eurogrid_pc_cb.jld")

println("# cycles: ", length(cb))

# display 10 first cycles
for i in 1:10
	println("cycle $i: ", cb[i])
end
println("cycle 13: ", cb[13])

println("cycle lengths: ", [length(c) for c in cb])
println("Max cycle length: ", maximum([length(c) for c in cb]))
println("Min cycle length: ", minimum([length(c) for c in cb]))

cl1 = filter(x -> length(x) == 1, cb)
println("cycle length 1: ", cl1)
 

c1 = filter(x -> 1 in x, cb)
lc = filter(x -> sum(x) == minimum(sum(c1)), c1)
println(lc)


