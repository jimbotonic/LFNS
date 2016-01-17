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

direct_cycles(cb)
