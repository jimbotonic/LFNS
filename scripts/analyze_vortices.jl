using Distances

include("../init.jl")
include("../data.jl")
include("../simulator.jl")
include("../metrics.jl")

n = 50
m = 50

# position -> state
states = load_serialized("states.jld")

cycles = Array{Array{Int64,1},1}()
bcycle =  get_sq_lattice_contour_cycle(n,m)
push!(cycles,bcycle)

for j in 1:10
	for i in 2:11
		k = (i-1)*n+j
		T = states[k].T
		@info("Vortex at position ($i,$j)")
		@info("norm: ", norm(T))
		@info("vorticity: ", vorticity(T,cycles))
		@info("---")
	end
end
