using Logging

@Logging.configure(level=DEBUG)

#addprocs(4)

@everywhere include("../data.jl")
@everywhere include("../solvers.jl")
@everywhere include("../simulator.jl")
@everywhere include("../metrics.jl")
@everywhere include("../plotly.jl")

# load graph
g = load_serialized(ARGS[1])
# load simulation states
states = load_serialized(ARGS[2])
# load cycle basis
cycles = load_serialized(ARGS[3])

# initialize simulation parameters
sb = 1.
max_iter = round(Int64,1e5)
epsilon = 1e-6
o_args = Dict{Symbol,Any}()
o_args[:h] = .12
s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)
sp = get_sparams(s)

max_v = 0.
mean_v = 0.

# compute vorticity vectors
for k in keys(states)
	T = states[k].T
	V = vorticity(T, cycles)
	mean_v += mean(V)
	max = maximum(V)
	if max > max_v
		max_v = max
	end
end

@info("Max vorticity: ", max_v)
@info("Mean vorticity: ", (mean_v/length(states)))

X = Float64[]
Y = Float64[]

# plot lambda 2
for k in keys(states)
	T = states[k].T
	V = vorticity(T, cycles)
	push!(X,k)
	push!(Y,mean(V))
	@debug("mean V: ", mean(V))
end

# plot data
plot_scatter_data(X,Y,"scatter","markers", "eurogrid_vorticity", None)
quit()

X = Float64[]
Y = Float64[]

# plot lambda 2
for k in keys(states)
	T = states[k].T
	l2 = get_lambda2(T, sp.Y)
	push!(X,k)
	push!(Y,l2)
	@debug("l2: $l2")
end

# plot data
plot_scatter_data(X,Y,"scatter","lines+markers", "eurogrid_lambda2", None)
