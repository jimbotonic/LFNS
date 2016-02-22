using Logging
#using ConfParser

@Logging.configure(level=DEBUG)

include("../data.jl")
include("../solvers.jl")
include("../simulator.jl")
include("../metrics.jl")
include("../plotly.jl")

# loading config file
#conf = ConfParse("../config.ini")

# load graph
g = load_serialized(ARGS[1])
# load simulation states
states = load_serialized(ARGS[2])

# initialize simulation parameters
sb = 1.
max_iter = round(Int64,1e5)
epsilon = 1e-8
o_args = Dict{Symbol,Any}()
o_args[:h] = 3e-2
s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)
sp = get_sparams(s)

X = Float64[]
Y = Float64[]

counter = 1
for T in states
	l2 = get_lambda2(T, sp.Y)
	push!(X,counter)
	push!(Y,l2)
	counter += 1
	@debug("l2: $l2 (counter=$counter)")
end

# plot data
plot_scatter_data(X,Y,"scatter","lines+markers", "eurogrid_lambda2", None)
