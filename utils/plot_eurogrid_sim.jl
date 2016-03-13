using Logging

@Logging.configure(level=DEBUG)

include("../data.jl")
include("../solvers.jl")
include("../simulator.jl")
include("../metrics.jl")
include("../plotly.jl")

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

# plot vorticity
for k in keys(states)
	T = states[k].T
	V = vorticity(T, cycles)
	push!(X,k)
	push!(Y,var(V))
	#@debug("mean V: ", mean(V))
end

# plot data
plot_scatter_data(X,Y,"scatter","markers", "eurogrid_vorticity", None)

X = Float64[]
Y = Float64[]

par = true

if !par
	# plot lambda 2
	for k in keys(states)
		T = states[k].T
		l2 = get_lambda2(T, sp.Y)
		push!(X,k)
		push!(Y,l2)
		#@debug("l2: $l2")
	end
else
	addprocs(40)

	@everywhere include("../data.jl")
	@everywhere include("../solvers.jl")
	@everywhere include("../simulator.jl")
	@everywhere include("../metrics.jl")
	@everywhere include("../plotly.jl")
	@everywhere function get_par_lambda2(k::Float64,T::Array{Float64,1},Y::SparseMatrixCSC{Complex{Float64},Int64})
		return k,get_lambda2(T,Y)
	end

	kstates = keys(states)
	tic()
	@sync results = pmap(get_par_lambda2,collect(kstates),Array{Float64,1}[states[k].T for k in kstates],SparseMatrixCSC{Complex{Float64},Int64}[sp.Y for k in kstates])	
	toc()

	for r in results
		push!(X,r[1])
		push!(Y,r[2])
	end		
end

# plot data
plot_scatter_data(X,Y,"scatter","markers", "eurogrid_lambda2", None)
