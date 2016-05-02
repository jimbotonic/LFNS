using Logging

@Logging.configure(level=DEBUG)

include("../data.jl")
include("../plotly.jl")

# stats associated to a given random P
type Stats
	# vortex stable & convergence
	Vst::Dict{Float64,Float64}
	# vortex inside & convergence
	Vin::Dict{Float64,Float64}
	# vortex outside & convergence
	Vout::Dict{Float64,Float64}
	# no convergence
	Div::Dict{Float64,Float64}
end

# load graph
stats = load_serialized(ARGS[1])

# sort dict keys
skeys = sort(collect(keys(stats.Vst)))

X = Float64[]
Y = Float64[]

# plot vorticity
for k in skeys
	push!(X,k)
	push!(Y,stats.Vst[k])
end

# plot data
plot_scatter_data(X,Y,"scatter","markers", "Vst_dist", None)

###

# sort dict keys
skeys = sort(collect(keys(stats.Vin)))

X = Float64[]
Y = Float64[]

# plot vorticity
for k in skeys
	push!(X,k)
	push!(Y,stats.Vin[k])
end

# plot data
plot_scatter_data(X,Y,"scatter","markers", "Vin_dist", None)

###

# sort dict keys
skeys = sort(collect(keys(stats.Vout)))

X = Float64[]
Y = Float64[]

# plot vorticity
for k in skeys
	push!(X,k)
	push!(Y,stats.Vout[k])
end

# plot data
plot_scatter_data(X,Y,"scatter","markers", "Vout_dist", None)

###

# sort dict keys
skeys = sort(collect(keys(stats.Div)))

X = Float64[]
Y = Float64[]

# plot vorticity
for k in skeys
	push!(X,k)
	push!(Y,stats.Div[k])
end

# plot data
plot_scatter_data(X,Y,"scatter","markers", "Div_dist", None)

