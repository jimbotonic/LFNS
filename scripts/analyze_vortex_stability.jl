using Logging, ConfParser

include("../init.jl")
include("../graphs.jl")
include("../simulator.jl")
include("../solvers.jl")
include("../data.jl")
include("../metrics.jl")

@Logging.configure(level=DEBUG)

# loading config file
conf = ConfParse("../config.ini")
parse_conf!(conf)

addprocs(45)

@everywhere begin
	include("../init.jl")
	include("../graphs.jl")
	include("../simulator.jl")
	include("../solvers.jl")
	include("../data.jl")
	include("../metrics.jl")

	# simple square lattice
	n = 49
	m = 49
end

# generate square lattice
g = generate_sq_lattice(n,m)

# initialize simulation parameters
sb = parse(Float64,retrieve(conf,"solvers","base_voltage"))
max_iter = round(Int,parse(Float64,retrieve(conf,"solvers","max_iter")))
epsilon = parse(Float64,retrieve(conf,"solvers","epsilon"))

o_args = Dict{Symbol,Any}()
o_args[:h] = parse(Float64,retrieve(conf,"rk","h"))
s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)

# get the contour cycle of the lattice
bcycle =  get_sq_lattice_contour_cycle(n,m)

# central face cycle
@everywhere begin 
	ccycle = Int[]
	ipos = ceil(Int,m/2)
	jpos = ceil(Int,n/2)

	push!(ccycle, (ipos-1)*n+jpos)
	push!(ccycle, (ipos-1)*n+jpos+1)
	push!(ccycle, ipos*n+jpos+1)
	push!(ccycle, ipos*n+jpos)
end

###
# create vortex
###

# vortex or antivortex in the middle
# choose middle square
i = 25
j = 25
T = create_vortex_on_sq_lattice(n,m,i,j)
#T = create_antivortex_on_sq_lattice(n,m,i,j)

# set new angles
change_T(s.g,T)

low = 0
step = 1e-2
high = 1.25

# number of random initial random distribution
np = 100

# set of alpha values
alphas = collect(low:step:high)

###
# parallelization
###

@everywhere has_moved = false

# solver callback function
@everywhere function callback_func(sp::SParams,n_iter::Int,error::Float64)
	global has_moved
	if !has_moved
		# the vortex has moved?
		if vorticity(sp.T,ccycle) == 0
			has_moved = true
		end
	end
	return true
end

@everywhere function get_simulation_state(s::Simulator,P_ref::Array{Float64,1},alpha::Float64)
	global has_moved = false
	P = P_ref*alpha
	change_P(s.g,P)
	state = simulation(s,callback_func)
	return alpha,state,has_moved
end

# stats associated to a given random P
type AStats
	# vortex stable & convergence
	Vst::Dict{Float64,Float64}
	# vortex inside & convergence
	Vin::Dict{Float64,Float64}
	# vortex outside & convergence
	Vout::Dict{Float64,Float64}
	# no convergence
	Div::Dict{Float64,Float64}
end

# initialize AStats
stats = AStats(Dict{Float64,Float64}(),Dict{Float64,Float64}(),Dict{Float64,Float64}(),Dict{Float64,Float64}())
for alpha in alphas
	stats.Vst[alpha] = 0.
	stats.Vin[alpha] = 0.
	stats.Vout[alpha] = 0.
	stats.Div[alpha] = 0.
end

# P -> (alpha -> symbol)
pstats = Dict{Array{Float64,1},Dict{Float64,Symbol}}()

# for the different initial P distributions
for i in 1:np
	# generate a random injection vector
	U = init_rand_dist(n*m)
	P_ref = init_P3(U)
	as = Dict{Float64,Symbol}()

	# compute the final state for different values of alpha in parallel
	@sync results = pmap(get_simulation_state,Simulator[s for alpha in alphas],Array{Float64,1}[P_ref for alpha in alphas],alphas)	
	nr = length(results)
	for r in results
		alpha = r[1]
		state = r[2]
		has_moved = r[3]
		if state.n_iter < max_iter 
			if !has_moved
				stats.Vst[alpha] += 1.
				as[alpha] = :Vst
			else
	 			v = vorticity(state.T,bcycle)
				if v == 1
					stats.Vin[alpha] += 1.
					as[alpha] = :Vin
				else
					stats.Vout[alpha] += 1.
					as[alpha] = :Vout
				end
			end
		else
			stats.Div[alpha] += 1.
			as[alpha] = :Div
		end 
	end		
	
	pstats[P_ref] = as
end

# normalize
for alpha in alphas
	stats.Vst[alpha] /= np
	stats.Vin[alpha] /= np
	stats.Vout[alpha] /= np
	stats.Div[alpha] /= np
end

println(stats)
serialize_to_file(stats, "AStats_$n-$m-$np-$low-$step-$high.jld")
serialize_to_file(pstats, "PStats_$n-$m-$np-$low-$step-$high.jld")

