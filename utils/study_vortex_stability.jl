using Logging, ConfParser

@Logging.configure(level=DEBUG)

# loading config file
conf = ConfParse("../config.ini")
parse_conf!(conf)

#addprocs(3)

@everywhere include("../init.jl")
@everywhere include("../graphs.jl")
@everywhere include("../simulator.jl")
@everywhere include("../solvers.jl")
@everywhere include("../data.jl")
@everywhere include("../metrics.jl")

# simple square lattice
n = 49
m = 49

# generate square lattice
@everywhere g = generate_sq_lattice(n,m)

# initialize simulation parameters
sb = parse(Float64,retrieve(conf,"solvers","base_voltage"))
max_iter = round(Int,parse(Float64,retrieve(conf,"solvers","max_iter")))
epsilon = parse(Float64,retrieve(conf,"solvers","epsilon"))

o_args = Dict{Symbol,Any}()
o_args[:h] = parse(Float64,retrieve(conf,"rk","h"))
@everywhere s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)

# get the contour cycle of the lattice
@everywhere bcycle =  get_sq_lattice_contour_cycle(n,m)

# central face cycle
@everywhere ccycle = Int[]
ipos = ceil(Int,m/2)
jpos = ceil(Int,n/2)

push!(ccycle, (ipos-1)*n+jpos)
push!(ccycle, (ipos-1)*n+jpos+1)
push!(ccycle, ipos*n+jpos+1)
push!(ccycle, ipos*n+jpos)

###
# create vortex
###

# choose central  square
i = 25
j = 25

T = create_vortex_on_sq_lattice(n,m,i,j)
#T = create_antivortex_on_sq_lattice2(n,m,i,j)
change_T(s.g,T)


low = 0
step = 1e-2
high = 0.1

# number of random initial random distribution
np = 100

###
# parallelization
###

@everywhere has_moved = false

# solver callback function
@everywhere function callback_func(sp::SParams,n_iter::Int,error::Float64)
	# the vortex has moved?
	v = vorticity(sp.T,ccycle)
	if v == 0
		has_moved = true
	end
	return true
end

@everywhere function get_simulation_state(s::Simulator,P_ref::Array{Float64,1},alpha::Float64)
	P = P_ref*alpha
	change_P(s.g,P)
	state = simulation(s,callback_func)
	return state,alpha,has_moved
end

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

alphas = collect(low:step:high)

# initialize stats
stats = Stats(Dict{Float64,Float64}(),Dict{Float64,Float64}(),Dict{Float64,Float64}(),Dict{Float64,Float64}())
for alpha in alphas
	stats.Vst[alpha] = 0.
	stats.Vin[alpha] = 0.
	stats.Vout[alpha] = 0.
	stats.Div[alpha] = 0.
end

# for the different initial P distributions
for i in 1:1
	# generate a random injection vector
	U = init_rand_dist(n*m)
	P_ref = init_P3(U)

	# compute the final state for different values of alpha in parallel
	@sync results = pmap(get_simulation_state,Simulator[s for alpha in alphas],Array{Float64,1}[P_ref for alpha in alphas],alphas)	
	nr = length(results)
	for r in results
		state = r[1]
		alpha = r[2]
		has_moved = r[3]
		if state.n_iter < max_iter 
			if !has_moved
				stats.Vst[alpha] += 1.
			else
	 			v = vorticity(state.T,bcycle)
				if v == 1
					stats.Vin[alpha] += 1.
				else
					stats.Vout[alpha] += 1.
				end
			end
		else
			stats.Div[alpha] += 1.
		end 
	end		

	# normalize
	for alpha in alphas
		stats.Vst[alpha] /= np
		stats.Vin[alpha] /= np
		stats.Vout[alpha] /= np
		stats.Div[alpha] /= np
	end
end

serialize_to_file(stats, "stats_$np-$low-$step-$high.jld")

