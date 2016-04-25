addprocs(40)

@everywhere include("../init.jl")
@everywhere include("../graphs.jl")
@everywhere include("../simulator.jl")
@everywhere include("../solvers.jl")
@everywhere include("../data.jl")
@everywhere include("../metrics.jl")

using Logging, ConfParser

@Logging.configure(level=DEBUG)

# loading config file
conf = ConfParse("../config.ini")
parse_conf!(conf)

# simple square lattice
n = 49
m = 49

# generate square lattice
@everywhere g = generate_sq_lattice(n,m)

# initialize simulation parameters
@everywhere sb = parse(Float64,retrieve(conf,"solvers","base_voltage"))
@everywhere max_iter = round(Int,parse(Float64,retrieve(conf,"solvers","max_iter")))
@everywhere epsilon = parse(Float64,retrieve(conf,"solvers","epsilon"))

@everywhere o_args = Dict{Symbol,Any}()
o_args[:h] = parse(Float64,retrieve(conf,"rk","h"))
@everywhere s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)

# get the contour cycle of the lattice
@everywhere bcycle =  get_sq_lattice_contour_cycle(n,m)

# central face cycle
@everywhere ccycle = Int[]
ipos = ceil(Int,m)
jpos = ceil(Int,n)

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

T = create_vortex_on_sq_lattice2(n,m,i,j)
#T = create_antivortex_on_sq_lattice2(n,m,i,j)
change_T(s.g,T)

# generat3dde a random injection vector
U = init_rand_dist(n*m)
@everywhere P_ref = init_P3(U)

low = 0
step = 1e-2
high = 2

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
	v == 0 && has_moved = true
	return true
end

@everywhere function get_simulation_state(s::Simulator,P_ref::Array{Float64,1},alpha::Float64)
	P = P_ref*alpha
	change_P(s.g,P)
	state = simulation(s,callback_func)
	return P,state,has_moved
end

# stats associated to a given random P
type Stats
	# vortex stable & convergence
	Vst::Float64
	# vortex inside & convergence
	Vin::Float64
	# vortex outside & convergence
	Vout::Float64
	# no convergence
	Div::Float64
end

PStats = Dict{Array{Float64,1},Stats}()
alphas = collect(low:step:high)
stats = Stats(0.,0.,0.,0.)

# for the different initial P distributions
for i in 1:np
	# compute the final state for different values of alpha in parallel
	@sync results = pmap(get_simulation_state,Array{Simulator,1}[s for alpha in alphas],Array{Array{Float64,1},1}[P_ref for alpha in alphas],alphas)	
	nr = length(r)
	for r in results
		P = r[1]
		state = r[2]
		has_moved = r[3]
		if state.n_iter < max_iter 
			if !has_moved
				stats.Vst += 1.
			else
	 			v = vorticity(state.T,bcycle)
				if v == 1
					stats.Vin += 1.
				else
					stats.Vout += 1.
				end
			end
		else
			stats.Div += 1.
		end 
	end		
	stats.Vst /= nr
	stats.Vin /= nr
	stats.Vout /= nr
	stats.Div /= nr

	PStats[P] = stats
end


