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


###
# create 2 vortices
###

@everywhere begin 
	# double vortex on the same row
	i1 = 25; j1 = 20
	i2 = 25; j2 = 31
	# double vortex on the same line
	#i1 = 20; j1 = 25
	#i2 = 31; j2 = 25
	# double vortex in the diagonal
	#i1 = 20; j1 = 20
	#i2 = 31; j2 = 31
end

T1 = create_vortex_on_sq_lattice(n,m,i1,j1)
T2 = create_vortex_on_sq_lattice(n,m,i2,j2)
#T2 = create_antivortex_on_sq_lattice(n,m,i2,j2)
T = T1+T2

# central face cycle
@everywhere begin 
	ccycle1 = Int[]
	ccycle2 = Int[]

	push!(ccycle1, (i1-1)*n+j1)
	push!(ccycle1, (i1-1)*n+j1+1)
	push!(ccycle1, i1*n+j1+1)
	push!(ccycle1, i1*n+j1)
	
	push!(ccycle2, (i2-1)*n+j2)
	push!(ccycle2, (i2-1)*n+j2+1)
	push!(ccycle2, i2*n+j2+1)
	push!(ccycle2, i2*n+j2)
end

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

@everywhere has_moved1 = false
@everywhere has_moved2 = false

# solver callback function
@everywhere function callback_func(sp::SParams,n_iter::Int,error::Float64)
	global has_moved1,has_moved2
	if !has_moved1
		# the vortex 1 has moved?
		if vorticity(sp.T,ccycle1) == 0
			has_moved1 = true
		end
	end
	if !has_moved2
		# the vortex 2 has moved?
		if vorticity(sp.T,ccycle2) == 0
			has_moved2 = true
		end
	end
	return true
end

@everywhere function get_simulation_state(s::Simulator,P_ref::Array{Float64,1},alpha::Float64)
	global has_moved1 = false
	global has_moved2 = false
	P = P_ref*alpha
	change_P(s.g,P)
	state = simulation(s,callback_func)
	return alpha,state,has_moved1,has_moved2
end

# stats associated to a given random P
type AStats2
	# vortex stable & convergence
	Vst1::Dict{Float64,Float64}
	Vst2::Dict{Float64,Float64}
	# vortices inside & convergence
	Vin::Dict{Float64,Float64}
	# 1 and 2 vortices outside & convergence
	Vout1::Dict{Float64,Float64}
	Vout2::Dict{Float64,Float64}
	# no convergence
	Div::Dict{Float64,Float64}
end

# initialize AStats
astats2 = AStats2(Dict{Float64,Float64}(),Dict{Float64,Float64}(),Dict{Float64,Float64}(),Dict{Float64,Float64}(),Dict{Float64,Float64}(),Dict{Float64,Float64}())

for alpha in alphas
	astats2.Vst1[alpha] = 0.
	astats2.Vst2[alpha] = 0.
	astats2.Vin[alpha] = 0.
	astats2.Vout1[alpha] = 0.
	astats2.Vout2[alpha] = 0.
	astats2.Div[alpha] = 0.
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
		has_moved1 = r[3]
		has_moved2 = r[4]
		if state.n_iter < max_iter 
			if !has_moved1 && !has_moved2
				astats2.Vst1[alpha] += 1.
				as[alpha] = :Vst1
			elseif (!has_moved1 && has_moved2) || (has_moved1 && !has_moved2)
				astats2.Vst2[alpha] += 1.
				as[alpha] = :Vst2
			else
	 			v = vorticity(state.T,bcycle)
				if v == 2
					astats2.Vin[alpha] += 1.
					as[alpha] = :Vin
				elseif v == 1
					astats2.Vout1[alpha] += 1.
					as[alpha] = :Vout1
				elseif v == 0
					astats2.Vout2[alpha] += 1.
					as[alpha] = :Vout2
				end
			end
		else
			astats2.Div[alpha] += 1.
			as[alpha] = :Div
		end 
	end		
	
	pstats[P_ref] = as
end

# normalize
for alpha in alphas
	astats2.Vst1[alpha] /= np
	astats2.Vst2[alpha] /= np
	astats2.Vin[alpha] /= np
	astats2.Vout1[alpha] /= np
	astats2.Vout2[alpha] /= np
	astats2.Div[alpha] /= np
end

println(astats2)
serialize_to_file(astats2, "AStats2_$n-$m-$np-$low-$step-$high.jld")
serialize_to_file(pstats, "PStats2_$n-$m-$np-$low-$step-$high.jld")

