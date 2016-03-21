using Logging, ConfParser

include("../init.jl")
include("../graphs.jl")
include("../simulator.jl")
include("../solvers.jl")
include("../data.jl")
include("../metrics.jl")

# loading config file
conf = ConfParse("../config.ini")
parse_conf!(conf)

# simple square lattice
n = 50
m = 50

g = generate_sq_lattice(n,m)

# initialize simulation parameters
sb = parse(Float64,retrieve(conf,"solvers","base_voltage"))
max_iter = round(Int,parse(Float64,retrieve(conf,"solvers","max_iter")))
epsilon = parse(Float64,retrieve(conf,"solvers","epsilon"))

o_args = Dict{Symbol,Any}()
o_args[:h] = parse(Float64,retrieve(conf,"rk","h"))
s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)

U = init_unif_dist(n*m)
alpha = 1e-2

# rescale distribution
P_ref = init_P3(U)
P = P_ref*alpha
change_P(s.g,P)

# dictionary computed states (pos n*i+j -> state)
states = Dict{Int,State}()

# get the contour cycle of the lattice
cycles = Array{Array{Int64,1},1}()
bcycle =  get_sq_lattice_contour_cycle(n,m)
push!(cycles,bcycle)

#for j in 1:(n-1)
#	for i in 2:m
for j in 20:20
	for i in 20:20
		@debug("Simulation: ", (n*(i-1)+j))
		T = create_vortex_on_sq_lattice(n,m,i,j)
		println("start vorticity: ", vorticity(T,cycles))
		change_T(s.g,T)
		state = simulation(s)
		println("end vorticity: ", vorticity(state.T,cycles))
		states[n*(i-1)+j] = state
	end
end

serialize_to_file(states, "states.jld")
