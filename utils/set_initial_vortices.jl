include("../init.jl")
include("../graphs.jl")
include("../simulator.jl")
include("../solvers.jl")
include("../data.jl")
include("../metrics.jl")

using Logging, ConfParser

@Logging.configure(level=DEBUG)

# loading config file
conf = ConfParse("../config.ini")
parse_conf!(conf)

# simple square lattice
n = 49
m = 49

g = generate_sq_lattice(n,m)

# initialize simulation parameters
sb = parse(Float64,retrieve(conf,"solvers","base_voltage"))
max_iter = round(Int,parse(Float64,retrieve(conf,"solvers","max_iter")))
epsilon = parse(Float64,retrieve(conf,"solvers","epsilon"))

o_args = Dict{Symbol,Any}()
o_args[:h] = parse(Float64,retrieve(conf,"rk","h"))
s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)

# initialize the injections with a uniform distribution
#U = init_unif_dist(n*m)
U = init_rand_dist(n*m)
alpha = 1e-1

# rescale distribution and set injections
P_ref = init_P3(U)
P = P_ref*alpha
change_P(s.g,P)

# dictionary computed states (pos n*i+j -> state)
states = Dict{Int,State}()

# get the contour cycle of the lattice
#cycles = Array{Array{Int64,1},1}()
bcycle =  get_sq_lattice_contour_cycle(n,m)
#push!(cycles,bcycle)

# solver callback function
function callback_func(sp::SParams,n_iter::Int,error::Float64)
	@info("vorticity (iteration: $n_iter, error: $error): ", vorticity(sp.T,bcycle))
	X,Y,Z = find_vortices_in_sq_lattice(n,m,sp.T)
	if length(X) > 0
		println(X)
		println(Y)
		println(Z)
		@info("position: ", X[1], Y[1])
	end
end

# choose middle square
i = 25
j = 25
T = create_vortex_on_sq_lattice2(n,m,i,j)

@debug("angle ($i,$j): ", T[(i-1)*n+j])
@debug("angle ($i,$(j+1)): ", T[(i-1)*n+j+1])
@debug("angle ($(i+1),$j): ", T[i*n+j])
@debug("angle ($(i+1),$(j+1)): ", T[i*n+j+1])

@info("start vorticity: ", vorticity(T,bcycle))
change_T(s.g,T)
state = simulation(s,callback_func)
@info("end vorticity: ", vorticity(state.T,bcycle))

