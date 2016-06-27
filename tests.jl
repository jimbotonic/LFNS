include("init.jl")
include("data.jl")
include("solvers.jl")
include("graphs.jl")

using Base.Test, Logging, Distances

@Logging.configure(level=INFO)

BASE_FOLDER = "./data/tests"

###
# test RK solver
###

@info("######## Testing RK solver")

p_fn = BASE_FOLDER * "/RK/UK1/P_in.csv"
y_fn = BASE_FOLDER * "/RK/UK1/Y_in.csv"
t_fn = BASE_FOLDER * "/RK/UK1/T_out.csv"

g = load_graph(p_fn,y_fn) 
T_out = collect(load_csv_data(t_fn)[1])

o_args = Dict{Symbol,Any}()
o_args[:h] = 3e-2
s = Simulator(g,RK_solver1,o_args,1.,1e-10,round(Int64,1e5))

# launch the simulation
tic()
state = simulation(s)
toc()

#@info("T_sim: ", state.T[1:20])
#@info("T_ref: ", T_out[1:20])

d = chebyshev(state.T, T_out)
@info("distance $d")
@info("# iter: ", state.n_iter)

@test_approx_eq_eps d 0. 1e-4

###
# test SD solver
###

@info("######## Testing SD solver")

p_fn = BASE_FOLDER * "/RK/UK1/P_in.csv"
y_fn = BASE_FOLDER * "/RK/UK1/Y_in.csv"
t_fn = BASE_FOLDER * "/RK/UK1/T_out.csv"

g = load_graph(p_fn,y_fn) 
T_out = collect(load_csv_data(t_fn)[1])
T_out = uniform_phase_shift(T_out)

o_args = Dict{Symbol,Any}()
o_args[:d] = 1
s = Simulator(g,SD_solver,o_args,1.,1e-8,round(Int64,1e5))

# launch the simulation
tic()
state = simulation(s)
toc()

state.T = uniform_phase_shift(state.T)

#@info("T_sim: ", state.T[1:20])
#@info("T_ref: ", T_out[1:20])

d = chebyshev(state.T, T_out)
@info("distance $d")
@info("# iter: ", state.n_iter)

@test_approx_eq_eps d 0. 1e-4

###
# test NR solver
###

@info("######## Testing NR solver")

sys_fn = BASE_FOLDER * "/NR/IEEE/ieee14cdf.txt"
T_fn = BASE_FOLDER * "/NR/IEEE/T_out.csv"
V_fn = BASE_FOLDER * "/NR/IEEE/V_out.csv"

T_ref = collect(load_csv_data(T_fn)[1])
V_ref = collect(load_csv_data(V_fn)[1])

g = load_IEEE_SLFD(sys_fn)

o_args = Dict{Symbol,Any}()
o_args[:g] = g
o_args[:bootstrap_iter] = 0
s = Simulator(g,NR_solver,o_args,100.,1e-8,15)

tic()
state = simulation(s)
toc()

error_V = chebyshev(state.V,V_ref)
state.T = state.T*180/pi
error_T = chebyshev(state.T,T_ref)

@info("Max error in V: $error_V")
@info("Max error in theta: $error_T")

@test_approx_eq_eps error_V 0. 1e-4
@test_approx_eq_eps error_T 0. 1e-4

###
# test lattice initialization
###

@info("######## Testing lattice initialization")

# simple square lattice
n = 15
m = 25

g = generate_sq_lattice(n,m)
nv = length(vertices(g))
ne = length(edges(g))

# n*m
@test nv == n*m
@test ne == (n-1)*m + (m-1)*n

# square lattice on the sphere
g = generate_sq_lattice_on_sphere(n)

for v in vertices(g)
	@test out_degree(v,g) == 4 
end

###
# test graph analysis functions on real graphs
###

@info("######## Testing contour cycle detection")
g = load_serialized("./data/ukgrid/ukgrid.jld")

# contour cycle starting from lowest node
bcycle = [1,5,4,6,7,38,39,40,41,37,35,117,49,50,51,52,54,68,69,77,78,80,81,89,95,97,96,100,101,102,104,106,108,109,107,94,92,91,90,87,86,83,82,75,74,60,59,31,29,28,22,20,16,15,14,13,12,10,119,11,3,2]

# correct lng of vertex 108 and 58
vs = vertices(g)
vs[108].lng = vs[107].lng
vs[58].lng = vs[59].lng

# remove vertices (cut component 109 + branch 84-83) 
ignore_vids = Set{Int}([84,110,111,112,113,114,115,116])

# get contour cycle
bcycle2 = get_geolocalized_graph_contour_cycle(g,ignore_vids)

@test bcycle == bcycle2 

###
# test vorticity functions
###

a = [pi/2,pi,3pi/2,0]
v1 = vorticity(a,[1,2,3,4])
v2 = vorticity2(a,[1,2,3,4])

@test v1 == -1
@test v2 == -1
