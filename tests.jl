include("init.jl")
include("data.jl")
include("solvers.jl")
include("graphs.jl")
include("metrics.jl")

using Base.Test, Logging, Distances

@Logging.configure(level=INFO)

BASE_FOLDER = "./data/tests"

###
# test RK solver
###

@info("############ Testing RK solver")
@info("######## UK grid without dissipation (RK solver)")

p_fn = BASE_FOLDER * "/RK/UK1/P_in.csv"
y_fn = BASE_FOLDER * "/RK/UK1/Y_in.csv"
t_fn = BASE_FOLDER * "/RK/UK1/T_out.csv"

g = load_graph(p_fn,y_fn) 
T_out = uniform_phase_shift(collect(load_csv_data(t_fn)[1]))

o_args = Dict{Symbol,Any}()
o_args[:h] = 1e-2
s = Simulator(g,RK_solver1,o_args,1.,1e-10,round(Int64,1e5))

# launch the simulation
@time state = simulation(s)

state.T = uniform_phase_shift(state.T)

#@info("T_sim: ", state.T[1:20])
#@info("T_ref: ", T_out[1:20])

d = chebyshev(state.T, T_out)
@info("distance $d")
@info("# iter: ", state.n_iter)

@test_approx_eq_eps d 0. 1e-4

@info("######## UK grid with dissipation (RK solver)")

p_fn = BASE_FOLDER * "/RK/UK2/P_losses.csv"
y_fn = BASE_FOLDER * "/RK/UK2/Y_in2.csv"
t_fn = BASE_FOLDER * "/RK/UK2/theta.csv"

g = load_graph(p_fn,y_fn) 
T_out = uniform_phase_shift(collect(load_csv_data(t_fn)[1]))

o_args = Dict{Symbol,Any}()
o_args[:h] = 1e-2
s = Simulator(g,RK_solver1,o_args,1.,1e-10,round(Int64,1e5))

# launch the simulation
@time state = simulation(s)

state.T = uniform_phase_shift(state.T)

#@info("T_sim: ", state.T[1:20])
#@info("T_ref: ", T_out[1:20])

d = chebyshev(state.T, T_out)
@info("distance $d")
@info("# iter: ", state.n_iter)

@test_approx_eq_eps d 0. 1e-4


## Too long to compute...

@info("######## Eurogrid without dissipation (RK solver)")

#g = load_jld_serialized("g", "./data/eurogrid/eurogrid.jld") 
#@info("# vertices: ", length(vertices(g)))

#pc_ids = Int64[v.id for v in get_principal_component(g)]
#pc = get_subgraph(g,pc_ids)
#serialize_to_jld(pc, "g", "./data/eurogrid/eurogrid_pc")

p_fn = BASE_FOLDER * "/RK/Eurogrid/P_in.csv"
t_fn = BASE_FOLDER * "/RK/Eurogrid/T_out_1e-10.csv"

# load principal component of eurogrid
g = load_jld_serialized("g", "./data/eurogrid/eurogrid_pc.jld") 
@info("# vertices: ", length(vertices(g)))

T_out = uniform_phase_shift(collect(load_csv_data(t_fn)[1]))

P = collect(load_csv_data(p_fn)[1])
change_P(g,P)

o_args = Dict{Symbol,Any}()
o_args[:h] = 1e-2
s = Simulator(g,RK_solver1,o_args,1.,1e-12,round(Int64,1e8))

# export admittance matrix for TC code
#vs = vertices(g)
#es = edges(g)
#cppFile = open("B.csv", "w")
#for edge in es
#	s = edge.source.id
#	t = edge.target.id
#
#	write(cppFile, "$(s-1) $(t-1) 1.0\n")
#	write(cppFile, "$(t-1) $(s-1) 1.0\n")
#end
#
#close(cppFile)

# launch the simulation
@time state = simulation(s)

state.T = uniform_phase_shift(state.T)

@info("T_sim: ", state.T[1:20])
@info("T_ref: ", T_out[1:20])

d = chebyshev(state.T, T_out)
@info("distance $d")
@info("# iter: ", state.n_iter)

export_csv_data(state.T, "T_out.csv")

@test_approx_eq_eps d 0. 1e-4




###
# test SD solver
###

@info("############ Testing SD solver")
@info("######## UK grid without dissipation (SD solver)")

p_fn = BASE_FOLDER * "/RK/UK1/P_in.csv"
y_fn = BASE_FOLDER * "/RK/UK1/Y_in.csv"
t_fn = BASE_FOLDER * "/RK/UK1/T_out.csv"

g = load_graph(p_fn,y_fn) 
T_out = uniform_phase_shift(collect(load_csv_data(t_fn)[1]))

o_args = Dict{Symbol,Any}()
o_args[:d] = 1
s = Simulator(g,SD_solver,o_args,1.,1e-10,round(Int64,1e5))

# launch the simulation
@time state = simulation(s)

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

@info("########### Testing NR solver")
@info("####### IEEE benchmark (NR solver)")

sys_fn = BASE_FOLDER * "/NR/IEEE/ieee14cdf.txt"
T_fn = BASE_FOLDER * "/NR/IEEE/T_out.csv"
V_fn = BASE_FOLDER * "/NR/IEEE/V_out.csv"

T_ref = uniform_phase_shift(collect(load_csv_data(T_fn)[1]))
V_ref = collect(load_csv_data(V_fn)[1])

g = load_IEEE_SLFD(sys_fn)

o_args = Dict{Symbol,Any}()
o_args[:g] = g
o_args[:bootstrap_iter] = 0
s = Simulator(g,NR_solver,o_args,100.,1e-10,15)

@time state = simulation(s)

error_V = chebyshev(state.V,V_ref)
# convert angles from radians to degrees
state.T = uniform_phase_shift(state.T*180/pi)
error_T = chebyshev(state.T,T_ref)

@info("Max error in V: $error_V")
@info("Max error in theta: $error_T")
@info("# iter: ", state.n_iter)

@test_approx_eq_eps error_V 0. 1e-4
@test_approx_eq_eps error_T 0. 1e-4

# Test on UK grid without dissipation, same setting as RK_solver and SD_solver
@info("######## UK grid without dissipation (NR solver)")

p_fn = BASE_FOLDER * "/RK/UK1/P_in.csv"
y_fn = BASE_FOLDER * "/RK/UK1/Y_in.csv"
t_fn = BASE_FOLDER * "/RK/UK1/T_out.csv"

g = load_graph(p_fn,y_fn)
# define node 33 as the slack bus
vertices(g)[33].bus_type = 3

T_ref = uniform_phase_shift(collect(load_csv_data(t_fn)[1]))
V_ref = ones(120)

o_args = Dict{Symbol,Any}()
o_args[:g] = g
o_args[:bootstrap_iter] = 0

s = Simulator(g,NR_solver,o_args,1.,1e-10,round(Int64,1e5))

# launch the simulation
@time state = simulation(s)

state.T = uniform_phase_shift(state.T)

#@info("T_sim: ", state.T[1:20])
#@info("T_ref: ", T_out[1:20])

error_V = chebyshev(state.V,V_ref)
error_T = chebyshev(state.T,T_ref)

@info("Max error in V: $error_V")
@info("Max error in theta: $error_T")
@info("# iter: ", state.n_iter)

@test_approx_eq_eps error_V 0. 1e-4
@test_approx_eq_eps error_T 0. 1e-4

# Test on UK grid with dissipation, same setting as for RK_solver
@info("######## UK grid with dissipation (NR solver)")

p_fn = BASE_FOLDER * "/RK/UK2/P_losses.csv"
y_fn = BASE_FOLDER * "/RK/UK2/Y_in2.csv"
t_fn = BASE_FOLDER * "/RK/UK2/theta.csv"

g = load_graph(p_fn,y_fn)
# defining node 33 as the slack bus
vertices(g)[33].bus_type = 3

T_ref = uniform_phase_shift(collect(load_csv_data(t_fn)[1]))
V_ref = ones(120)

o_args = Dict{Symbol,Any}()
o_args[:g] = g
o_args[:bootstrap_iter] = 0
s = Simulator(g,NR_solver,o_args,1.,1e-10,round(Int64,1e5))

# launch the simulation
@time state = simulation(s)

state.T = uniform_phase_shift(state.T)

#@info("T_sim: ", state.T[1:20])
#@info("T_ref: ", T_out[1:20])

error_V = chebyshev(state.V,V_ref)
error_T = chebyshev(state.T,T_ref)

@info("Max error in V: $error_V")
@info("Max error in theta: $error_T")
@info("# iter: ", state.n_iter)

@test_approx_eq_eps error_V 0. 1e-4
@test_approx_eq_eps error_T 0. 1e-4

@info("######## Eurogrid without dissipation (NR solver)")

p_fn = BASE_FOLDER * "/RK/Eurogrid/P_in.csv"
t_fn = BASE_FOLDER * "/RK/Eurogrid/T_out_1e-10.csv"

# load principal component of eurogrid
g = load_jld_serialized("g", "./data/eurogrid/eurogrid_pc.jld") 
@info("# vertices: ", length(vertices(g)))

T_out = uniform_phase_shift(collect(load_csv_data(t_fn)[1]))

# define node 1 as the slack bus
vertices(g)[1].bus_type = 3

P = collect(load_csv_data(p_fn)[1])
change_P(g,P)

o_args = Dict{Symbol,Any}()
o_args[:g] = g
o_args[:bootstrap_iter] = 0
s = Simulator(g,NR_solver,o_args,1.,1e-10,round(Int64,1e5))

# launch the simulation
@time state = simulation(s)

state.T = uniform_phase_shift(state.T)
error_T = chebyshev(state.T,T_out)

@info("T_sim: ", state.T[1:20])
@info("T_ref: ", T_out[1:20])

@info("Max error in theta: $error_T")
@info("# iter: ", state.n_iter)

change_T(g,state.T)

o_args = Dict{Symbol,Any}()
o_args[:h] = 1e-2
s = Simulator(g,RK_solver1,o_args,1.,1e-6,round(Int64,1e6))

@time state = simulation(s)

state.T = uniform_phase_shift(state.T)
error_T = chebyshev(state.T,T_out)

@info("T_sim: ", state.T[1:20])
@info("T_ref: ", T_out[1:20])

@info("Max error in theta: $error_T")
@info("# iter: ", state.n_iter)

#@test_approx_eq_eps error_T 0. 1e-3

###
# test lattice initialization
###

@info("########### Testing lattice initialization")

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

@info("########### Testing contour cycle detection")
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
