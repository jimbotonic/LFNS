include("data.jl")
include("solvers.jl")
include("graphs.jl")

using Base.Test, Logging

@Logging.configure(level=INFO)

BASE_FOLDER = "./data/tests"
# test RK solver

###
# test RG solver
###

p_fn = BASE_FOLDER * "/RK/UK1/P_in.csv"
y_fn = BASE_FOLDER * "/RK/UK1/Y_in.csv"
t_fn = BASE_FOLDER * "/RK/UK1/T_out.csv"

g = load_graph(p_fn,y_fn) 
T_out = collect(load_csv_data(t_fn)[1])

o_args = Dict{Symbol,Any}()
o_args[:d] = 1e-2
s = Simulator(g,SD_solver,o_args,1.,1e-6,round(Int64,1e1))

# launch the simulation
simulation(s)
state = s.states[1]

@info("T_sim: ", state.T[1:20])
@info("T_ref: ", T_out[1:20])

d = euclidean(state.T, T_out)
@info("distance $d")
@info("# iter: ", state.n_iter)

#@test_approx_eq_eps d 0. 1e-4


###
# test NR solver
###

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

simulation(s)
state = s.states[1]

error_V = maximum(abs(state.V-V_ref))
state.T = state.T*180/pi
error_T = maximum(abs(state.T-T_ref))
@info("Max error in V: $error_V")
@info("Max error in theta: $error_T")

@test_approx_eq_eps error_V 0. 1e-6
@test_approx_eq_eps error_T 0. 1e-6
