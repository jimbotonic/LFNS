include("data.jl")
include("solvers.jl")
include("graphs.jl")

using Base.Test

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
s = Simulator(g,SD_solver,o_args,1.,1e-6,round(Int64,1e5))

# launch the simulation
simulation(s)
state = s.states[1]

println("T_sim: ", state.T[1:20])
println("T_ref: ", T_out[1:20])

d = euclidean(state.T, T_out)
println("distance $d")
