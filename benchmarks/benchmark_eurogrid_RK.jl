using Logging, ConfParser, Graphs, Distributions, PyPlot

include("../../graphs.jl")
include("../../data.jl")
include("../../simulator.jl")
include("../../solvers.jl")
include("../../metrics.jl")
include("../../init.jl")

@Logging.configure(level=INFO)

epsilon = 1e-7

g0 = load_jld_serialized("g","./data/eurogrid/eurogrid.jld")
b_pc = get_principal_component(g)
vids = Array{Int64,1}()
for b in b_pc
	push!(vids,b.id)
end
g = get_subgraph(g0,vids,true)

P = readcsv("g","./data/eurogrid/eurogrid_P_benchmark.dat")

o_args = Dict{Symbol,Any}()
o_args[:h] = 1e-2

s = Simulator(g,RK_solver1,o_args,1.,epsilon,round(Int64,1e6))

tic()
state = simulation(s)
toc()


