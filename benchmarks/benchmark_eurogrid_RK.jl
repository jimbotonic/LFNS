using Logging, ConfParser, Graphs, Distributions, PyPlot

include("../graphs.jl")
include("../data.jl")
include("../simulator.jl")
include("../solvers.jl")
include("../metrics.jl")
include("../init.jl")

@Logging.configure(level=INFO)

epsilon = 1e-7
max_iter = round(Int64,1e2)

g0 = load_jld_serialized("g","../data/eurogrid/eurogrid.jld")
b_pc = get_principal_component(g0)
vids = Array{Int64,1}()
for b in b_pc
	push!(vids,b.id)
end
g = get_subgraph(g0,vids,true)

P = collect(readcsv("../data/eurogrid/eurogrid_P_benchmark.dat"))
change_P(g,P)

o_args = Dict{Symbol,Any}()
o_args[:h] = 1e-2

s = Simulator(g,RK_solver1,o_args,1.,epsilon,max_iter)

tic()
state = simulation(s)
toc()


