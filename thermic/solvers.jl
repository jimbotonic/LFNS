using Distances

include("simulator.jl")

## INPUT
# callback_func: callback function to be called at each iteration of the solver
#
## OUTPUT
function Euler_boiler_solver(sp::SParams,callback_func::Function)
	# return State(Float64[],sp.T,Tdot,n_iter)
end

## INPUT
# callback_func: callback function to be called at each iteration of the solver
#
## OUTPUT
function Euler_heatpump_solver(sp::SParams,callback_func::Function)
	# return State(Float64[],sp.T,Tdot,n_iter)
end
