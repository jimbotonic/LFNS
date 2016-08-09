include("thermic.jl")

# thermic solver parameters
type SParams
	# set of nearby houses
	neighborhood::Array{House,1}
	# convergence criteria
	iter_max::Int64	
	# solver additional optional arguments
	o_args::Dict{Symbol,Any}

	# default constructor
	function SParams(neighborhood::Array{House,1},iter_max::Int64,o_args::Dict{Symbol,Any})
		return new(neighborhood,iter_max,o_args)
	end
end 

# state of the system at a given index (or time)
type State
	# house powers
	P::Array{Float64,1}
	# temperatures (heatpump or boiler)
	K::Array{Float64,1}
	# state optional data
	o_data::Dict{Symbol,Any}
	
	# default constructor
	function State(P::Array{Float64,1},K::Array{Float64,1},o_data::Dict{Symbol,Any})
		return new(P,K,o_data)
	end
end

type Simulator
	# set of houses
	neighborhood::Array{House,1}
	# solver
	solver::Function
	# solver additional optional arguments
	o_args::Dict{Symbol,Any}
	# convergence criteria
	iter_max::Int64	
	# states history
	states::Array{State,1}
	
	# default constructor
	function Simulator(neighborhood::Array{House,1},solver::Function,o_args::Dict{Symbol,Any},iter_max::Int64,states::Array{State,1})
		new(neighborhood,solver,o_args,iter_max,states)
	end	
end

# lead a simulation 
#
# callback_func: callback function to be called at each iteration of the solver
function simulation(s::Simulator,callback_func::Function)
	sp = get_sparams(s)
	state = s.solver(sp,callback_func)
	push!(s.states,state)
	return state
end

# initialize the solvers parameters
function get_sparams(s::Simulator)
	return SParams(s.neighborhood,s.iter_max,s.o_args)
end
