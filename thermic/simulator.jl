include("thermic.jl")

# state of the system at a given index (or time)
type State
	# thermic device powers
	P::Array{Float64,1}
	# temperatures (heatpump or boiler)
	K::Array{Float64,1}
	# state optional data
	o_data::Dict{Symbol,Any}

	# default constructor
	function State(P::Array{Float64,1},K::Array{Float64,1},o_data::Dict{Symbol,Any})
		return new(P,K,o_data)
	end
	
	function State(P::Array{Float64,1},K::Array{Float64,1})
		return new(P,K,Dict{Symbol,Any}())
	end
end

# thermic solver parameters
type SParams
	# set of nearby houses
	neighborhood::Array{House,1}
	# controller function
	controller::Function
	# solver additional optional arguments
	o_args::Dict{Symbol,Any}
	# states history
	states::Array{State,1}
	# convergence criteria
	iter_max::Int64	
	# time step (time is measured in hour by default and delta is equal to 1 minute (i.e. 1/60 [h]))
	delta::Float64

	# default constructor
	function SParams(neighborhood::Array{House,1},controller::Function,o_args::Dict{Symbol,Any},iter_max::Int64,delta::Float64=1/60)
		return new(neighborhood,controller,o_args,State[],iter_max,delta)
	end
end 

type Simulator
	# set of houses
	neighborhood::Array{House,1}
	# solver & controller functions
	solver::Function
	controller::Function
	# solver additional optional arguments
	o_args::Dict{Symbol,Any}
	# convergence criteria
	iter_max::Int64	
	# time step (time is measured in hour by default and delta is equal to 1 minute (i.e. 1/60 [h]))
	delta::Float64
	
	# default constructor
	function Simulator(neighborhood::Array{House,1},solver::Function,controller::Function,o_args::Dict{Symbol,Any},iter_max::Int64,delta::Float64=1/60)
		new(neighborhood,solver,controller,o_args,iter_max,delta)
	end	
end

# lead a simulation 
function simulation(s::Simulator)
	sp = get_sparams(s)
	s.solver(sp)
	return sp.states
end

# initialize the solvers parameters
function get_sparams(s::Simulator)
	return SParams(s.neighborhood,s.controller,s.o_args,s.iter_max,s.delta)
end
