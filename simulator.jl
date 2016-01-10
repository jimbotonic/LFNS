include("graphs.jl")

abstract  Change

# bus change
abstract BChange <: Change
# line change
type LChange <: Change

# change powers
type PChange <: BChange
	# index -> new value
	nvs::Dict{Int,Float64}	
end

# change of angle 
type TChange <: BChange
	# index -> new value
	nvs::Dict{Int,Float64}	
end

# change of admittance
type YChange <: LChange
	i::Int
	j::Int
	nv::Complex{Float64}
end

# solver parameters
type SParams
	# powers
	V::Array{Float64,1}
	# angles
	T::Array{Float64,1}
	# admittance matrix
	Y::Array{Complex{Float64},2}
	P::Array{Float64,1}
	Q::Array{Float64,1}
	# convergence criteria
	epsilon::Float64
	iter_max::Int64	
	# solver additional optional arguments
	o_args::Dict{AbstractString,Any}

	# default constructor
	function SParams(V::Array{Float64,1},T::Array{Float64,1},Y::Array{Complex{Float64},2},P::Array{Float64,1},Q::Array{Float64,1},epsilon::Float64,iter_max::Int64,o_args::Dict{AbstractString,Any})
		return new SParams(V,T,Y,P,Q,epsilon,iter_max,o_args)
	end
end 

# state of the system at a given index (or time)
type State
	V::Array{Float64,1}
	T::Array{Float64,1}
	Tdot::Array{Float64,1}
	n_iter::Int64	
	# state optional data
	o_data::Dict{AbstractString,Any}
	
	# default constructor
	function State(V::Array{Float64,1},T::Array{Float64,1},Tdot::Array{Float64,1},n_iter::Int64,o_data::Dict{AbstractString,Any})
		return new State(V,T,Tdot,n_iter,o_data)
	end
	
	function State(V::Array{Float64,1},T::Array{Float64,1},Tdot::Array{Float64,1},n_iter::Int64)
		return new State(V,T,Tdot,n_iter,Dict{AbstractString,Any}())
	end
end

type Simulator
	# graph being studied (specifiy implicitly Y,P,Q,V,T) 
	g::Graphs.AbstractGraph{Bus,Line}
	# solver
	solver::Function
	# solver additional optional arguments
	o_args::Dict{AbstractString,Any}
	# base line voltage
	sb::Float64
	# convergence criteria
	epsilon::Float64
	iter_max::Int64	
	# list of changes
	changes::Array{Change,1}	
	# states history
	states::Array{State,1}
	
	# default constructor
	function Simulator(g::Graphs.AbstractGraph{Bus,Line},solver::Function,o_args::Dict{AbstractString,Any},sb::Float64,epsilon::Float64,iter_max::Int64,changes::Array{Change,1},states::Array{State,1})
		return new Simulator(g,solver,o_args,sb,epsilon,iter_max,changes,states)
	end
	
	function Simulator(g::Graphs.AbstractGraph{Bus,Line},solver::Function,o_args::Dict{AbstractString,Any},sb::Float64,epsilon::Float64,iter_max::Int64)
		return new Simulator(g,solver,o_args,sb,epsilon,iter_max,Array{Change,1}(),Array{State,1}())
	end
end


# lead a simulation
function simulation(s::Simulator)
	# initialize simulation data
	sp = get_sparams(s)
	state = s.solver(sp)
	push!(s.states,state)
	# iterate over all changes and store state at each iteration
end
