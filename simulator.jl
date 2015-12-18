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

# state of the system at a given index
type State
	i::Int
	V::Array{Float64,1}
	T::Array{Float64,1}
	Tdot::Array{Float64,1}
end

type Simulator
	# initial conditions (specifiy implicitly Y,P,Q,V,T) 
	g::Graphs.AbstractGraph{Bus,Line}
	# solver
	solver::Function
	# solver additional arguments
	s_args::Dict{AbstractString,Any}
	# convergence criteria
	epsilon::Float64
	iter_max::Int64	
	# list of changes
	changes::Array{Change,1}	
	# states history
	states::Array{State,1}
end
