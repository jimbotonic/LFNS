include("graphs.jl")

abstract  Change

# vertex change
type VChange <: Change
	# index -> new value
	nvs::Dict{Int,Float64}	
end

# edge change
type EChange <: Change
	i::Int
	j::Int
	nv::Float64
end 

# state of the system at a given index
type State
	i::Int
	V::Array{Float64,1}
	T::Array{Float64,1}
	Tdot::Array{Float64,1}
end

type Simulator
	# initial conditions (specifiy P,V,T,Y) 
	g::Graphs.AbstractGraph{Bus,Line}
	# convergence criteria
	epsilon::Float64
	iter_max::Int64	
	# list of changes
	changes::Array{Change,1}	
	# states history
	states::Array{State,1}
end
