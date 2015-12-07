# initial conditions
type ICConditions
	P::Array{Float64,1}
	V::Array{Float64,1}
	T::Array{Float64,1}
	Y::Array{Complex{Float64},2}
end

type Changes
	
end

type Store
	T
end

type Simulator
	# initial conditions
	IC::IConditions
	# list of changes
	
	# store
	s::Store	
end
