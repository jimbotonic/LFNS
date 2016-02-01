include("graphs.jl")

abstract  Change

# bus change
abstract BChange <: Change
# line change
abstract LChange <: Change

# change powers
type PChange <: BChange
	# array of new values
	nvs::Array{Float64,1}	
end

# change of angle 
type TChange <: BChange
	# array of new values
	nvs::Array{Float64,1}	
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
	o_args::Dict{Symbol,Any}

	# default constructor
	function SParams(V::Array{Float64,1},T::Array{Float64,1},Y::Array{Complex{Float64},2},P::Array{Float64,1},Q::Array{Float64,1},epsilon::Float64,iter_max::Int64,o_args::Dict{Symbol,Any})
		return new(V,T,Y,P,Q,epsilon,iter_max,o_args)
	end
end 

# state of the system at a given index (or time)
type State
	V::Array{Float64,1}
	T::Array{Float64,1}
	Tdot::Array{Float64,1}
	n_iter::Int64	
	# state optional data
	o_data::Dict{Symbol,Any}
	
	# default constructor
	function State(V::Array{Float64,1},T::Array{Float64,1},Tdot::Array{Float64,1},n_iter::Int64,o_data::Dict{Symbol,Any})
		return new(V,T,Tdot,n_iter,o_data)
	end
	
	function State(V::Array{Float64,1},T::Array{Float64,1},Tdot::Array{Float64,1},n_iter::Int64)
		return new(V,T,Tdot,n_iter,Dict{Symbol,Any}())
	end
end

type Simulator
	# graph being studied (specifiy implicitly Y,P,Q,V,T) 
	g::Graphs.AbstractGraph{Bus,Line}
	# solver
	solver::Function
	# solver additional optional arguments
	o_args::Dict{Symbol,Any}
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
	function Simulator(g::Graphs.AbstractGraph{Bus,Line},solver::Function,o_args::Dict{Symbol,Any},sb::Float64,epsilon::Float64,iter_max::Int64,changes::Array{Change,1},states::Array{State,1})
		return new(g,solver,o_args,sb,epsilon,iter_max,changes,states)
	end
	
	function Simulator(g::Graphs.AbstractGraph{Bus,Line},solver::Function,o_args::Dict{Symbol,Any},sb::Float64,epsilon::Float64,iter_max::Int64)
		return new(g,solver,o_args,sb,epsilon,iter_max,Array{Change,1}(),Array{State,1}())
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

# initialize the solvers parameters (admittance matrix, active/reactive  injection vectors, angles, ...)
# Sb: base power (for converting in p.u.)
function get_sparams(s::Simulator)
	vs = vertices(s.g)
	es = edges(s.g)
	n = length(vs)
	#Y = zeros(Complex{Float64},n,n)
	Y = spzeros(Complex{Float64},n,n)

    	for edge in es
        	if edge.line_status
			y = edge.admittance
			if edge.line_type == 0
				y_sh = edge.sh_susceptance*im
				Y[edge.source.id, edge.source.id] += y + y_sh/2
				Y[edge.target.id, edge.target.id] += y + y_sh/2
				Y[edge.source.id, edge.target.id] += -y
				Y[edge.target.id, edge.source.id] += -y
			elseif edge.line_type == 1
				#y_sh = (edge.sh_conductance + edge.sh_susceptance*im)
				y_sh = edge.sh_susceptance*im
				Y[edge.source.id, edge.source.id] += y*abs(edge.s_ratio)^2 
				Y[edge.target.id, edge.target.id] += abs(edge.t_ratio)^2*(y + y_sh)
				Y[edge.source.id, edge.target.id] += -y*edge.s_ratio*edge.t_ratio
				Y[edge.target.id, edge.source.id] += -y*conj(edge.s_ratio*edge.t_ratio)
			end
        	end
    	end

	# add shunt susceptance to each node
	Y = Y + diagm(Float64[v.sh_susceptance for v in vs])*im
    
	# injections
	S = Complex{Float64}[(-v.load + v.generation)/s.sb for v in vs]
	P = Float64[real(s) for s in S]
	Q = Float64[imag(s) for s in S]
	T = Float64[v.angle for v in vs]
	V = Float64[v.init_voltage for v in vs]

	return SParams(V,T,Y,P,Q,s.epsilon,s.iter_max,s.o_args)
end

# change the value of P
function change_P(g::Graphs.AbstractGraph{Bus,Line},P::Array{Float64,1})
	vs = vertices(g)
	for i in 1:length(vs)
		if P[i] >= 0
			vs[i].generation = P[i]
		else
			vs[i].load = -P[i]
		end
	end
end
