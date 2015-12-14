using Graphs

# vertex type
type Bus
	id::Int64
	name::AbstractString
	# 0:PQ, 1:Q \theta, 2: PV, 3: V\theta
	bus_type::Int
	init_voltage::Float64
	final_voltage::Float64
	# base voltage in KV
	base_voltage::Float64
	angle::Float64
	# ENTSOE convention: Re(load) > 0, Re(generation) < 0
	load::Complex{Float64}
	generation::Complex{Float64}
	Q_min::Float64
	Q_max::Float64
	P_min::Float64
	P_max::Float64
	sh_conductance::Float64
	sh_susceptance::Float64
	
	# default constructor
	function Bus(id::Int64, name::AbstractString, bus_type::Int, init_voltage::Float64, final_voltage::Float64, base_voltage::Float64, angle::Float64, load::Complex{Float64}, generation::Complex{Float64}, Q_min::Float64, Q_max::Float64, P_min::Float64, P_max::Float64, sh_conductance::Float64, sh_susceptance::Float64)
		new(id, name, bus_type, init_voltage, final_voltage, base_voltage, angle, load, generation, Q_min, Q_max, P_min, P_max, sh_conductance, sh_susceptance)
	end

	# simple constructor
	function Bus(id::Int64, t::Float64, p::Float64)
		if p > 0
			# load
			new(id, "$id", 2, 1., 1., 1., t, p, 0., 0., 0., 0., 0., 0., 0.)
		else
			# generation
			new(id, "$id", 2, 1., 1., 1., t, 0., p, 0., 0., 0., 0., 0., 0.)
		end
	end 
end

# edge type
type Line
	source::Bus
	target::Bus
	# 0: normal, 1: transformer 
	line_type::Int
	# default on 
	# 0,1,2-> ON
	# 7,8,9 -> OFF  
	line_status::Bool
	# Z = R + iX, Y = 1/Z
	admittance::Complex{Float64}
	sh_susceptance::Float64
	# 1 for a standard power line
	# transformer on the source side (IEEE)
	s_ratio::Complex{Float64}
	# transformer on the target side (ENTSOE)
	t_ratio::Complex{Float64}

	# default constructor
	function Line(source::Bus, target::Bus, line_type::Int, line_status::Bool, admittance::Complex{Float64}, sh_susceptance::Float64, s_ratio::Complex{Float64}, t_ratio::Complex{Float64})
		new(source, target, line_type, line_status, admittance, sh_susceptance, s_ratio, t_ratio)
	end

	# simple constructor
	function Line(source::Bus, target::Bus, y::Complex{Float64})
		new(source, target, 0, true, y, 0., 1., 1.)
	end
end

# redefine some of the base functions
Graphs.source(e::Line, g::Graphs.AbstractGraph{Bus,Line}) = e.source
Graphs.target(e::Line, g::Graphs.AbstractGraph{Bus,Line}) = e.target
Graphs.revedge(e::Line) = Line(e.target, e.source, e.line_type, e.line_status, e.admittance, e.sh_susceptance, e.s_ratio, e.t_ratio)

# get the connected component ids containing the slack bus
function get_slack_component_ids(g::Graphs.AbstractGraph{Bus,Line})
	cs = connected_components(g)
	for c in cs
		if in(3,Int[v.bus_type for v in c])
			return Int[v.id for v in c]
		end
	end
end

# initialize the admittance matrix and the active/reactive  injection vectors
# Sb: base power (for converting in p.u.)
function generate_YPQTV(g::Graphs.AbstractGraph{Bus,Line}, Sb::Float64=100.)
	vs = vertices(g)
	es = edges(g)
	n = length(vs)
	Y = zeros(Complex{Float64}, n,n)

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
	S = Complex{Float64}[-(v.load + v.generation)/Sb for v in vs]
	P = Float64[real(s) for s in S]
	Q = Float64[imag(s) for s in S]
	T = Float64[v.angle for v in vs]
	V = Float64[v.init_voltage for v in vs]

	return Y,P,Q,T,V
end

# generate a double cycle with a bus where p is injected
function generate_double_cycle(l::Int,c:Int,r::Int,p::Float64)
	vs = Bus[]
	es = Line[]

	# left branch 1:l
	for i in 1:l
		# t=0., p=0.
		push!(vs,Bus(i,0.,0.))
	end
	for i in 2:l
		push!(es, Line(vs[i-1],vs[i],1im))
	end
	# central branch (l+1):(c+l)
	for i in (l+1):(c+l)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+2):(c+l)
		push!(es, Line(vs[i-1],vs[i],1im))
	end
	# right branch (l+c+1):(c+l+r)
	for i in (l+c+1):(c+l+r)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+c+2):(c+l+r)
		push!(es, Line(vs[i-1],vs[i],1im))
	end

	b_in = Bus(l+c+r+1,0.,p)
	b_out = Bus(l+c+r+2,0.,-p)

	push!(es, Line(vs[l+c+r+1],vs[1],1im))
	push!(es, Line(vs[l+c+r+1],vs[l+1],1im))
	push!(es, Line(vs[l+c+r+1],vs[l+c+1],1im))
	
	push!(es, Line(vs[l+c+r+2],vs[l],1im))
	push!(es, Line(vs[l+c+r+2],vs[l+c],1im))
	push!(es, Line(vs[l+c+r+2],vs[l+c+r],1im))
	
	return graph(vs, es, is_directed=false)
end
