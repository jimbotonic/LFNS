include("graphs.jl")

# TODO: see generate_PY in graphs.jl
#
# initialize the admittance matrix and the active/reactive  injection vectors
# Sb: base power (for converting in p.u.)
function init_NR_data(g::Graphs.AbstractGraph{Bus,Line}, Sb::Float64=100.)
	vs = vertices(g)
	es = edges(g)
	n = length(vs)
	Y = zeros(Complex{Float64}, n,n)

    	for edge in es
        	if edge.line_status
			#z = edge.resistance + edge.reactance*im
			y = edge.admittance
			if edge.line_type == 0
				ysh = edge.sh_susceptance*im
				Y[edge.source.id, edge.source.id] += y + ysh/2
				Y[edge.target.id, edge.target.id] += y + ysh/2
				Y[edge.source.id, edge.target.id] += -y
				Y[edge.target.id, edge.source.id] += -y
			elseif edge.line_type == 1
				#y = (edge.sh_conductance + edge.sh_susceptance*im)
				y = edge.sh_susceptance*im
				Y[edge.source.id, edge.source.id] += y*abs(edge.s_ratio)^2 
				Y[edge.target.id, edge.target.id] += abs(edge.t_ratio)^2*(y + ysh)
				Y[edge.source.id, edge.target.id] += -y*edge.s_ratio*edge.t_ratio
				Y[edge.target.id, edge.source.id] += -y*conj(edge.s_ratio*edge.t_ratio)
			end
        	end
    	end

	# add shunt susceptance to each node
	Y = Y + diagm(Float64[v.sh_susceptance for v in vs])*im
    
	# injections
	S0 = Complex{Float64}[-(v.load + v.generation)/Sb for v in vs]
	P0 = Float64[real(s) for s in S0]
	Q0 = Float64[imag(s) for s in S0]

	return Y,P0,Q0
end



