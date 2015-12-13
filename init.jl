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
			z = edge.resistance + edge.reactance*im
			if edge.line_type == 0
				y = edge.sh_susceptance*im
				Y[edge.source.id, edge.source.id] += 1/z + y/2
				Y[edge.target.id, edge.target.id] += 1/z + y/2
				Y[edge.source.id, edge.target.id] += -1/z
				Y[edge.target.id, edge.source.id] += -1/z
			elseif edge.line_type == 1
				#y = (edge.sh_conductance + edge.sh_susceptance*im)
				y = edge.sh_susceptance*im
				Y[edge.source.id, edge.source.id] += (1/z)*abs(edge.s_ratio)^2 
				Y[edge.target.id, edge.target.id] += abs(edge.t_ratio)^2*(1/z + y)
				Y[edge.source.id, edge.target.id] += -(1/z)*edge.s_ratio*edge.t_ratio
				Y[edge.target.id, edge.source.id] += -(1/z)*conj(edge.s_ratio*edge.t_ratio)
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



