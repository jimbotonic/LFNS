
# initialize the admittance matrix and the active/reactive  injection vectors
# Sb: base power (for converting in p.u.)
function init_NR_data(nodes, edges, Sb::Float64=100.)
	n = length(nodes)
	Y = zeros(Complex{Float64}, n,n)

    	for edge in edges
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
	Y = Y + diagm(Float64[node.sh_susceptance for node in nodes])*im
    
	# injections
	S0 = Complex{Float64}[-(n.load + n.generation)/Sb for n in nodes]
	P0 = Float64[real(s) for s in S0]
	Q0 = Float64[imag(s) for s in S0]

	return Y,P0,Q0
end

# load P0 and Y matrix from the specified files 
#
# CSV files with no-header and comma-separated are expected
# P0_fn contains one float per line
# Y_fn: node_id1, node_id2, G_value (0 in the non-dissipative case), B value
function load_RK_data(P0_fn::AbstractString, Y_fn::AbstractString)
	P0_df = load_csv_data(P0_fn)
	Y_df = load_csv_data(Y_fn)

	# the size of the network is assumed to be the # of rows in P0 
	P0 = collect(P0_df[1])
	n = length(P0)
	Y = zeros(Complex{Float64},n,n)
	# TO FINISH !!
	return Y,P0
end
