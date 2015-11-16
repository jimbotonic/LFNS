# initialize the admittance matrix and the active/reactive  injection vectors
function init_simulation_data(nodes, edges, Sb::Float64=100.)
	n = length(nodes)
	Y = zeros(Complex{Float64}, n,n)

    	for edge in edges
        	if edge.line_status
			z = edge.resistance + edge.reactance*im
			if edge.line_type == 0
				y = 1e-6*edge.sh_susceptance*im
				Y[edge.source_id, edge.source_id] += 1/z + y/2
				Y[edge.target_id, edge.target_id] += 1/z + y/2
				Y[edge.source_id, edge.target_id] += -1/z
				Y[edge.target_id, edge.source_id] += -1/z
			elseif edge.line_type == 1
				#y = 1e-6*(edge.sh_conductance + edge.sh_susceptance*im)
				y = 1e-6*edge.sh_susceptance*im
				Y[edge.source_id, edge.source_id] += 1/z
				Y[edge.target_id, edge.target_id] += (1/z)*abs(edge.ratio)^2 + y
				Y[edge.source_id, edge.target_id] += -(1/z)*edge.ratio
				Y[edge.target_id, edge.source_id] += -(1/z)*conj(edge.ratio)
			end
        	end
    	end
    
	# injections
	S0 = Complex{Float64}[-(n.load + n.generation)/Sb for n in nodes]
	P0 = Float64[real(s) for s in S0]
	Q0 = Float64[imag(s) for s in S0]

	return Y,P0,Q0
end

# Newton-Raphson solver for Flow Data networks
#
## INPUT
# Y: admittance square matrix (nxn square matrix)
# V: inital voltages
# T: initial thetas
# P0: initial active powers
# Q0: initial reactive powers
# PQ_ids: array of PQ bus ids (m-dimensional vector with m<n)
# slack_id: slack bus id
#
## OUTPUT
# V: updated voltages
# T: updated thetas
# n_iter: # of iterations before convergence
function NR_solver(Y, V, T, P0, Q0, PQ_ids, slack_id, epsilon::Float64=1e-4, iter_max::Int=50)
    	n = size(Y)[1]
	# compute Y element-wise absolute values
    	Y_abs = abs(Y)
	# compute Y element-wise angle values
    	Y_angle = angle(Y)
	# ids from 1:n except slack_id
	ids = collect(1:n)
	deleteat!(ids, slack_id)
    	
	error  = epsilon
    	n_iter = 1

    	while(error >= epsilon && n_iter < iter_max)
		# pre-computations 1
		M1 = V*V'*Y_abs 
		M2 = repmat(T,1,n)-repmat(T',n,1)-Y_angle
		M3 = M1*sin(M2)
		M4 = M1*cos(M2)
		V1 = diag(Y_abs)*sin(diag(Y_angle))*V^2
		V2 = diag(Y_abs)*cos(diag(Y_angle))*V^2

		# compute P and Q (n-dimensional vectors)
		# 
		# S_i = P_i + j Q_i
		# P_i = Re[Vi (\sum Y_{ij}^* V_j^*)]
		P = sum(M4,2)
		Q = sum(M3,2)

		# pre-computations 2
		V3 = -Q-V1
		V4 = Q-V1
		V5 = P+V2
		V6 = P-V2
		
		# compare P and Q with their "known" value to get dPQ ((n-1+m)-dimensional vector)
		dP = P0[ids] - P[ids]
		dQ = Q0[PQ_ids] - Q[PQ_ids]
		dPQ = [dP;dQ]
		
		# compute the jacobian matrix
		# see Power System Analysis, Bergen/Vital, p345-347
		# L: (n-1)x(n-1), O: mxm, N: (n-1)xm M: mx(n-1), J:(n-1+m)x(n-1xm)
		L = M3[ids,ids] - diagm(diag(M3[ids,ids])) + diagm(V3[ids])
		O = M3[PQ_ids,PQ_ids] - diagm(diag(M3[PQ_ids,PQ_ids])) + diagm(V4[PQ_ids])
		N = M4 - diagm(diag(M4)) + diagm(V5)
		N = N[ids,PQ_ids]
		M = -M4 - diagm(diag(-M4)) + diagm(V6)
		M = M[PQ_ids,ids]
		J = [L N;M O]
		# solve JX = dPQ
		X = J\dPQ
		
		# update V and theta
		T[ids] = T[ids] + X[1:n-1]
		V[PQ_ids] = V[PQ_ids] + X[n:end]
		error = max(abs(dPQ))
		n_iter = n_iter + 1
    	end

	return V,T,n_iter
end
