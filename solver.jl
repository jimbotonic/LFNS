# Newton-Raphson solver for Flow Data networks
#
## INPUT
# Y: admittance square matrix
# V: inital voltages
# T: initial thetas
# P0: initial active powers
# Q0: initial reactive powers
# PQ_ids: array of PQ buses ids
# slack_id: slack bus id
#
## OUTPUT
# V: updated voltages
# T: updated thetas
# n_iter: # of iterations before convergence
function NR_solver(Y, V, T, P0, Q0, PQ_ids, slack_id, epsilon::Float=1e-4, iter_max::Int=50)
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

		# compute P and Q
		P = sum(M4,2)
		Q = sum(M3,2)

		# pre-computations 2
		V3 = -Q-V1
		V4 = Q-V1
		V5 = P+V2
		V6 = P-V2
		
		# compare P and Q with their "known" value
		dP = P0[ids] - P[ids]
		dQ = Q0[PQ_ids] - Q[PQ_ids]
		dPQ = [dP;dQ]
		
		# compute the jacobian matrix
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
