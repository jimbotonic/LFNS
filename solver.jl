using Distances

# initialize the admittance matrix and the active/reactive  injection vectors
function init_NR_data(nodes, edges, Sb::Float64=100.)
	n = length(nodes)
	Y = zeros(Complex{Float64}, n,n)

    	for edge in edges
        	if edge.line_status
			z = edge.resistance + edge.reactance*im
			if edge.line_type == 0
				y = edge.sh_susceptance*im
				Y[edge.source_id, edge.source_id] += 1/z + y/2
				Y[edge.target_id, edge.target_id] += 1/z + y/2
				Y[edge.source_id, edge.target_id] += -1/z
				Y[edge.target_id, edge.source_id] += -1/z
			elseif edge.line_type == 1
				#y = (edge.sh_conductance + edge.sh_susceptance*im)
				y = edge.sh_susceptance*im
				Y[edge.source_id, edge.source_id] += (1/z)*abs(edge.ratio)^2 + y
				Y[edge.target_id, edge.target_id] += 1/z
				Y[edge.source_id, edge.target_id] += -(1/z)*edge.ratio
				Y[edge.target_id, edge.source_id] += -(1/z)*conj(edge.ratio)
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

# Newton-Raphson solver for Flow Data networks
#
## INPUT
# V: inital voltages
# T: initial thetas
# Y: admittance square matrix (nxn square matrix of complex numbers)
# P0: initial active powers
# Q0: initial reactive powers
# PQ_ids: array of PQ bus ids (m-dimensional vector with m<n)
# slack_id: slack bus id
#
## OUTPUT
# V: updated voltages
# T: updated thetas
# n_iter: # of iterations before convergence
function NR_solver(V::Array{Float64,1}, T::Array{Float64,1}, Y::Array{Complex{Float64},2}, P0::Array{Float64,1}, Q0::Array{Float64,1}, PQ_ids::Array{Int64,1}, slack_id::Int64, epsilon::Float64=1e-6, iter_max::Int64=50)
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
		M1 = V*V'.*Y_abs 
		M2 = repmat(T,1,n)-repmat(T',n,1)-Y_angle
		M3 = M1.*sin(M2)
		M4 = M1.*cos(M2)
		V1 = diag(Y_abs).*sin(diag(Y_angle)).*V.^2
		V2 = diag(Y_abs).*cos(diag(Y_angle)).*V.^2

		# compute P and Q (n-dimensional vectors)
		# 
		# S_i = P_i + j Q_i
		# P_i = Re[Vi (\sum Y_{ij}^* V_j^*)]
		P = vec(sum(M4,2))
		Q = vec(sum(M3,2))

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
		# NB: to keep only the diagonal elements, use diagm(collect(diag(M))) OR triu(tril(M))
		N = M4 - diagm(diag(M4)) + diagm(V5)
		N = N[ids,PQ_ids]
		M = -M4 - diagm(diag(-M4)) + diagm(V6)
		M = M[PQ_ids,ids]
		J = [L N;M O]
		# solve JX = dPQ
		X = J\dPQ
		
		# update V and theta
		T[ids] += X[1:n-1]
		V[PQ_ids] += X[n:end]
		error = maximum(abs(dPQ))
		n_iter += 1
    	end

	return V,T,n_iter
end

# RK4 Runge-Kutta method (used to solve y_dot = f(t,y), y(t_0) = y_0)
#
# f: function of t and y
# t: time
# y: function of t
# h: step size
#
# time dependent version
#function RK4(f)
#        return   (t, y, h)-> 
#               ( (k1   )-> 
#               ( (k2   )-> 
#               ( (k3   )-> 
#               ( (k4   )->( k1 + 2*k2 + 2*k3 + k4 ) / 6 
#               )( h * f( t + h  , y + k3   ) )
#               )( h * f( t + h/2, y + k2/2 ) )
#               )( h * f( t + h/2, y + k1/2 ) )
#               )( h * f( t      , y        ) )
#end
#
# varargs time-independent version
function RK4(f)
        return   (y, h, v...)-> 
               ( (k1   )-> 
               ( (k2   )-> 
               ( (k3   )-> 
               ( (k4   )-> (( k1 + 2*k2 + 2*k3 + k4 )/6, k1) 
               )( h * f( y + k3  , v... ) )
               )( h * f( y + k2/2, v... ) )
               )( h * f( y + k1/2, v... ) )
               )( h * f( y       , v... ) )
end

# right-hand side of the differential equation 
function f1(T::Array{Float64,1}, V::Array{Float64,1}, Y::Array{Complex{Float64},2}, P0::Array{Float64,1})
	M1 = V*V'
	M2 = G.*M1
	M3 = B.*M1
	V1 = diag(M2)
	V2 = cos(T)
	V3 = sin(T)
	# set diagonal elements to 0
	M2 = M2 - diagm(diag(M2))
	M3 = M3 - diagm(diag(M3))

	return (P0 - V1 + (V2.*(-M2*V2 + M3*V3) - v3.*(M2*V3 + M3*V2))) 
end

# Runge-Kutta solver for Flow Data networks
#
## INPUT
# T: initial thetas
# h: step size (default value 1e-2)
# Tdot: initial theta derivatives (often vector of 0s, i.e., "flat start")
# V: initial voltages
# Y: admittance square matrix (nxn square matrix of complex numbers)
# P0: initial active powers
#
## OUTPUT
# T: updated thetas
# Tdot: updated theta_dots
# n_iter: # of iterations before convergence
function RK_solver1(T::Array{Float64,1}, h::Float64, V::Array{Float64,1}, Y::Array{Complex{Float64},2}, P0::Array{Float64,1}, epsilon::Float64=1e-6, iter_max::Int64=1e4)
	dU = RK4(f1)
	nTdot = zeros(Float64, length(n))

	n_iter = 1
	while n_iter < iter_max
		(dT, Tdot) = dU(T, h, V, Y, P0)
		if chebyshev(nTdot, Tdot) < epsilon
			break
		end
		nTdot = copy(Tdot)
		T += dT
		n_iter += 1
	end

	return T,Tdot,n_iter
end
