using Distances

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

# Gauss-Seidel solver for Flow Data networks
# used to initiate NR solver when flat start diverges
function GS_solver(V0::Array{Float64,1}, T::Array{Float64,1}, Y::Array{Complex{Float64},2}, P0::Array{Float64,1}, Q0::Array{Float64,1}, PQ_ids::Array{Int64,1}, slack_id::Int64, iter::Int64=5)
	n = length(V0)
	is_PQ = falses(n,1)
	is_PQ[PQ_ids] = true
	set = setdiff(1:n,slack_id)
	# define the voltages of the non slack buses
	V = V0.*exp(T*im)
	Vnew = V
	for ii in 1:iter
		for i in set
			if(is_PQ[i] == false)
				Q = imag(Y[i,:]*V)
				Vtemp = 1/Y[i,i]*(conj((P0[i]+Q*im)/V[i]) - Y[i,1:i-1]*Vnew[1:i-1] - Y[i,i+1:end]*V[i+1:end])
				Vnew[i] = V0[i]*exp(angle(Vtemp[1])*im)
			else
				Q = Q0[i];
				Vtemp = 1/Y[i,i]*(conj((P0[i]+Q*im)/V[i]) - Y[i,1:i-1]*Vnew[1:i-1] - Y[i,i+1:end]*V[i+1:end])
				Vnew[i] = Vtemp[1]
			end
		end
		V = Vnew
	end
	Vf = zeros(n,1)
	theta = zeros(n,1)
	Vf = abs(Vnew)
	T = angle(Vnew)
	return Vf, T
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
# returns the value of one RK4 iteration + k1 oefficient
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
	M2 = real(Y).*M1 #M2ij=Gij*Vi*Vj
	M3 = imag(Y).*M1 #M3ij=Bij*Vi*Vj
	V1 = diag(M2)
	V2 = cos(T)
	V3 = sin(T)
	# set diagonal elements to 0
	M2 = M2 - diagm(diag(M2))
	M3 = M3 - diagm(diag(M3))

	return (P0 - V1 + (V2.*(-M2*V2 + M3*V3) - V3.*(M2*V3 + M3*V2))) 
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
	nTdot = zeros(Float64, length(T))
	Tdot = zeros(Float64, length(T))

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

# Steepest descent solver for Flow Data networks
#
## INPUT
# T: initial thetas
# Y: admittance square matrix (nxn square matrix of complex numbers)
# P0: initial active powers
#
## OUTPUT
# T: updated thetas
# n_iter: # of iterations before convergence
# delta: norm of the last gradient
function SD_solver(T::Array{Float64,1}, Y::Array{Complex{Float64},2}, P::Array{Float64,1}, epsilon::Float64=1e-6, iter_max::Int64=1e4, del::FLoat64=1e-2)
	n_iter = 0
	n = length(T)

	# We only use the susceptive part of the admittance matrix, with zero diagonal elements	
	K = imag(Y-diagm(diag(Y)))
	
	dT = T*ones(1,n)-ones(n,1)*T'
	
	# f0 is the potential in the phase space whose extremas are solutions of the PFEs
	f0 = sum(P.*T) + .5*sum(K.*cos(dT))
	
	# nabla is the gradient of this potential
	nabla = P - sum(K.*sin(dT),2)
	delta = norm(nabla)
	
	# Follow the path of highest gradient, until the correction is less than delta, to reach the top of a hill
	while delta > epsilon && n_iter < iter_max
		n_iter += 1
		A = copy(T)
		T += nabla*del
		f1 = sum(P.*T) + .5*sum(K.*cos(T*ones(1,n)-ones(n,1)*T'))
		while f1 < f0 && norm(f1-f0) > epsilon
			T = A
			del = del/2
			T += nabla*del
			f1 = sum(P.*T) + .5*sum(K.*cos(T*ones(1,n)-ones(n,1)*T'))
		end
		dT = T*ones(1,n)-ones(n,1)*T'
		f0 = f1
		nabla = P - sum(K.*sin(dT),2)
		delta = norm(nabla)
	end
	
	# rotate all angles by using the last one as the reference
	# all angles belong to [-pi,pi] afterward
	T = mod(T-T[end]+pi,2*pi)-pi
	
	return T,n_iter,delta
end
