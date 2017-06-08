using Distances

include("simulator.jl")

# shift vector angles by using the first one as the reference
# all angles belong to [-pi,pi] afterward
function uniform_phase_shift(T::Array{Float64,1})
	return mod(T-T[1]+pi,2*pi)-pi
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
function NR_solver(sp::SParams)
	vs = vertices(sp.o_args[:g])
	n = length(vs) 
	
	# node ids whose bus type is 0
	PQ_ids = Int64[v.id for v in filter(v -> v.bus_type == 0, vs)]

	# bus types of the slack component 	
	bus_type = Int64[v.bus_type for v in vs]
	# PQ bus positions in the slack component
	PQ_pos = findin(bus_type, 0)
	# slack position in the slack component
	slack_pos = findin(bus_type,3)[1]
	
	# bootstrap simulation
	sp.o_args[:PQ_pos] = PQ_pos
	sp.o_args[:slack_pos] = slack_pos
	GS_solver(sp)

	# compute Y element-wise absolute values
    	Y_abs = abs(sp.Y)
	# compute Y element-wise angle values
    	Y_angle = angle(sp.Y)
	
	# node id whose bus type is 3
	slack_id = filter(v -> v.bus_type == 3, vs)[1].id
	# ids from 1:n except slack_id
	ids = collect(1:n)
	deleteat!(ids, slack_id)
    	
	error  = sp.epsilon
    	n_iter = 1

    	while(error >= sp.epsilon && n_iter < sp.iter_max)
		# pre-computations 1
		M1 = sp.V*sp.V'.*Y_abs 
		M2 = repmat(sp.T,1,n)-repmat(sp.T',n,1)-Y_angle
		M3 = M1.*sin(M2)
		M4 = M1.*cos(M2)
		V1 = diag(Y_abs).*sin(diag(Y_angle)).*sp.V.^2
		V2 = diag(Y_abs).*cos(diag(Y_angle)).*sp.V.^2

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
		dP = sp.P[ids] - P[ids]
		dQ = sp.Q[PQ_ids] - Q[PQ_ids]
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
		sp.T[ids] += X[1:n-1]
		sp.V[PQ_ids] += X[n:end]
		error = norm(dPQ,Inf)
		n_iter += 1
    	end

	return State(sp.V,sp.T,Float64[],n_iter)
end

# Gauss-Seidel solver for Flow Data networks
# used to initiate NR solver when flat start diverges
function GS_solver(sp::SParams)
	n = length(sp.V)
	is_PQ = falses(n,1)
	is_PQ[sp.o_args[:PQ_pos]] = true
	set = setdiff(1:n,sp.o_args[:slack_pos])
	# define the voltages of the non slack buses
	V = sp.V.*exp(sp.T*im)
	# copy V
	Vnew = copy(V)
	for ii in 1:sp.o_args[:bootstrap_iter]
		for i in set
			# PV bus
			if (is_PQ[i] == false)
				Q = imag(sp.Y[i,:]*V)
				Vtemp = 1/sp.Y[i,i]*(conj((sp.P[i]+Q*im)/V[i]) - sp.Y[i,1:i-1]*Vnew[1:i-1] - sp.Y[i,i+1:end]*V[i+1:end])
				Vnew[i] = V[i]*exp(angle(Vtemp[1])*im)
			# PQ bus
			else
				Vtemp = 1/sp.Y[i,i]*(conj((sp.P[i]+sp.Q[i]*im)/V[i]) - sp.Y[i,1:i-1]*Vnew[1:i-1] - sp.Y[i,i+1:end]*V[i+1:end])
				Vnew[i] = Vtemp[1]
			end
		end
		# copy Vnew
		V = copy(Vnew)
	end
	# update powers and angles vectors
	sp.V = abs(Vnew)
	sp.T = angle(Vnew)
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
               ( (k4   )-> (( k1 + 2*k2 + 2*k3 + k4 )/6, k1/h) 
               )( h * f( y + k3  , v... ) )
               )( h * f( y + k2/2, v... ) )
               )( h * f( y + k1/2, v... ) )
               )( h * f( y       , v... ) )
end

# right-hand side of the differential equation 
function f1(T::Array{Float64,1}, V1::Array{Float64,1}, M1::SparseMatrixCSC{Float64,Int64}, M2::SparseMatrixCSC{Float64,Int64}, P::Array{Float64,1})
	V2 = cos(T)
	V3 = sin(T)
	
#	@info(M1[1,1:2])
#	@info(M2[1,1:2])
	
	return (P - V1 + (V2.*(-M1*V2 + M2*V3) - V3.*(M1*V3 + M2*V2))) 
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
function RK_solver1(sp::SParams)
	dU = RK4(f1)
	oTdot = zeros(Float64, length(sp.T))
	Tdot = zeros(Float64, length(sp.T))
	
	M = sp.V'.*sp.Y.*sp.V
	
	M = M - spdiagm(diag(M))
	M1 = real(M)
	M2 = imag(M)
	# vector of the sum of M1 lines
	# NB: -sum(M1,2) returns an Array{Float64,1} in v0.45 and Array{Float64,2} in v0.5
	V1 = -sum(M1,2)[:,1]
	
	(dT, Tdot) = dU(sp.T, sp.o_args[:h], V1, M1, M2, sp.P)
	sp.T += dT
	n_iter = 2
	
	while n_iter < sp.iter_max
		oTdot = copy(Tdot)
		(dT, Tdot) = dU(sp.T, sp.o_args[:h], V1, M1, M2, sp.P)
		error1 = norm(Tdot-oTdot,Inf)
		error2 = var(Tdot)
		# the simulation has converged if all the theta derivatives have not changed between the last two iterations (error1 < epsilon) and all the theta derivatives are the same (error2 < epsilon)
		if error1 < sp.epsilon && error2 < sp.epsilon
			break
		end
		#@debug("# iter $n_iter with max velocity=$error1")
		#@debug("# iter $n_iter with max velocity change=$error2")
		sp.T += dT
		n_iter += 1
	end
	
	@info("$n_iter iterations")
	return State(Float64[],sp.T,Tdot,n_iter)
end


# Same as RK_solver1 but does not print the number of iteration

function RK_solver_quiet(sp::SParams)
	dU = RK4(f1)
	oTdot = zeros(Float64, length(sp.T))
	Tdot = zeros(Float64, length(sp.T))
	
	M = sp.V'.*sp.Y.*sp.V
	
	M = M - spdiagm(diag(M))
	M1 = real(M)
	M2 = imag(M)
	# vector of the sum of M1 lines
	# NB: -sum(M1,2) returns an Array{Float64,1} in v0.45 and Array{Float64,2} in v0.5
	V1 = -sum(M1,2)[:,1]
	
	(dT, Tdot) = dU(sp.T, sp.o_args[:h], V1, M1, M2, sp.P)
	sp.T += dT
	n_iter = 2
	
	while n_iter < sp.iter_max
		oTdot = copy(Tdot)
		(dT, Tdot) = dU(sp.T, sp.o_args[:h], V1, M1, M2, sp.P)
		error1 = norm(Tdot-oTdot,Inf)
		error2 = var(Tdot)
		# the simulation has converged if all the theta derivatives have not changed between the last two iterations (error1 < epsilon) and all the theta derivatives are the same (error2 < epsilon)
		if error1 < sp.epsilon && error2 < sp.epsilon
			break
		end
		#@debug("# iter $n_iter with max velocity=$error1")
		#@debug("# iter $n_iter with max velocity change=$error2")
		sp.T += dT
		n_iter += 1
	end
	
	return State(Float64[],sp.T,Tdot,n_iter)
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
# callback_func: callback function to be called at each iteration of the solver
#
## OUTPUT
# T: updated thetas
# Tdot: updated theta_dots
# n_iter: # of iterations before convergence
function RK_solver1(sp::SParams,callback_func::Function)
	dU = RK4(f1)
	oTdot = zeros(Float64, length(sp.T))
	Tdot = zeros(Float64, length(sp.T))
	
	M = sp.V'.*sp.Y.*sp.V
	
	M = M - spdiagm(diag(M))
	M1 = real(M)
	M2 = imag(M)
	# vector of the sum of M1 lines
	# NB: -sum(M1,2) returns an Array{Float64,1} in v0.45 and Array{Float64,2} in v0.5
	V1 = -sum(M1,2)[:,1]

	(dT, Tdot) = dU(sp.T, sp.o_args[:h], V1, M1, M2, sp.P)
	sp.T += dT
	n_iter = 2
	
	while n_iter < sp.iter_max
		oTdot = copy(Tdot)
		(dT, Tdot) = dU(sp.T, sp.o_args[:h], V1, M1, M2, sp.P)
		error1 = norm(Tdot-oTdot,Inf)
		error2 = var(Tdot)
		# the simulation has converged if all the theta derivatives have not changed between the last two iterations (error1 < epsilon) and all the theta derivatives are the same (error2 < epsilon)
		if error1 < sp.epsilon && error2 < sp.epsilon
			break
		end
		#@debug("# iter $n_iter with max velocity=$error1")
		#@debug("# iter $n_iter with max velocity change=$error2")
		sp.T += dT
		n_iter += 1
		go_on = callback_func(sp,n_iter,[error1,error2])
		!go_on && break
	end
	
	@info("$n_iter iterations")
	return State(Float64[],sp.T,Tdot,n_iter)
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
function SD_solver(sp::SParams)
	n_iter = 0
	n = length(sp.T)
	del = sp.o_args[:d]

	# We only use the susceptive part of the admittance matrix, with zero diagonal elements	
	K = imag(sp.Y-spdiagm(diag(sp.Y)))
	ST = sin(sp.T)
	CT = cos(sp.T)

	V1 = K*CT
	V2 = K*ST
		
	# f0 is the potential in the phase space whose extremas are solutions of the PFEs
	f0 = -sum(sp.P.*sp.T) - .5*sum(CT'*V1 + ST'*V2)

	# nabla is the gradient of this potential
	nabla = -sp.P + ST.*V1 - CT.*V2
	delta = norm(nabla,Inf)
	
	# Follow the path of most negative gradient, until the correction is less than delta, to reach the bottom of a well
	while delta > sp.epsilon && n_iter < sp.iter_max
		n_iter += 1
		A = copy(sp.T)
		sp.T -= nabla*del
		ST = sin(sp.T)
		CT = cos(sp.T)
		f1 = -sum(sp.P.*sp.T) - .5*sum(CT'*K*CT + ST'*K*ST)
		while f1 > f0 && norm(f1-f0) > sp.epsilon
			sp.T = copy(A)
			del = del/2
			sp.T -= nabla*del
			ST = sin(sp.T)
			CT = cos(sp.T)
			f1 = -sum(sp.P.*sp.T) - .5*sum(CT'*K*CT + ST'*K*ST)
		end
		f0 = copy(f1)
		nabla = -sp.P + ST.*(K*CT) - CT.*(K*ST)
		delta = norm(nabla,Inf)
		@debug("# iter $n_iter with error=$delta")
	end

	if n_iter == sp.iter_max
		sp.T = zeros(n)
		@info(" Did not converge.")
	end

	o_data = Dict{Symbol,Any}()
	o_data[:delta] = delta
	return State(Float64[],sp.T,Float64[],n_iter,o_data)
end

