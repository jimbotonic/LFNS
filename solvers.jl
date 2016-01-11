using Distances

include("simulator.jl")

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
    	#n = size(sp.Y)[1]
	
	# node ids whose bus type is 0
	PQ_ids = Int64[v.id for v in filter(v -> v.bus_type == 0, vs)]
	# set PQ bus voltages to 1 pu
	sp.V[PQ_ids] = 1.
	# node id whose bus type is 3
	slack_id = filter(v -> v.bus_type == 3, vs)[1].id

	# set of node ids which are part of the connected component containing the slack bus 
	sc_ids = get_slack_component_ids(sp.g)
	sp.Y = sp.Y[sc_ids, sc_ids]
	sp.V = sp.V[sc_ids]
	sp.T = sp.T[sc_ids]
	sp.P = sp.P[sc_ids]
	sp.Q = sp.Q[sc_ids]

	# bus types of the slack component 	
	bus_type = Int64[v.bus_type for v in vs]
	# PQ bus positions in the slack component
	PQ_pos = findin(bus_type[sc_ids], 0)
	# slack position in the slack component
	slack_pos = findin(bus_type[sc_ids],3)[1]
	
	# bootstrap simulation
	sp.o_args[:PQ_pos] = PQ_pos
	sp.o_args[:slack_pos] = slack_pos
	GS_solver(sp)

	# compute Y element-wise absolute values
    	Y_abs = abs(sp.Y)
	# compute Y element-wise angle values
    	Y_angle = angle(sp.Y)
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
		error = maximum(abs(dPQ))
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
               ( (k4   )-> (( k1 + 2*k2 + 2*k3 + k4 )/6, k1) 
               )( h * f( y + k3  , v... ) )
               )( h * f( y + k2/2, v... ) )
               )( h * f( y + k1/2, v... ) )
               )( h * f( y       , v... ) )
end

# right-hand side of the differential equation 
function f1(T::Array{Float64,1}, V::Array{Float64,1}, Y::Array{Complex{Float64},2}, P::Array{Float64,1})
	M1 = V*V'
	M2 = real(Y).*M1 #M2ij=Gij*Vi*Vj
	M3 = imag(Y).*M1 #M3ij=Bij*Vi*Vj
	V1 = diag(M2)
	V2 = cos(T)
	V3 = sin(T)
	# set diagonal elements to 0
	M2 = M2 - diagm(diag(M2))
	M3 = M3 - diagm(diag(M3))

	return (P - V1 + (V2.*(-M2*V2 + M3*V3) - V3.*(M2*V3 + M3*V2))) 
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
	nTdot = zeros(Float64, length(sp.T))
	Tdot = zeros(Float64, length(sp.T))

	n_iter = 1
	while n_iter < sp.iter_max
		(dT, Tdot) = dU(sp.T, sp.o_args["h"], sp.V, sp.Y, sp.P)
		if chebyshev(nTdot, Tdot) < sp.epsilon
			break
		end
		nTdot = copy(Tdot)
		sp.T += dT
		n_iter += 1
	end

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
	K = imag(sp.Y-diagm(diag(sp.Y)))
	dT = sp.T*ones(1,n)-ones(n,1)*sp.T'
	
	# f0 is the potential in the phase space whose extremas are solutions of the PFEs
	f0 = -sum(sp.P.*sp.T) - .5*sum(K.*cos(dT))
	
	# nabla is the gradient of this potential
	nabla = -sp.P + sum(K.*sin(dT),2)
	delta = norm(nabla)
	
	# Follow the path of most negative gradient, until the correction is less than delta, to reach the bottom of a well
	while delta > sp.epsilon && n_iter < sp.iter_max
		n_iter += 1
		A = copy(sp.T)
		T -= nabla*del
		f1 = -sum(sp.P.*sp.T) - .5*sum(K.*cos(sp.T*ones(1,n)-ones(n,1)*sp.T'))
		while f1 > f0 && norm(f1-f0) > sp.epsilon
			sp.T = copy(A)
			del = del/2
			sp.T -= nabla*del
			f1 = -sum(sp.P.*sp.T) - .5*sum(K.*cos(sp.T*ones(1,n)-ones(n,1)*sp.T'))
		end
		del = sp.o_args["d"]
		dT = sp.T*ones(1,n)-ones(n,1)*sp.T'
		f0 = copy(f1)
		nabla = -sp.P + sum(K.*sin(dT),2)
		delta = norm(nabla)
	end
	
	# rotate all angles by using the last one as the reference
	# all angles belong to [-pi,pi] afterward
	sp.T = mod(sp.T-sp.T[end]+pi,2*pi)-pi

	o_data = Dict{Symbol,Any}()
	o_data[:delta] = delta
	return State(Float64[],sp.T,Float64[],n_iter,o_data)
end

# Stability matrix
#
## INPUT
# T: thetas
# Y: admittance matrix
# approximation_level: defines which approximation is used for the stability
#	1: linearized equations, i.e. the sines are linearized and G=0.
#	2: lossless case, i.e. G=0
#	3: lossy case without reactive power
# TODO	>=4: full case with reactive power
#
## OUTPUT
# M: stability matrix
#
function get_stability_matrix(T::Array{Float64,1}, Y::Array{Complex{Float64},2}, approx_level::Int64)
	B = imag(Y) 
	G = real(Y)
	
	n = length(T)
	dT = T*ones(1,n)-ones(n,1)*T'
	
	# 1: linearized equations, i.e. the sines are linearized and G=0.
	if approx_level == 1
		M = B.*dT
	# 2: lossless case, i.e. G=0
	elseif approx_level == 2
		M = B.*cos(dT)
	# 3: lossy case without reactive power
	elseif approx_level == 3
		M = B.*cos(dT) + G.*sin(dT)
	# 4: full case with reactive power
	else
		# TODO
	end
	
	M = M-diagm(collect(sum(M,1)))
	
	return M
end

# get lambda 2
#
## INPUT
# T: thetas
# Y: admittance matrix
#
## OUTPUT
# l2: greatest non-null real part of the stability matrix eigenvalues
function get_lambda2(T::Array{Float64,1}, Y::Array{Complex{Float64},2}, epsilon::Float64=1e-12)
	# get stability matrix
	M = get_stability_matrix(T,Y)
	evs = eigvals(M)
	l1 = evs[end]
	
	if abs(l1) > epsilon
		l2 = l1
	else
		l2 = evs[end-1]
	end
	
	return l2
end
