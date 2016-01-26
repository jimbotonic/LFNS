# compute the vorticity of the graph
function vorticity(T::Array{Float64,1}, cycles::Array{Array{Int64,1},1})
	X = Float64[]
	for c in cycles
		c = copy(cycles[i])
		push!(c,c[1])
		s = 0.
		for j in 1:(length(c)-1)
			s += mod(c[j] - c[j+1] + pi, 2*pi) - pi
		end
		push!(X,abs(s))
	end
	return X
end

# computes the winding number along the given sequence of nodes
# it is assumed that the given sites form a cycle
#
## INPUT
# T: thetas
#
## OUPUT
# q: winding number
#
function winding_number(T::Array{Float64,1})
	# dT[i] = T[i-1]-T[i]
	dT = [T[end]-T[1],]
	n = length(T)
	
	for i in 2:n
		push!(dT,T[i-1]-T[i])
	end
	
	dT = mod(dT+pi,2*pi)-pi
	q = sum(dT)/(2*pi)
	
	return q
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
function get_stability_matrix(T::Array{Float64,1}, Y::Array{Complex{Float64},2}, approx_level::Int64=3)
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

# computes the flows at each node and on every line
#
## INPUT
# T: thetas
# Y: admittance matrix
# 
## OUTPUT
# P: vector of power balance at each node
# F: matrix of power flows on every line
# L: matrix of losses on every line
function flow_test(T::Array{Float64,1}, Y::Array{Complex{Float64},2})
	n = length(T)
	
	YY = Y.*(1-eye(n))
	YY = YY-diagm(collect(sum(YY,2)))
	
	G = real(YY)
	B = imag(YY)
	
	dT = T*ones(1,n)-ones(n,1)*T'
	F = B.*sin(dT)+G.*cos(dT)
	P = collect(sum(F.*(1-eye(n)),2))
	L = F+F'
	
	return P,F,L
end
