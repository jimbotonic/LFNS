# compute the vorticity of the graph
function vorticity(T::Array{Float64,1}, cycles::Array{Array{Int64,1},1})
	V = Float64[]
	for cy in cycles
		c = copy(cy)
		push!(c,c[1])
		s = 0.
		for j in 1:(length(c)-1)
			s += mod(T[c[j]] - T[c[j+1]] + pi, 2*pi) - pi
		end
		# entries are multiples of 2*pi
		push!(V,round(Int,abs(s/(2*pi))))
	end
	return V
end

# compute the vorticity of a given cycle
function vorticity(T::Array{Float64,1}, cycle::Array{Int64,1})
	c = copy(cycle)
	push!(c,c[1])
	s = 0.
	for j in 1:(length(c)-1)
		s += mod(T[c[j]] - T[c[j+1]] + pi, 2*pi) - pi
	end
	return round(Int,abs(s/(2*pi)))
end

# find the positions of the vortices
#
# n: height of the lattice
# m: width of the lattice
function find_vortices_in_sq_lattice(n::Int,m::Int,T::Array{Float64,1})
	X = Int[]
	Y = Int[]
	V = Float64[]
	for i in 1:(n-1)
		for j in 1:(m-1)
			cycle = Int[]
			push!(cycle, (i-1)*n+j)
			push!(cycle, (i-1)*n+j+1)
			push!(cycle, i*n+j+1)
			push!(cycle, i*n+j)
			v = vorticity(T,cycle)
			if v > 0
				push!(X,j)	
				push!(Y,i)	
				push!(V,v)	
			end
		end	
	end
	return X,Y,V
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
function get_stability_matrix(T::Array{Float64,1}, Y::SparseMatrixCSC{Complex{Float64},Int64}, approx_level::Int64=3)
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
function get_lambda2(T::Array{Float64,1}, Y::SparseMatrixCSC{Complex{Float64},Int64}, epsilon::Float64=1e-12)
	# get stability matrix
	M = get_stability_matrix(T,Y)
	n = size(M)[1]
	M = Symmetric(M)
	evs = eigvals(M,(n-1:n))
	l1 = evs[end]
	
	if abs(l1) > epsilon
		return l1
	else
		return evs[end-1]
	end
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
function flow_test(T::Array{Float64,1}, Y::SparseMatrixCSC{Complex{Float64},Int64})
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

