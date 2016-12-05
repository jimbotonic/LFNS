# compute distance between 2 vectors of angles
#
# norm_type: typically 2 or Inf
function distance(T_ref::Array{Float64,1}, T::Array{Float64,1}, norm_type::Real)
	return norm(mod(T_ref[i] - T[i] + pi, 2pi) - pi, norm_type)
end

# estimate slope of function f at x numerically
#
# NB: f is assumed to be continuous and differentiable 
function slope(f::Function, x::Float64)
	# compute smallest possible h based on machine precision
	h = sqrt(eps(typeof(x)))
	return (f(x+h)-f(x))/h
end

# compute the vorticity of the graph
function vorticity(T::Array{Float64,1}, cycles::Array{Array{Int64,1},1})
	V = Float64[]
	for c in cycles
		s = 0.
		for j in 1:(length(c)-1)
			s += mod(T[c[j]] - T[c[j+1]] + pi, 2pi) - pi
		end
		s += mod(T[c[end]] - T[c[1]] + pi, 2pi) - pi
		# entries are multiples of 2*pi
		push!(V,round(Int,s/(2pi)))
	end
	return V
end

# compute the vorticity of a given cycle
function vorticity(T::Array{Float64,1}, cycle::Array{Int64,1})
	s = 0.
	for j in 1:(length(cycle)-1)
		s += mod(T[cycle[j]] - T[cycle[j+1]] + pi, 2pi) - pi
	end
	s += mod(T[cycle[end]] - T[cycle[1]] + pi, 2pi) - pi
	return round(Int,s/(2pi))
end

# compute the vorticity of a given cycle
function vorticity2(T::Array{Float64,1}, cycle::Array{Int64,1})
	s = 0.
	for j in 1:(length(cycle)-1)
		# d in [-2pi,2pi]
		d = mod(T[cycle[j]], 2pi) - mod(T[cycle[j+1]], 2pi)
		if d > pi 
			d -= 2pi
		elseif d < -pi
			d += 2pi
		end
		# d in [-pi,pi]
		s += d
	end
	d = mod(T[cycle[end]], 2pi) - mod(T[cycle[1]], 2pi)
	if d > pi 
		d -= 2pi
	elseif d < -pi
		d += 2pi
	end
	# d in [-pi,pi]
	s += d
	return round(Int,s/(2pi))
end

# find the positions of the vortices
#
# matrix-like coordinates:
# -> n: height of the lattice
# -> m: width of the lattice
function find_vortices_in_sq_lattice(n::Int,m::Int,T::Array{Float64,1})
	X = Int[]
	Y = Int[]
	V = Float64[]
	for i in 1:(n-1)
		for j in 1:(m-1)
			cycle = Int[]
			push!(cycle, (i-1)*m+j)
			push!(cycle, (i-1)*m+j+1)
			push!(cycle, i*m+j+1)
			push!(cycle, i*m+j)
			v = vorticity(T,cycle)
			if v > 0
				push!(X,i)	
				push!(Y,j)	
				push!(V,v)	
			end
		end	
	end
	return X,Y,V
end

# get the value of the potential
#
# matrix-like coordinates:
# -> n: height of the lattice
# -> m: width of the lattice
function get_potential_in_sq_lattice(n::Int,m::Int,P::Array{Float64,1},T::Array{Float64,1})
	u = 0.
	# compute \sum P_i * theta_i
	for i in 1:n
		for j in 1:m
			pos = (i-1)*m+j
			v -= P[pos]*T[pos]
		end
	end
	# 4 corners
	pos = 1
	v -= cos(T[pos+1]-T[pos])
	v -= cos(T[pos+m]-T[pos])
	pos = m
	v -= cos(T[pos-1]-T[pos])
	v -= cos(T[pos+m]-T[pos])
	pos = (n-1)*m + 1 
	v -= cos(T[pos+1]-T[pos])
	v -= cos(T[pos-m]-T[pos])
	pos = n*m
	v -= cos(T[pos-1]-T[pos])
	v -= cos(T[pos-m]-T[pos])
	# first line
	for i in 1:1
		for j in 2:(m-1)
			pos = (i-1)*m+j
			v -= cos(T[pos-1]-T[pos])
			v -= cos(T[pos+1]-T[pos])
			v -= cos(T[pos+m]-T[pos])
		end
	end
	# last line
	for i in n:n
		for j in 2:(m-1)
			pos = (i-1)*m+j
			v -= cos(T[pos-1]-T[pos])
			v -= cos(T[pos+1]-T[pos])
			v -= cos(T[pos-m]-T[pos])
		end
	end
	# left column
	for i in 2:(n-1)
		for j in 1:1
			pos = (i-1)*m+j
			v -= cos(T[pos-m]-T[pos])
			v -= cos(T[pos+m]-T[pos])
			v -= cos(T[pos+1]-T[pos])
		end
	end
	# right column
	for i in 2:(n-1)
		for j in m:m
			pos = (i-1)*m+j
			v -= cos(T[pos-m]-T[pos])
			v -= cos(T[pos+m]-T[pos])
			v -= cos(T[pos-1]-T[pos])
		end
	end
	# inside of the lattice
	for i in 2:(n-1)
		for j in 2:(m-1)
			pos = (i-1)*m+j
			v -= cos(T[pos-1]-T[pos])
			v -= cos(T[pos+1]-T[pos])
			v -= cos(T[pos-m]-T[pos])
			v -= cos(T[pos+m]-T[pos])
		end
	end

	return v
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

function get_stability_matrix(g::Graphs.AbstractGraph{Bus,Line}, approx_level::Int64=3)
	Y = get_admittance_matrix(g)
	T = get_angles(g)
	
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

# get the eigenvalues
#
## INPUT
# T: thetas
# Y: admittance matrix
#
## OUTPUT
# eigv: array of the eigenvalues
function get_eigenvalues(T::Array{Float64,1}, Y::SparseMatrixCSC{Complex{Float64},Int64})
	# get stability matrix
	M = get_stability_matrix(T,Y)
	n = size(M)[1]
	eigv = eigvals(M)
	return eigv
end	

function get_eigenvalues(g::Graphs.AbstractGraph{Bus,Line})
	# get stability matrix
	M = get_stability_matrix(g)
	n = size(M)[1]
	eigv = eigvals(M)
	return eigv
end	


# get the real part of lambda 2
#
## INPUT
# T: thetas
# Y: admittance matrix
#
## OUTPUT
# l2: greatest non-null real part of the stability matrix eigenvalues
function get_lambda2(T::Array{Float64,1}, Y::SparseMatrixCSC{Complex{Float64},Int64}, epsilon::Float64=1e-12)
	evs = sort(real(get_eigenvalues(T,Y)))
	l1 = evs[end]

	if abs(l1) > epsilon
		return l1
	else
		return evs[end-1]
	end
end

function get_lambda2(g::Graphs.AbstractGraph{Bus,Line}, epsilon::Float64=1e-12)
	evs = sort(real(get_eigenvalues(g)))
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
# Y: admittance matrix with positive imaginary part and negative real part for out-diagonal elements !
# 
## OUTPUT
# P: vector of power balance at each node
# F: matrix of power flows on every line
# L: matrix of losses on every line
function flow_test(T::Array{Float64,1}, Y::SparseMatrixCSC{Complex{Float64},Int64})
	n = length(T)
	YY = Y - spdiagm(diag(Y))
	G = -real(YY)
	B = imag(YY)
	dT = T*ones(1,n) - ones(n,1)*T'

	F = B.*sin(dT) + G.*(1 - cos(dT))
	P = collect(sum(F,2))
	L = F + F'
	
	return P,F,L
end

