using Distributions

include("graphs.jl")

# initialize a random distribution whose entries are drawn uniformaly at random
function init_unif_dist(n::Int64)
	return Float64[1/n for i in 1:n]
end

# initialize a random distribution whose entries are drawn uniformaly at random
function init_rand_dist(n::Int64)
	U = rand(Uniform(),n)
	return U/sum(U)
end

# initialize P vector with entry values in [-1, 1] and absolute value of greatest value equal to 1
#
# NB: U is assumed to be a distribution
function init_P1(U::Array{Float64,1})
	P = U*2-1
	# make sure that sum(P)=0
	P -= mean(P)
	P *= 1/maximum(abs(P))
	return P
end

# initialize P vector with entry values in [-1, 1] and absolute value of greatest value equal to 1
#
# switch entries sign at random
# NB: U is assumed to be a distribution
function init_P2(U::Array{Float64,1})
	n = length(U)
	# switch entry signs randomly
	P = rand([-1,1],n).*U
	# make sure that sum(P)=0
	P -= mean(P)
	P *= 1/maximum(abs(P))
	return P
end

# initialize P vector with entry values in [-1, 1] and absolute value of greatest value equal to 1
#
# switch entries sign such that the highest value is positive, the 2nd highest value is negative, the 3rd is positive, ...
# NB: U is assumed to be a distribution
function init_P3(U::Array{Float64,1})
	n = length(U)
	V = sortrows([U 1:n])
	V[:,1] = V[:,1].*(2*mod(1:n,2)-1)
	V = [V[:,2] V[:,1]]
	V = sortrows(V)
	P = V[:,2]
	# make sure that sum(P)=0
	P -= mean(P)
	P *= 1/maximum(abs(P))
	return P
end

# generate a cycle with one producer at vertex 1 and one consumer at a chosen vertex
#
## INPUT
# n: length of the cycle
# ic: index of the vertex of the consumer
# p: produced/consumed power
function generate_cycle(n::Int,ic::Int,p::Float64)
	vs = Bus[]
	es = Line[]
	
	push!(vs,Bus(1,0.,p))
	for i in 2:n
		if i==ic
			push!(vs,Bus(i,0.,-p))
		else
			push!(vs,Bus(i,0.,0.))
		end
		push!(es,Line(i-1,vs[i-1],vs[i],-1.im))
	end
	push!(es,Line(n,vs[n],vs[1],-1.im))
	
	return graph(vs, es, is_directed=false)
end

# generate a flat square lattice
#
## INPUT
# n,m: width and height of the lattice
function generate_sq_lattice(n::Int,m::Int)
	vs = Bus[]
	es = Line[]	
	ecounter = 1

	for i in 1:n*m
		push!(vs,Bus(i,0.,0.))
	end

	# rows 
	for i in 0:(m-1)
		for j in 1:(n-1)
			push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(m-2)
			push!(es, Line(ecounter,vs[j*n+i],vs[(j+1)*n+i],-1.im))
			ecounter += 1
		end
	end
	
	return graph(vs, es, is_directed=false)
end

# generate a flat square lattice on a cylinder
#
## INPUT
# n,m: width and height of the lattice
function generate_sq_lattice_on_cylinder(n::Int,m::Int)
	vs = Bus[]
	es = Line[]	
	ecounter = 1

	for i in 1:n*m
		push!(vs,Bus(i,0.,0.))
	end

	# rows 
	for i in 0:(m-1)
		for j in 1:(n-1)
			push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(m-2)
			push!(es, Line(ecounter,vs[j*n+i],vs[(j+1)*n+i],-1.im))
			ecounter += 1
		end
	end
	
	# stitch left and right sides
	for i in 0:(m-1)
		push!(es, Line(ecounter,vs[i*n+1],vs[n*(i+1)],-1.im))
		ecounter += 1
	end

	return graph(vs, es, is_directed=false)
end

# generate a flat square lattice on a donut
#
## INPUT
# n,m: width and height of the lattice
function generate_sq_lattice_on_donut(n::Int,m::Int)
	vs = Bus[]
	es = Line[]	
	ecounter = 1

	for i in 1:n*m
		push!(vs,Bus(i,0.,0.))
	end

	# rows 
	for i in 0:(m-1)
		for j in 1:(n-1)
			push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(m-2)
			push!(es, Line(ecounter,vs[j*n+i],vs[(j+1)*n+i],-1.im))
			ecounter += 1
		end
	end
	
	# stitch left and right sides
	for i in 0:(m-1)
		push!(es, Line(ecounter,vs[i*n+1],vs[n*(i+1)],-1.im))
		ecounter += 1
	end
	
	# stitch up and bottom sides
	for i in 1:n
		push!(es, Line(ecounter,vs[i],vs[m*(n-1)+i],-1.im))
		ecounter += 1
	end

	return graph(vs, es, is_directed=false)
end

# generate a square lattice on the sphere
#
## INPUT
# n: size of the cubic faces
function generate_sq_lattice_on_sphere(n::Int)
	vs = Bus[]
	es = Line[]	
	ecounter = 1

	# add 4 faces
	for i in 1:4n^2
		push!(vs,Bus(i,0.,0.))
	end

	# connect the 4 faces
	#
	# rows 
	for i in 0:(4n-1)
		for j in 1:(n-1)
			push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(4n-2)
			push!(es, Line(ecounter,vs[j*n+i],vs[(j+1)*n+i],-1.im))
			ecounter += 1
		end
	end
	
	# add 2 more faces
	for i in 1:2*n^2
		push!(vs,Bus(i,0.,0.))
	end

	# connect the 2 faces
	#
	# rows 
	for i in 0:(n-1)
		for j in 1:(n-1)
			push!(es, Line(ecounter,vs[4n^2+i*n+j],vs[4n^2+i*n+j+1],-1.im))
			ecounter += 1
			push!(es, Line(ecounter,vs[5n^2+i*n+j],vs[5n^2+i*n+j+1],-1.im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(n-2)
			push!(es, Line(ecounter,vs[4n^2+j*n+i],vs[4n^2+(j+1)*n+i],-1.im))
			ecounter += 1
			push!(es, Line(ecounter,vs[5n^2+j*n+i],vs[5n^2+(j+1)*n+i],-1.im))
			ecounter += 1
		end
	end

	# stitch face 5 and 6 with the 4 first
	for i in 1:n
		# up
		push!(es, Line(ecounter,vs[4n^2+i],vs[2n^2+(i-1)*n+1],-1.im))
		ecounter += 1
		push!(es, Line(ecounter,vs[5n^2+i],vs[3n^2-(i-1)*n],-1.im))
		ecounter += 1
		# left
		push!(es, Line(ecounter,vs[4n^2+(i-1)*n+1],vs[2n^2-i*n+1],-1.im))
		ecounter += 1
		push!(es, Line(ecounter,vs[5n^2+(i-1)*n+1],vs[3n^2+i*n],-1.im))
		ecounter += 1
		# right
		push!(es, Line(ecounter,vs[4n^2+i*n],vs[3n^2+(i-1)*n+1],-1.im))
		ecounter += 1
		push!(es, Line(ecounter,vs[5n^2+i*n],vs[2n^2-(i-1)*n],-1.im))
		ecounter += 1
		# bottom
		push!(es, Line(ecounter,vs[5n^2-n+i],vs[n^2-i*n+1],-1.im))
		ecounter += 1
		push!(es, Line(ecounter,vs[6n^2-n+i],vs[n^2-(i-1)*n],-1.im))
		ecounter += 1
	end
	
	return graph(vs, es, is_directed=false)
end

# generate a double cycle with a bus where p is injected
# producer and consumer are located at the degree 3 vertices
#
## INPUT
# l,c,r: number of VERTICES on each branch of the double cycle
# p: produced/consumed power
function generate_double_cycle(l::Int,c::Int,r::Int,p::Float64)
	vs = Bus[]
	es = Line[]	
	ecounter = 1

	# left branch 1:l
	for i in 1:l
		# t=0., p=0.
		push!(vs,Bus(i,0.,0.))
	end
	for i in 2:l
		push!(es, Line(ecounter,vs[i-1],vs[i],-1.im))
		ecounter += 1
	end
	# central branch (l+1):(c+l)
	for i in (l+1):(c+l)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+2):(c+l)
		push!(es, Line(ecounter,vs[i-1],vs[i],-1.im))
		ecounter += 1
	end
	# right branch (l+c+1):(c+l+r)
	for i in (l+c+1):(c+l+r)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+c+2):(c+l+r)
		push!(es, Line(ecounter,vs[i-1],vs[i],-1.im))
		ecounter += 1
	end

	# 2 last remaining vertices
	push!(vs,Bus(l+c+r+1,0.,p))
	push!(vs, Bus(l+c+r+2,0.,-p))

	push!(es, Line(ecounter,vs[l+c+r+1],vs[1],-1.im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+1],vs[l+1],-1.im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+1],vs[l+c+1],-1.im))
	ecounter += 1

	push!(es, Line(ecounter,vs[l+c+r+2],vs[l],-1.im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+2],vs[l+c],-1.im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+2],vs[l+c+r],-1.im))

	return graph(vs, es, is_directed=false)
end
