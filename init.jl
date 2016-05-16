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

# initialize P vector with entry values in [-1, 1] and absolute value of greatest value equal to 1
#
# switch sign of entries of odd position
# NB: U is assumed to be a distribution
function init_P4(U::Array{Float64,1})
	P = copy(U)
	for i in 1:length(P)
		if isodd(i)
			P[i] = -P[i]
		end
	end
	# make sure that sum(P)=0
	P -= mean(P)
	P *= 1/maximum(abs(P))
	return P
end

# generate a ring with one producer at vertex 1 and one consumer at a chosen vertex
#
## INPUT
# n: length of the cycle
# ic: index of the vertex of the consumer
# p: produced/consumed power
function generate_ring(n::Int,ic::Int,p::Float64)
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
# n: width of the lattice
# m: height of the lattice
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

# initialize a vector T to create vortex on square lattice
# create a vortex on a square at the specified position
# initialize all other angles by using atan
#
# INPUT
# matrix-like coordinates:
# -> n,m: (height,width) of the lattice
# -> i,j: (row,column) coordinates of the vortex center
function create_vortex_on_sq_lattice(n::Int,m::Int,i::Int,j::Int)
	T = zeros(Float64,n*m)
	for p in 1:(n*m)
		# get coordinates of the current point 
		x = mod(p,m)
		if x == 0
			x = m
		end
		y = ceil(Int,p/m)
		# center of the vortex
		cx = j + 1/2
		cy = i + 1/2
		ba = atan(abs(y-cy)/abs(x-cx))
		# NO
		# NB: each square has a 1 unit length
		if x <= cx && y <= cy
			T[p] = -ba + pi
		# NE
		elseif x >= cx && y <= cy
			T[p] = ba 
		# SE
		elseif x >= cx && y >= cy
			T[p] = -ba
		# SO
		elseif x <= cx && y >= cy
			T[p] = ba - pi
		end
	end	

	return T
end

# initialize a vector T to create vortex on square lattice
# create a vortex on a square at the specified position
# initialize all other angles by using atan
#
# INPUT
# matrix-like coordinates:
# -> n,m: (height,width) of the lattice
# -> i,j: (row,column) coordinates of the vortex center
function create_antivortex_on_sq_lattice(n::Int,m::Int,i::Int,j::Int)
	T = zeros(Float64,n*m)
	for p in 1:(n*m)
		# get coordinates of the current point 
		x = mod(p,m)
		if x == 0
			x = m
		end
		y = ceil(Int,p/m)
		# center of the vortex
		cx = j + 1/2
		cy = i + 1/2
		ba = atan(abs(y-cy)/abs(x-cx))
		# NO
		# NB: each square has a 1 unit length
		if x <= cx && y <= cy
			T[p] = ba - pi
		# NE
		elseif x >= cx && y <= cy
			T[p] = -ba
		# SE
		elseif x >= cx && y >= cy
			T[p] = ba
		# SO
		elseif x <= cx && y >= cy
			T[p] = -ba + pi
		end
	end	

	return T
end

# generate the contour cycle of a square lattice
#
# n: width of the lattice
# m: height of the lattice
function get_sq_lattice_contour_cycle(n::Int,m::Int)
	bcycle = Array{Int64,1}()

	for j in 1:n
		push!(bcycle,j)
	end	
	for i in 2:m
		push!(bcycle,i*n)
	end	
	for j in 1:(n-1)
		push!(bcycle,m*n-j)
	end	
	for i in 2:m
		push!(bcycle,(m-i)*n+1)
	end	
	return bcycle
end

# generate the contour cycle of a graph from (lat,lng) pairs
#
# mode A: start from leftest node and visit nodes clockwise
# -> angles are measured for y=0+ axis
# -> from leftest node, choose next child so that angle is maximum in upper-right quadrant
# -> for next nodes, choose next child so that:
# 	-> angle is maximum if first crossed edge has been visited
#	-> angle is maximum before an already visited edge is crossed
#
# mode B: start from leftest node and visit nodes anticlockwise
# -> angles are measured with respect to y=0- axis
# -> from leftest node, choose next child so that angle is maximum in lower-right quadrant
# -> for next nodes, choose next child so that:
# 	-> angle is maximum if first crossed edge has been visited
#	-> angle is maximum before an already visited edge is crossed
function get_graph_contour_cycle(g::Graphs.AbstractGraph{Bus,Line})
	bcycle = Array{Int64,1}()

	X = Float64[]
	Y = Float64[]
	X2 = Float64[]
	Y2 = Float64[]
	ids = Set{Int}()

	for v in vertices(g)
		push!(X,v.lat)
		push!(Y,v.lng)
		# ignore sources and sinks
		if out_degree(v,g) > 1
			push!(ids,v)
			push!(X2,v.lat)
			push!(Y2,v.lng)
		end
	end

	# get leftest vertex
	lvi = indmin(X2)
	push!(bcycle,lvi)

	cv = lvi
	while true
		nei = out_neighbors(cv,g)
		
	end
			
	return bcycle
end

# generate a flat square lattice on a cylinder
#
## INPUT
# n: width of the lattice
# m: height of the lattice
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
# n: width of the lattice
# m: height of the lattice
function generate_sq_lattice_on_torus(n::Int,m::Int)
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
	for i in (4n^2+1):(4n^2+2*n^2)
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

	# stitch face 5 and 6 with the 4 first faces
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

	# stitch face 1 (up) with face 4 (bottom)
	for i in 1:n
		push!(es, Line(ecounter,vs[i],vs[4n^2-n+i],-1.im))
		ecounter += 1
	end
	
	return graph(vs, es, is_directed=false)
end

# generate a flat triangular lattice of equilateral shape 
# with corners at positions (west, north-east, south-east)
# numbering of the sites is from bottom to top and from left to right
#
## INPUT
# n: length of the side of the lattice
function generate_tri_lattice(n::Int)
	vs = Bus[]
	es = Line[]	
	vcounter = 0
	ecounter = 0
	
	sites = zeros(Int(n*(n+1)/2),3)
	
	for k in 0:(n-1)
		for l in 0:k
			vcounter += 1
			sites[vcounter,:] = [vcounter, l, k-l]
			push!(vs,Bus(vcounter,0.,0.))
		end
	end
	
	for i in 1:Int((n*(n+1)/2))
		if sites[i,3] > 0
			ecounter += 1
			push!(es, Line(ecounter,vs[i],vs[i+1],-1.im))
			#push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.im))
		end
		if (sites[i,2]+sites[i,3]) < n-1
			ecounter += 1
			push!(es, Line(ecounter,vs[i],vs[Int(sum(sites[i,1:3]))+1],-1.im))
			ecounter += 1
			push!(es, Line(ecounter,vs[i],vs[Int(sum(sites[i,1:3]))+2],-1.im))
		end
	end
		
	return graph(vs, es, is_directed=false)
end

# initialize a vector T to create cortex on triangular lattice
# 
# INPUT
# n: length of the side of the lattice
# i,j: (column number, triangle number in the column from bottom) coordinates of the triangle carrying the vortex
# 1 <= i <= n-1, 1 <= j <= 2*i-1
function create_vortex_on_tri_lattice(n::Int,i::Int,j::Int)
	T = zeros(Float64,Int(n*(n+1)/2))
	
	xcoords = Float64[]
	ycoords = Float64[]

	# get the coordinates of the nodes in a orthonormal coordinate system	
	for k in 0:(n-1)
		for l in 0:k
			push!(xcoords,k*cos(pi/6))
			push!(ycoords,(2*l-k)*sin(pi/6))
		end
	end

	# get the index of one of the nodes of the triangle carrying the vortex
	idx = Int(i*(i+1)/2) + 1
	
	if j%2 == 1
		idx += Int((j-1)/2)
		xvortex = xcoords[idx] + sqrt(3)/3
		yvortex = ycoords[idx]
	else
		idx += Int((j-2)/2)
		xvortex = xcoords[idx] + sqrt(3)/6
		yvortex = ycoords[idx] + .5
	end	

	# compute angles
	dx = xcoords - xvortex
	dy = ycoords - yvortex
	
	t = (dx.<0)
	T = atan(dy./dx) + t*pi
	T = mod(T+pi,2*pi) - pi
	
	return T
end

# generate the contour cycle of a triangular lattice
function get_tri_lattice_contour_cycle(n::Int)
	bcycle = [1,]
	elcycb = [1,]
	
	x = 1
	
	for i in 2:n
		x += i
		unshift!(elcycb,x)
		push!(bcycle,elcycb[2]+1)
	end
	
	append!(bcycle,elcycb)
	
	return bcycle
	
end

# generate a double ring with a bus where p is injected
# producer and consumer are located at the degree 3 vertices
#
## INPUT
# l,c,r: number of VERTICES on each branch of the double cycle
# p: produced/consumed power
function generate_double_ring(l::Int,c::Int,r::Int,p::Float64)
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
