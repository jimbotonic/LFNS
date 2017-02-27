using Distributions

include("graphs.jl")

# initialize a uniform distribution whose entries are all equal
function init_unif_dist(n::Int64)
	return Float64[1/n for i in 1:n]
end

# initialize a random distribution whose entries are drawn uniformly at random
function init_rand_dist(n::Int64)
	U = rand(Uniform(),n)
	return U/sum(U)
end

# generate random unit-length (norm-2) n-dimensional vector with positive entries
function get_rand_unit_vector_N2(n::Int)
	R = rand(n)
	return R./sqrt(sumabs2(R,1))
end

# generate random unit-length (norm-Inf) n-dimensional vector with positive entries
function get_rand_unit_vector_Ninf(n::Int)
	R = rand(n)
	return R./maximum(R)
end

# generate random with entries drawn uniformly around 0 in [-hwidth,hwidth]
function get_rand_unif_vector(n::Int, hwidth::Number=2.)
	T = Float64[]
	for i in 1:n
		push!(T,2*hwidth*rand(Uniform())-hwidth)
	end
	return T
end

# perturbe T vector (T_ref += alpha*(2pi*rand()-pi)
# the vector entries are distributed uniformly in [-pi,pi]
function perturbe_T_Unif(T_ref::Array{Float64,1}, alpha::Float64=1.)
	T = copy(T_ref)
	n = length(T)
	for i in 1:n
		T[i] += alpha*(2pi*rand(Uniform())-pi)
	end
	return T
end

# perturbe T vector (T_ref += alpha*(2pi*R2-pi))
function perturbe_T_N2(T_ref::Array{Float64,1}, alpha::Float64=1.)
	T = copy(T_ref)
	n = length(T)
	# get random vector of norm-2 1 (i.e. R belongs to the hypercube of side 1)
	R = get_rand_unit_vector_N2(n)
	for i in 1:n
		# the perturbative vector belongs to the hypersphere of radius pi*alpha centered at the origin
		T[i] += (R[i]*2pi-pi)*alpha
	end
	return T
end

# perturbe T vector (T_ref += alpha*(2pi*Rinf-pi))
function perturbe_T_Ninf(T_ref::Array{Float64,1}, alpha::Float64=1.)
	T = copy(T_ref)
	n = length(T)
	# get random vector of norm-inf 1 (i.e. R belongs to the hypercube of side 1)
	R = get_rand_unit_vector_Ninf(n)
	for i in 1:n
		# the perturbative vector belongs to the hypercube of side 2*pi*alpha centered at the origin
		T[i] += (R[i]*2pi-pi)*alpha
	end
	return T
end

# initialize P vector with entry values in [-1, 1] and absolute value of greatest value equal to 1
#
# NB: U is assumed to be a distribution
function init_P1(U::Array{Float64,1})
	P = U*2-1
	# make sure that sum(P)=0
	P -= mean(P)
	# |P|inf = 1
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

# generate a complete graph of size n

function generate_complete_graph(n::Int)
	vs = Bus[]
	es = Line[]
	push!(vs,Bus(1,0.,0.))
	line_idx = 0
	for i in 2:n
		push!(vs,Bus(i,0.,0.))
		for j in 1:i-1
			line_idx += 1
			push!(es,Line(line_idx,vs[j],vs[i],-1.0im))
		end
	end
	
	return graph(vs,es,is_directed=false)
end

# generate a ring 
#
## INPUT
# n: length of the cycle
function generate_ring(n::Int)
	vs = Bus[]
	es = Line[]
	push!(vs,Bus(1,0.,0.))
	for i in 2:n
		push!(vs,Bus(i,0.,0.))
		push!(es,Line(i-1,vs[i-1],vs[i],-1.0im))
	end
	
	push!(es,Line(n,vs[n],vs[1],-1.0im))
	
	vs[1].bus_type = 3
	
	return graph(vs, es, is_directed=false)
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
#	Bus(1, "1", 3, 1., 1., 1., 0., load::Complex{Float64}, generation::Complex{Float64}, Q_min::Float64, Q_max::Float64, P_min::Float64, P_max::Float64, sh_conductance::Float64, sh_susceptance::Float64, lng::Float64, lat::Float64)
	for i in 2:n
		if i==ic
			push!(vs,Bus(i,0.,-p))
		else
			push!(vs,Bus(i,0.,0.))
		end
		push!(es,Line(i-1,vs[i-1],vs[i],-1.0im))
	end
	push!(es,Line(n,vs[n],vs[1],-1.0im))
	
	vs[1].bus_type = 3
		
	return graph(vs, es, is_directed=false)
end

# genereate a graph from its adjacency matrix
#
## INPUT
# adj: adjacency matrix
# p: power injections
function generate_graph_from_adj(adj::Array{Int,2},p::Array{Float64,1}=zeros(size(adj)[1]))
	vs = Bus[]
	es = Line[]
	n = size(adj)[1]
	nl = Int(0)
	
	for i in 1:n
		push!(vs,Bus(i,0.,p[i]))
	end
	
	for i in 1:n-1
		for j in i+1:n
			if adj[i,j] == 1
				nl += 1
				push!(es,Line(nl,vs[i],vs[j],-1.0im))
			end
		end
	end
	
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
			push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.0im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(m-2)
			push!(es, Line(ecounter,vs[j*n+i],vs[(j+1)*n+i],-1.0im))
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
#
# NB: returned angles belong to [-pi,pi]
function create_vortex_on_sq_lattice(n::Int,m::Int,i::Int,j::Int)
	T = zeros(Float64,n*m)
	for p in 1:(n*m)
		# get coordinates of the current point 
		# NB: here (x,y): (column,row)
		x = mod(p,m)
		if x == 0
			x = m
		end
		y = ceil(Int,p/m)
		# center of the vortex
		# NB: each square has a 1 unit length
		cx = j + 1/2
		cy = i + 1/2
		ba = atan(abs((y-cy)/(x-cx)))
		# Q2 (NE)
		if x < cx && y <= cy
			T[p] = -ba + pi
		# Q1 (NO)
		elseif x >= cx && y <= cy
			T[p] = ba 
		# Q4 (SO)
		elseif x >= cx && y > cy
			T[p] = -ba
		# Q3 (SE) (x < cx && y > cy)
		else
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
#
# NB: returned angles belong to [-pi,pi]
function create_antivortex_on_sq_lattice(n::Int,m::Int,i::Int,j::Int)
	T = zeros(Float64,n*m)
	for p in 1:(n*m)
		# get coordinates of the current point 
		# NB: here (x,y): (column,row)
		x = mod(p,m)
		if x == 0
			x = m
		end
		y = ceil(Int,p/m)
		# center of the vortex
		# NB: each square has a 1 unit length
		cx = j + 1/2
		cy = i + 1/2
		ba = atan(abs((y-cy)/(x-cx)))
		# Q2 (NE)
		if x < cx && y <= cy
			T[p] = ba - pi
		# Q1 (NO)
		elseif x >= cx && y <= cy
			T[p] = -ba
		# Q4 (SO)
		elseif x >= cx && y > cy
			T[p] = ba
		# Q3 (SE) (x < cx && y > cy)
		else
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

# compute the angle (mode A, i.e., with respect to y=0,x>0 axis) for the edge parent node - child node 
#
# NB: returned angle belongs to [-pi,pi[
function compute_edge_angle(px::Float64, py::Float64, cx::Float64, cy::Float64)
	#@debug("computing angle ($px,$py)-($cx,$cy)")
	dx = cx-px
	dy = cy-py
	t = atan(abs(dy/dx))
	# Q1
	if dx >= 0 && dy >= 0
		return t
	# Q2
	elseif dx < 0 && dy >= 0
		return -t+pi
	# Q3
	elseif dy < 0 && dx <= 0
		return t-pi
	# Q4 (dy < 0 && dx > 0)
	else
		return -t
	end
end

# compute the anti angle (mode A, i.e., with respect to y=0,x>0 axis) for the edge parent node - child node 
#
# NB: returned angle belongs to [-pi,pi[
function compute_anti_edge_angle(px::Float64, py::Float64, cx::Float64, cy::Float64)
	#@debug("computing anti angle ($px,$py)-($cx,$cy)")
	dx = cx-px
	dy = cy-py
	t = atan(abs(dy/dx))
	# Q1
	if dx >= 0 && dy >= 0
		return -t
	# Q2
	elseif dx < 0 && dy >= 0
		return t-pi
	# Q3
	elseif dy < 0 && dx <= 0
		return -t+pi
	# Q4 (dy < 0 && dx > 0)
	else
		return t
	end
end

# generate the contour cycle of a graph from (lat,lng) pairs
#
# NB: sink nodes need to be treated separately
#
# mode A: start from lowest node and visit nodes clockwise
# -> angles are measured anticlockwise with respect to y=0,x>0 axis
# -> for choosing next nodes, choose next child so that:
# 	-> angle is maximum if first crossed edge has been visited
#	-> angle is maximum before an already visited edge is crossed
#
# mode B: start from highest node and visit nodes anticlockwise
# -> angles are measured clockwise with respect to y=0,x<0 axis
# -> for choosing next nodes, apply same rule as in mode A
function get_geolocalized_graph_contour_cycle(g::Graphs.AbstractGraph{Bus,Line},ignore_vids::Set{Int}=Set{Int}())
	bcycle = Array{Int64,1}()

	X = Float64[]
	Y = Float64[]
	X2 = Float64[]
	Y2 = Float64[]
	nvids = Set{Int}()
	on_vid = Dict{Int,Int}()
	oid_v = Dict{Int,Bus}()

	# collect geolocalization data
	counter = 1
	for v in vertices(g)
		push!(X,v.lng)
		push!(Y,v.lat)
		# ignore sinks
		if out_degree(v,g) > 1
			push!(nvids,v.id)
			on_vid[v.id] = counter
			oid_v[v.id] = v
			push!(X2,v.lng)
			push!(Y2,v.lat)
			counter += 1
		end
	end

	# get lowest vertex (old index)
	lvi = indmin(Y)
	# get vertex from old index
	cv = oid_v[lvi]
	push!(bcycle,cv.id)
	while true
		#@debug("exploring vertex ", cv.id)
		# get children of current node
		nei = out_neighbors(cv,g)
		#@debug("nei: ", [v.id for v in nei])
		# no-sink children ids (old indices)
		cids = collect(setdiff(intersect(Set([v.id for v in nei]),nvids),ignore_vids))
		#@debug("cids: ", cids)
		# get parent node coordinates
		px = cv.lng
		py = cv.lat
		# compute angles for all no-sink children
		cangles = Float64[compute_edge_angle(px,py,X2[on_vid[o]],Y2[on_vid[o]]) for o in cids]
		# get index of the min angle
		mii = indmin(cangles)
		# first node was already explored
		if length(bcycle) > 1
			# get angle corresponding to current cv-pv edge
			pangle = cangles[find(x->x==pv.id,cids)[1]]
			#@debug("cv-pv edge angle: $pangle")
			#@debug("child edges angle: ", cangles)
			#@debug("min child edges angle: ", cangles[mii])
			# OR if the parent-child edge is the first one clockwise, select the maximum angle
			if cangles[mii] == pangle
				# get index of the max angle
				mai = indmax(cangles)
			# otherwise select the maximum angle lower than parent-child edge angle
			else
				# select maximum of angles < parent-child edge angle
				mfca = maximum(filter(x-> x < pangle,cangles))
				# get index of the max possible angle
				mai = find(x->x==mfca,cangles)[1]
			end
		else
			# get index of the max angle
			mai = indmax(cangles)
		end
		# set previous parent id for the next loop
		pv = cv
		# get new vertex
		cv = oid_v[cids[mai]]
		#@debug("parent  vertex ", pv.id ," (", pv.lng, " ", pv.lat, ")")
		#@debug("current vertex ", cv.id ," (", cv.lng, " ", cv.lat, ")")
		# if the cycle is closed
		if cv.id == bcycle[1]
			break
		else
			# add new node to the cycle
			push!(bcycle,cv.id)
		end
		#@debug("----------")
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
			push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.0im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(m-2)
			push!(es, Line(ecounter,vs[j*n+i],vs[(j+1)*n+i],-1.0im))
			ecounter += 1
		end
	end
	
	# stitch left and right sides
	for i in 0:(m-1)
		push!(es, Line(ecounter,vs[i*n+1],vs[n*(i+1)],-1.0im))
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
			push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.0im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(m-2)
			push!(es, Line(ecounter,vs[j*n+i],vs[(j+1)*n+i],-1.0im))
			ecounter += 1
		end
	end
	
	# stitch left and right sides
	for i in 0:(m-1)
		push!(es, Line(ecounter,vs[i*n+1],vs[n*(i+1)],-1.0im))
		ecounter += 1
	end
	
	# stitch up and bottom sides
	for i in 1:n
		push!(es, Line(ecounter,vs[i],vs[m*(n-1)+i],-1.0im))
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
			push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.0im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(4n-2)
			push!(es, Line(ecounter,vs[j*n+i],vs[(j+1)*n+i],-1.0im))
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
			push!(es, Line(ecounter,vs[4n^2+i*n+j],vs[4n^2+i*n+j+1],-1.0im))
			ecounter += 1
			push!(es, Line(ecounter,vs[5n^2+i*n+j],vs[5n^2+i*n+j+1],-1.0im))
			ecounter += 1
		end
	end

	# columns
	for i in 1:n
		for j in 0:(n-2)
			push!(es, Line(ecounter,vs[4n^2+j*n+i],vs[4n^2+(j+1)*n+i],-1.0im))
			ecounter += 1
			push!(es, Line(ecounter,vs[5n^2+j*n+i],vs[5n^2+(j+1)*n+i],-1.0im))
			ecounter += 1
		end
	end

	# stitch face 5 and 6 with the 4 first faces
	for i in 1:n
		# up
		push!(es, Line(ecounter,vs[4n^2+i],vs[2n^2+(i-1)*n+1],-1.0im))
		ecounter += 1
		push!(es, Line(ecounter,vs[5n^2+i],vs[3n^2-(i-1)*n],-1.0im))
		ecounter += 1
		# left
		push!(es, Line(ecounter,vs[4n^2+(i-1)*n+1],vs[2n^2-i*n+1],-1.0im))
		ecounter += 1
		push!(es, Line(ecounter,vs[5n^2+(i-1)*n+1],vs[3n^2+i*n],-1.0im))
		ecounter += 1
		# right
		push!(es, Line(ecounter,vs[4n^2+i*n],vs[3n^2+(i-1)*n+1],-1.0im))
		ecounter += 1
		push!(es, Line(ecounter,vs[5n^2+i*n],vs[2n^2-(i-1)*n],-1.0im))
		ecounter += 1
		# bottom
		push!(es, Line(ecounter,vs[5n^2-n+i],vs[n^2-i*n+1],-1.0im))
		ecounter += 1
		push!(es, Line(ecounter,vs[6n^2-n+i],vs[n^2-(i-1)*n],-1.0im))
		ecounter += 1
	end

	# stitch face 1 (up) with face 4 (bottom)
	for i in 1:n
		push!(es, Line(ecounter,vs[i],vs[4n^2-n+i],-1.0im))
		ecounter += 1
	end
	
	return graph(vs, es, is_directed=false)
end

# generate a flat triangular lattice of equilateral shape 
# with corners at positions (west, north-east, south-east)
# numbering of the sites is from bottom to top and from left to right
#
## INPUT
# n: number of nodes on the side of the lattice
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
			push!(es, Line(ecounter,vs[i],vs[i+1],-1.0im))
			#push!(es, Line(ecounter,vs[i*n+j],vs[i*n+j+1],-1.0im))
		end
		if (sites[i,2]+sites[i,3]) < n-1
			ecounter += 1
			push!(es, Line(ecounter,vs[i],vs[Int(sum(sites[i,1:3]))+1],-1.0im))
			ecounter += 1
			push!(es, Line(ecounter,vs[i],vs[Int(sum(sites[i,1:3]))+2],-1.0im))
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
		xvortex = xcoords[idx] - .5*tan(pi/6)
		yvortex = ycoords[idx] + .5
	else
		idx += Int((j)/2)
		xvortex = xcoords[idx] - .5*cos(pi/6)
		yvortex = ycoords[idx] 
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
# l,c,r: # of vertices on each branch of the double cycle
# p: produced/consumed power
function generate_double_ring(l::Int,c::Int,r::Int,p::Float64=0)
	vs = Bus[]
	es = Line[]	
	ecounter = 1

	# left branch 1:l
	for i in 1:l
		# t=0., p=0.
		push!(vs,Bus(i,0.,0.))
	end
	for i in 2:l
		push!(es, Line(ecounter,vs[i-1],vs[i],-1.0im))
		ecounter += 1
	end
	# central branch (l+1):(c+l)
	for i in (l+1):(c+l)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+2):(c+l)
		push!(es, Line(ecounter,vs[i-1],vs[i],-1.0im))
		ecounter += 1
	end
	# right branch (l+c+1):(c+l+r)
	for i in (l+c+1):(c+l+r)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+c+2):(c+l+r)
		push!(es, Line(ecounter,vs[i-1],vs[i],-1.0im))
		ecounter += 1
	end

	# 2 last remaining vertices
	push!(vs,Bus(l+c+r+1,0.,p))
	push!(vs, Bus(l+c+r+2,0.,-p))

	push!(es, Line(ecounter,vs[l+c+r+1],vs[1],-1.0im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+1],vs[l+c+1],-1.0im))
	ecounter += 1

	push!(es, Line(ecounter,vs[l+c+r+2],vs[l],-1.0im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+2],vs[l+c+r],-1.0im))
	ecounter += 1

	if c > 0
		push!(es, Line(ecounter,vs[l+c+r+1],vs[l+1],-1.0im))
		ecounter += 1
		push!(es, Line(ecounter,vs[l+c+r+2],vs[l+c],-1.0im))
	else
		push!(es, Line(ecounter,vs[l+c+r+1],vs[l+c+r+2],-1.0im))
	end
	return graph(vs, es, is_directed=false)
end
