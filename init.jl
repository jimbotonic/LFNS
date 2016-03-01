using Distributions

include("graphs.jl")

# initialize a uniform random distribution
function init_unif_dist(n::Int64)
	U = rand(Uniform(),n)
	return U/sum(U)
end

# initialize P vector with entry values taken uniformally at random in [-1, 1]
#
# U is assumed to be a distribution
function init_P1(U::Array{Float64,1})
	return (U*2-1)
end

# initialize P vector with entry values taken uniformally at random in [-1, 1]
#
# switch entries sign at random
# U is assumed to be a distribution
# absolute value of P's greatest value is 1
function init_P2(U::Array{Float64,1})
	n = length(U)
	# switch entry signs randomly
	P = rand([-1,1],n).*U
	P -= mean(P)
	P *= 1/maximum(abs(P))
	return P
end


# initialize P vector with entry values taken uniformally at random in [-1, 1]
#
# switch entries sign such that the highest value is positive, the 2nd highest value is negative, the 3rd is positive, ...
# U is assumed to be a distribution
# absolute value of P's greatest value is 1
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
# N: length of the cycle
# ic: index of the vertex of the consumer
# p: produced/consumed power
function generate_cycle(N::Int,ic::Int,p::Float64)
	vs = Bus[]
	es = Line[]
	
	push!(vs,Bus(1,0.,p))
	for i in 2:N
		if i==ic
			push!(vs,Bus(i,0.,-p))
		else
			push!(vs,Bus(i,0.,0.))
		end
		push!(es,Line(i-1,vs[i-1],vs[i],1.im))
	end
	push!(es,Line(N,vs[N],vs[1],1.im))
	
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
		push!(es, Line(ecounter,vs[i-1],vs[i],1.im))
		ecounter += 1
	end
	# central branch (l+1):(c+l)
	for i in (l+1):(c+l)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+2):(c+l)
		push!(es, Line(ecounter,vs[i-1],vs[i],1.im))
		ecounter += 1
	end
	# right branch (l+c+1):(c+l+r)
	for i in (l+c+1):(c+l+r)
		push!(vs,Bus(i,0.,0.))
	end
	for i in (l+c+2):(c+l+r)
		push!(es, Line(ecounter,vs[i-1],vs[i],1.im))
		ecounter += 1
	end

	# 2 last remaining vertices
	push!(vs,Bus(l+c+r+1,0.,p))
	push!(vs, Bus(l+c+r+2,0.,-p))

	push!(es, Line(ecounter,vs[l+c+r+1],vs[1],1.im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+1],vs[l+1],1.im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+1],vs[l+c+1],1.im))
	ecounter += 1

	push!(es, Line(ecounter,vs[l+c+r+2],vs[l],1.im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+2],vs[l+c],1.im))
	ecounter += 1
	push!(es, Line(ecounter,vs[l+c+r+2],vs[l+c+r],1.im))

	return graph(vs, es, is_directed=false)
end
