include("graphs.jl")

function init_unif_dist(size::Int64)
	A = (rand(size)*2-1)
	return (A-sum(A)/size)
end

# initialize P vector with entry values taken uniformally at random in [-alpha*max_value, alpha*max_value]
function init_P1(A::Array{Float64,1}, alpha::Float64, max_value::Float64)
	return A*alpha*max_value
end

# generate a cycle with one producer at vertex 1 and one consumer at a chosen vertex
## INPUT
# N: length of the cycle
# ic: index of the vertex of the consumer
# p: produced/consumed power
#
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
## INPUT
# l,c,r: number of VERTICES on each branch of the double cycle
# p: produced/consumed power
#
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
