

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
