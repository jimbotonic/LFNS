function init_unif_dist(size::Int64)
	A = (rand(size)*2-1)
	return (A-sum(A)/size)
end

# initialize P vector with entry values taken uniformally at random in [-alpha*max_value, alpha*max_value]
function init_P1(A::Array{Float64,1}, alpha::Float64, max_value::Float64)
	return A*alpha*max_value
end
