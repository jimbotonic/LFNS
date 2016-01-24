
# initialize P vector with entry values taken uniformally at random in [-alpha*max_value, alpha*max_value]
function init_P1(size::Int64, alpha::Float64, max_value::Float64)
	P = (rand(size)*2-1)*alpha*max_value
	return (P-sum(P)/size)
end
