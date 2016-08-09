include("thermic.jl")

# Boiler callback function 
function boiler_callback_func(control_func::Function,sp::SParams,i::Int64) 
	# control
	control_func(sp,i)

	# save state
	P = Float64[boiler.is_on?boiler.power:0. for boiler in sp.neighborhood]
	K = Float64[boiler.temp for boiler in sp.neighborhood]
	push!(sp.states,State(P,K)	
end
