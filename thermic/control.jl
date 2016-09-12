include("thermic.jl")

# Boiler callback function 
function th_device_controller(sp::SParams,i::Int64) 
	# thermostat control
	for house in sp.neighborhood
		if house.thermic_device.temp < (house.thermic_device.ref_temp - house.thermic_device.comfort_temp_delta)
			house.thermic_device.is_on = true
		end
		if house.thermic_device.temp > (house.thermic_device.ref_temp + house.thermic_device.comfort_temp_delta)
			house.thermic_device.is_on = false
		end
	end

	# save state
	P = Float64[house.thermic_device.is_on?house.thermic_device.power:0. for house in sp.neighborhood]
	K = Float64[house.thermic_device.temp for house in sp.neighborhood]
	push!(sp.states,State(P,K))	
end
