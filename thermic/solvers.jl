using Distances

include("thermic.jl")
include("simulator.jl")

# right-hand side function of boiler solver
function rhs_boiler(house::House,consumption::Float64)
	boiler = house.thermic_device
	# debit of hot water (replaced by cold water)
	a= (consumption/boiler.volume)*(boiler.cold_temp-boiler.temp)
	# contribution of heating
	b=(boiler.is_on?boiler.power:0.)/(boiler.volume*boiler.thermal_capacity)
	# thermal losses to the house
	c=(boiler.thermal_conductivity/(boiler.thermal_capacity*boiler.volume))*(house.temp-boiler.temp)
	return a+b+c
end

## INPUT
# sp: simulation parameters
#
## OUTPUT
function euler_boiler_solver(sp::SParams)
	for i in 1:sp.iter_max
		for house in sp.neighborhood
			house.thermic_device.temp += sp.delta * rhs_boiler(house, house.profile.hot_water_consumption.values[i,1])
		end
		# call back function (control and monitoring)
		sp.controller(sp,i)
	end
end

## INPUT
# sp: simulation parameters
#
## OUTPUT
function euler_heatpump_solver(sp::SParams)
end
