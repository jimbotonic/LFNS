using Distances

include("thermic.jl")
include("simulator.jl")

# right-hand side function of boiler solver
function rhs_boiler(house::House,consumption::Float64)
	boiler = house.boiler
	# debit of hot water (replaced by cold water)
	a= (consumption/boiler.volume)*(boiler.cold_temp-boiler.temp)
	# contribution of heating
	b=(boiler.is_on?boiler.power:0.)/(boiler.volume*boiler.thermal_capacity)
	# thermal losses to the house
	c=(boiler.thermal_conductivity/(boiler.thermal_capacity*boiler.volume))*(house.temp-boiler.temp)
	return a+b+c
end

## INPUT
# callback_func: callback function to be called at each iteration of the solver
#
## OUTPUT
function Euler_boiler_solver(sp::SParams,callback_func::Function)
	for i in 1:sp.iter_max
		for house in sp.neighborhood
			house.boiler.temp += sp.delta*rhs_boiler(house,house.profile.hot_water_consumption[i])
		end
		# call back (control and monitoring)
		callback_func(sp.controller,sp,i)
	end
end

## INPUT
# callback_func: callback function to be called at each iteration of the solver
#
## OUTPUT
function Euler_heatpump_solver(sp::SParams,callback_func::Function)
	# return State(Float64[],sp.T,Tdot,n_iter)
end
