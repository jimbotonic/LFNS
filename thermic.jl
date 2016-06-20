include("graphs.jl")

# boiler settings
type BoilerSettings
	# discrete time at which the boiler is switched on
	switch_on_time::Int
end

# Boiler 
type Boiler
	is_on::Bool
	settings::BoilerSettings
	temp::Float64
	ref_temp::Float64
	comfort_temp_delta::Float64
	# temperature of the cold water
	cold_water_temp::Float64
	volume::Float64
	thermal_capacity::Float64
	thermal_conductivity::Float64
	power::Float64

	# constructor
	# volume, power, is_on, ref_temp, confort_temp_delta, cold_water_temp, thermal_conductivity
	# thermal capacity = 4185 J/(Kg x K) = 1.1625 Wh / (L x K)
	# temp -> random temp in the comfort zone
end

type HeatPump
	is_on::Bool
	temp::Float64
	ref_temp::Float64
	comfort_temp_delta::Float64
	volume::Float64
	thermal_capacity::Float64
	thermal_conductivity::Float64
	# "coefficient of performance" (thermal power = cop * power)
	cop::Float64
	power::Float64
	
	# constructor
	# thermal_capacity, thermal_conductivity, power, is_on, ref_temp, comfort_temp_delta, 
	# temp -> random temp in the comfort zone
end

# House profile time series
type HouseProfile
	# domestic appliances load
	domestic_app_load::Float64[]
	# photovoltaic production
	pv_production::Float64[]
	# hot water consumption
	hot_water_consumption::Float64[]
end

type House <: Bus
	profile::HouseProfile
	temp::Float64
	# windows surface
	win_surface::Float64
	pv_surface::Float64
	# house boiler
	boiler::Boiler
	# house heat pump
	heat_pump::HeatPump
	# are the blinds down?
	are_blinds_down::Bool
	
	# constructor
	# temp, win_surface,pv_surface, boiler, heat_pump
	# are_blinds_down = True
end

# weather time series
type Weather
	external_temp::Float64[]
	radiation_90::Float64[]
	radiation_40::Float64[]
end
