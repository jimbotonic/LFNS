using TimeSeries

# NB: by default, powers are expressed in Watt, times in hours, time steps in minute
# NB: arrays of temperatures are denoted by K
# NB: a house is assumed to have only one thermic device

# thermal capacity = 4185 J/(Kg x K) = 1.1625 Wh / (L x K)
WATER_THERMAL_CAPACITY = 1.1625

abstract ThermicDevice

# Boiler 
type Boiler <: ThermicDevice
	is_on::Bool
	temp::Float64
	ref_temp::Float64
	comfort_temp_delta::Float64
	# temperature of the cold water
	cold_temp::Float64
	volume::Float64
	power::Float64
	# thermal capacity = 4185 J/(Kg x K) = 1.1625 Wh / (L x K)
	thermal_conductivity::Float64
	thermal_capacity::Float64

	# constructor
	# volume, power, is_on, ref_temp, confort_temp_delta, cold_water_temp, thermal_conductivity
	# temp -> random temp in the comfort zone
	function Boiler(is_on::Bool,temp::Float64,ref_temp::Float64,comfort_temp_delta::Float64,cold_temp::Float64,volume::Float64,power::Float64,thermal_conductivity::Float64)
		new(is_on,temp,ref_temp,comfort_temp_delta,cold_temp,volume,power,thermal_conductivity,WATER_THERMAL_CAPACITY)
	end
end

type HeatPump <: ThermicDevice
	is_on::Bool
	temp::Float64
	ref_temp::Float64
	comfort_temp_delta::Float64
	# "coefficient of performance" (thermal power = cop * power)
	cop::Float64
	power::Float64
	
	# constructor
	# power, is_on, ref_temp, comfort_temp_delta, 
	# temp -> random temp in the comfort zone
	function HeatPump(is_on::Bool,temp::Float64,ref_temp::Float64,comfort_temp_delta::Float64,cop::Float64,power::Float64)
		new(is_on,temp,ref_temp,comfort_temp_delta,cop,power)
	end
end

# House profile time series
type HouseProfile
	# domestic appliances load
	domestic_app_load::TimeArray
	# hot water consumption
	hot_water_consumption::TimeArray

	# default constructor
	function HouseProfile(domestic_app_load::TimeArray,hot_water_consumption::TimeArray)
		new(domestic_app_load,hot_water_consumption)
	end
end

type House 
	profile::HouseProfile
	temp::Float64
	# windows surface
	win_surface::Float64
	pv_surface::Float64
	# house thermal properties
	thermal_capacity::Float64
	thermal_conductivity::Float64
	# house boiler / heat pump
	thermic_device::ThermicDevice
	# are the blinds down?
	are_blinds_down::Bool
	
	# constructor
	# temp, win_surface,pv_surface, boiler, heat_pump
	# are_blinds_down = True
	function House(profile::HouseProfile,temp::Float64,win_surface::Float64,pv_surface::Float64,thermal_capacity::Float64,thermal_conductivity::Float64,thermic_device::ThermicDevice,are_blinds_down::Bool)
		new(profile,temp,win_surface,pv_surface,thermal_capacity,thermal_conductivity,thermic_device,are_blinds_down)
	end
end

# weather time series
type Weather
	external_temp::TimeArray
	radiation_90::TimeArray
	radiation_40::TimeArray
end


