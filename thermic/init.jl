include("thermic.jl")
include("../data.jl")

###
# default init params (House)
###

# unit: Watt . HOUR / Kelvin
MIN_TH_CAPACITY = 13000.
MAX_TH_CAPACITY = 27000.

# unit: W / Kelvin
MIN_TH_CONDUCTIVITY = 200.
MAX_TH_CONDUCTIVITY = 400.

# unit: Celsius degree
MIN_TEMP_REF = 21.
MAX_TEMP_REF = 23.

# unit: Celsius degree
COMFORT_TEMP_DELTA = 1.5

# unit: square meters
MIN_WIN_SURFACE = 15.
MAX_WIN_SURFACE = 20.

###
# default init params (Boiler)
###

IS_ON_FRACTION = 0.3

# unit: Celsius degree
REF_TEMP = 55. 
COMFORT_TEMP_DELTA = 5.
COLD_TEMP = 10.
# unit: Liter
VOLUME_VALUES = Float64[200,250,300,350]
# unit: Watt / Kelvin (from TD J. Mayor)
TH_CONDUCTIVITY_VALUES = Float64[0.98,1.16,1.34,1.52]

# unit: Watt / Liter
MIN_POWER = 10.
MAX_POWER = 20.


REF_TIME = DateTime(2000,01,01,00,00,00)

function get_rand_value(min_value::Float64,max_value::Float64)
	return min_value+rand(Float64)*(max_value-min_value)
end

# NB: pv_surface in square meters is the same for all houses
function load_neighborhood(n::Int, pv_surface::Float64, house_profiles::Array{HouseProfile,1})
	houses = House[]
	for i in 1:n
		# thermic device (boilers)
		if rand() < IS_ON_FRACTION
			is_on = true
		else
			is_on = false
		end
		rint = rand(collect(1:length(VOLUME_VALUES)))
		volume = VOLUME_VALUES[rint]
		th_conductivity = TH_CONDUCTIVITY_VALUES[rint]
		power = volume*get_rand_value(MIN_POWER,MAX_POWER)
		temp = get_rand_value(REF_TEMP-COMFORT_TEMP_DELTA,REF_TEMP+COMFORT_TEMP_DELTA)
		boiler = Boiler(is_on,temp,REF_TEMP,COMFORT_TEMP_DELTA,COLD_TEMP,volume,power,th_conductivity)

		# house 
		temp = get_rand_value(MIN_TEMP_REF,MAX_TEMP_REF)
		win_surface = get_rand_value(MIN_WIN_SURFACE,MAX_WIN_SURFACE)
		th_capacity = get_rand_value(MIN_TH_CAPACITY,MAX_TH_CAPACITY)
		th_conductivity = get_rand_value(MIN_TH_CONDUCTIVITY,MAX_TH_CONDUCTIVITY)
		are_blinds_down = true
			
		house = House(house_profiles[i],temp,win_surface,pv_surface,th_capacity,th_conductivity,boiler,are_blinds_down)
		push!(houses,house)
	end
	
	return houses
end

function load_house_profiles(fn_domestic_app_load::AbstractString,fn_hot_water_consumption::AbstractString)
	house_profiles = HouseProfile[]

	domestic_app_loads = load_csv_data(fn_domestic_app_load)
	hot_water_consumptions = load_csv_data(fn_hot_water_consumption)
	s = size(domestic_app_loads)
	dates = DateTime[REF_TIME+Dates.Minute(i-1) for i in 1:s[1]]

	for j in 1:s[2]
		dal_ts = TimeArray(dates,domestic_app_loads[:,j])
		hwc_ts = TimeArray(dates,hot_water_consumptions[:,j])
		push!(house_profiles,HouseProfile(dal_ts,hwc_ts))
	end

	return house_profiles
end

function load_weather_ts(fn_temp::AbstractString,fn_radiation90::AbstractString,fn_radiation40::AbstractString)

end
