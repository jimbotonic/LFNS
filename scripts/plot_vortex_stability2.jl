using Logging

include("../data.jl")
include("../plotly.jl")

@Logging.configure(level=DEBUG)

# stats associated to a given random P
type AStats2
	# vortex stable & convergence
	Vst1::Dict{Float64,Float64}
	Vst2::Dict{Float64,Float64}
	# vortices inside & convergence
	Vin::Dict{Float64,Float64}
	# 1 and 2 vortices outside & convergence
	Vout1::Dict{Float64,Float64}
	Vout2::Dict{Float64,Float64}
	# no convergence
	Div::Dict{Float64,Float64}
end

# load graph
astats = load_serialized(ARGS[1])
pstats = load_serialized(ARGS[2])

plot_astats = false
plot_pstats = true

###
# Plotting AStats
###
if plot_astats
	trace1 = scatter_trace(Float64[], Float64[], "Vst1")
	trace2 = scatter_trace(Float64[], Float64[], "Vst2")
	trace3 = scatter_trace(Float64[], Float64[], "Vin")
	trace4 = scatter_trace(Float64[], Float64[], "Vout1")
	trace5 = scatter_trace(Float64[], Float64[], "Vout2")
	trace6 = scatter_trace(Float64[], Float64[], "Div")

	skeys = sort(collect(keys(astats.Vst1)))
	for k in skeys
		push!(trace1.X,k)
		push!(trace1.Y,astats.Vst1[k])
	end

	skeys = sort(collect(keys(astats.Vst2)))
	for k in skeys
		push!(trace2.X,k)
		push!(trace2.Y,astats.Vst2[k])
	end

	skeys = sort(collect(keys(astats.Vin)))
	for k in skeys
		push!(trace3.X,k)
		push!(trace3.Y,astats.Vin[k])
	end

	skeys = sort(collect(keys(astats.Vout1)))
	for k in skeys
		push!(trace4.X,k)
		push!(trace4.Y,astats.Vout1[k])
	end

	skeys = sort(collect(keys(astats.Vout2)))
	for k in skeys
		push!(trace5.X,k)
		push!(trace5.Y,astats.Vout2[k])
	end

	skeys = sort(collect(keys(astats.Div)))
	for k in skeys
		push!(trace6.X,k)
		push!(trace6.Y,astats.Div[k])
	end

	atraces = Array{scatter_trace,1}()
	push!(atraces,trace1)
	push!(atraces,trace2)
	push!(atraces,trace3)
	push!(atraces,trace4)
	push!(atraces,trace5)
	push!(atraces,trace6)

	# plot data
	plot_scatter_data(atraces,"scatter","markers", "A2", None)
end

function display_pstats(title::AbstractString,d::Dict{Float64,Int})
	sd = sort(collect(d))
	freqs = Int[sd[i][2] for i in 1:length(sd)]
	freqs /= sum(freqs)
	alphas = Float64[sd[i][1] for i in 1:length(sd)]
	esp = dot(freqs,alphas)
	variance = dot(freqs,alphas.^2)-esp^2
	println("### Stats for: ", title)
	println("-> mean: ", esp)
	println("-> std : ", sqrt(variance))
end

###
# PStats
###
if plot_pstats
	counter = 1
	
	avst1 = Dict{Float64,Int}()
	avst2 = Dict{Float64,Int}()
	avin = Dict{Float64,Int}()
	avout1 = Dict{Float64,Int}()
	avout2 = Dict{Float64,Int}()
	adiv = Dict{Float64,Int}()
	
	for P in keys(pstats)
		
		as = pstats[P]
		skas = sort(collect(keys(as)))
		
		for a in skas
			s = as[a]
			# checking level
			if s == :Vst1
				if haskey(avst1,a)
					avst1[a] += 1
				else
					avst1[a] = 1
				end
			elseif s == :Vst2
				if haskey(avst2,a)
					avst2[a] += 1
				else
					avst2[a] = 1
				end
			elseif s == :Vin
				if haskey(avin,a)
					avin[a] += 1
				else
					avin[a] = 1
				end
			elseif s == :Vout1
				if haskey(avout1,a)
					avout1[a] += 1
				else
					avout1[a] = 1
				end
			elseif s == :Vout2
				if haskey(avout2,a)
					avout2[a] += 1
				else
					avout2[a] = 1
				end
			elseif s == :Div
				if haskey(adiv,a)
					adiv[a] += 1
				else
					adiv[a] = 1
				end
			end
		end
	end
	display_pstats("Vst1", avst1)
	display_pstats("Vst2", avst2)
	display_pstats("Vin", avin)
	display_pstats("Vout1", avout1)
	display_pstats("Vout2", avout2)
	display_pstats("Div", adiv)
end
