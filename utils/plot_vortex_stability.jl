using Logging

include("../data.jl")
include("../plotly.jl")

@Logging.configure(level=DEBUG)

# stats associated to a given random P
type Stats
	# vortex stable & convergence
	Vst::Dict{Float64,Float64}
	# vortex inside & convergence
	Vin::Dict{Float64,Float64}
	# vortex outside & convergence
	Vout::Dict{Float64,Float64}
	# no convergence
	Div::Dict{Float64,Float64}
end

# load graph
astats = load_serialized(ARGS[1])
pstats = load_serialized(ARGS[2])

plot_astats = false
plot_pstats = true

###
# AStats
###
if plot_astats
	trace1 = scatter_trace(Float64[], Float64[], "Vst")
	trace2 = scatter_trace(Float64[], Float64[], "Vin")
	trace3 = scatter_trace(Float64[], Float64[], "Vout")
	trace4 = scatter_trace(Float64[], Float64[], "Div")

	# sort dict keys
	skeys = sort(collect(keys(astats.Vst)))

	# plot vorticity
	for k in skeys
		push!(trace1.X,k)
		push!(trace1.Y,astats.Vst[k])
	end

	# sort dict keys
	skeys = sort(collect(keys(astats.Vin)))

	# plot vorticity
	for k in skeys
		push!(trace2.X,k)
		push!(trace2.Y,astats.Vin[k])
	end

	# sort dict keys
	skeys = sort(collect(keys(astats.Vout)))

	# plot vorticity
	for k in skeys
		push!(trace3.X,k)
		push!(trace3.Y,astats.Vout[k])
	end

	# sort dict keys
	skeys = sort(collect(keys(astats.Div)))

	# plot vorticity
	for k in skeys
		push!(trace4.X,k)
		push!(trace4.Y,astats.Div[k])
	end

	atraces = Array{scatter_trace,1}()
	push!(atraces,trace1)
	push!(atraces,trace2)
	push!(atraces,trace3)
	push!(atraces,trace4)

	# plot data
	plot_scatter_data(atraces,"scatter","markers", "A", None)
end


function display_pstats(title::AbstractString,d::Dict{Float64,Int})
	sd = sort(collect(d))
	esd = 0.
	vsd = 0.
	ntot = sum(Int[sd[i][2] for i in 1:length(sd)])
	for p in sd
		esd += p[1]*(p[2]/ntot)
		vsd += p[1]^2*(p[2]/ntot)
	end
	dev = sqrt(vsd - esd^2)
	println("### Stats for: ", title)
	println("-> esperance: ", esd)
	println("-> deviation: ", dev)
end

###
# PStats
###
if plot_pstats
	counter = 1
	
	avst = Dict{Float64,Int}()
	avin = Dict{Float64,Int}()
	avout = Dict{Float64,Int}()
	adiv = Dict{Float64,Int}()
	cvst = Array{Float64,1}()
	cvin = Array{Float64,1}()
	cvout = Array{Float64,1}()
	backtracks = Array{Int,1}()
	
	for P in keys(pstats)
		level = 1
		backtracking = 0
		
		as = pstats[P]
		skas = sort(collect(keys(as)))
		
		v1_min = -1.
		v2_min = -1.
		v3_min = -1.
		v4_min = -1.
		v1_max = 0.
		v2_max = 0.
		v3_max = 0.
		v4_max = 0.
		for a in skas
			s = as[a]
			# checking level
			if s == :Vst
				if haskey(avst,a)
					avst[a] += 1
				else
					avst[a] = 1
				end
				clevel = 1
				if v1_min == -1
					v1_min = a
				end
				if a > v1_max
					v1_max = a
				end
			elseif s == :Vin
				if haskey(avin,a)
					avin[a] += 1
				else
					avin[a] = 1
				end
				clevel = 2
				if v2_min == -1
					v2_min = a
				end
				if a > v2_max
					v2_max = a
				end
			elseif s == :Vout
				if haskey(avout,a)
					avout[a] += 1
				else
					avout[a] = 1
				end
				clevel = 3
				if v3_min == -1
					v3_min = a
				end
				if a > v3_max
					v3_max = a
				end
			elseif s == :Div
				if haskey(adiv,a)
					adiv[a] += 1
				else
					adiv[a] = 1
				end
				clevel = 4
				if v4_min == -1
					v4_min = a
				end
				if a > v4_max
					v4_max = a
				end
			end
			# updating level if necessary
			if clevel < level 
				#@info("backtracking -> $s ($counter)")
				backtracking += 1
			elseif clevel > level
				level = clevel
			end
		end
		@info("# backtracking $backtracking")
		@info("Vst $v1_min / $v1_max")
		@info("Vin $v2_min / $v2_max")
		@info("Vout $v3_min / $v3_max")
		@info("Div $v4_min / $v4_max")
		push!(cvst, mean([v1_min,v1_max]))
		if v2_min != -1
			push!(cvin, mean([v2_min,v2_max]))
		end
		if v3_min != -1
			push!(cvout, mean([v3_min,v3_max]))
		end
		push!(backtracks, backtracking)
		counter += 1
	end
	@info("AVG # backtracking: ", mean(backtracks))
	display_pstats("Vst", avst)
	display_pstats("Vin", avin)
	display_pstats("Vout", avout)
	display_pstats("Div", adiv)
	@info("Vst center mean/deviation: ", mean(cvst), "/", std(cvst))
	@info("Vin center mean/deviation: ", mean(cvin), "/", std(cvin))
	@info("Vout center mean/deviation: ", mean(cvout), "/", std(cvout))
end
