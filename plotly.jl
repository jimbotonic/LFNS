using Plotly

### fitting

# return the linear fitting parameters (normal equation)
function get_linear_fit(aX,aY)
	# add a column of 1's to X matrix
	X = hcat(aX,ones(size(aX,1)))
	return pinv(X'*X)*X'*aY
end

# remove null entries before taking the log
function get_filtered_log(X,Y,minx=None,maxx=None)
	nX = Float64[]
	nY = Float64[]
	for i in 1:length(X)
		x = X[i]
		y = Y[i]

		nentry = false
		if x > 0 && y > 0
			if minx != None
				if maxx != None
					if x >= minx && x <= maxx
						nentry = true
					end
				else
					if x >= minx
						nentry = true
					end
				end
			else
				nentry = true
			end
		end
		if nentry
			push!(nX,log(x))
			push!(nY,log(y))
		end
	end
	return nX,nY
end

### layouts

function get_layout(title::String,xname::String,yname::String,xtype::String,ytype::String)
	return [
  "title" => title,
  "xaxis" => [
      "title" => xname,
      "type" => xtype, 
      "autorange" => true
	    ], 
  "yaxis" => [
      "title" => yname,
      # linear, log, date, category
      "type" => ytype, 
      "autorange" => true
		]
	]
end

### plotting functions

function plot_heatmap(bin_matrix, filename::String, log_scale=false, layout=None)
	if log_scale 
		bin_matrix = log(bin_matrix)
	end
	data = [
	  	[
	        "z" => bin_matrix,
		"type" => "heatmap"
		]
	]

	if layout != None
		response = Plotly.plot([data], ["layout" => layout, "filename" => filename, "fileopt" => "overwrite"])
	else
		response = Plotly.plot([data], ["filename" => filename, "fileopt" => "overwrite"])
	end
	plot_url = response["url"]
end

# send scatter plot data (one trace)
function plot_scatter_data(X::Array{Float64,1},Y::Array{Float64,1}, ptype::String, mode::String, filename::String, layout=None)
	data = [
	  [
	    "x" => X, 
	    "y" => Y, 
	    # scatter: y=f(x) points
	    # histogram2d: 2D histogram
	    "type" => ptype,
	    # markers, lines+markers
	    "mode" => mode
	  ]
	]

	if layout != None
		response = Plotly.plot([data], ["layout" => layout, "filename" => filename, "fileopt" => "overwrite"])
	else
		response = Plotly.plot([data], ["filename" => filename, "fileopt" => "overwrite"])
	end
	plot_url = response["url"]
end

# scatter trace type
type scatter_trace
	X::Array{Float64,1}
	Y::Array{Float64,1}
	name::String
end

# generate plotting data from an array of scatter traces
function get_scatter_data(traces::Array{scatter_trace,1},ptype::String,mode::String)
	data = Dict{String,Any}[]
	for t in traces
		ptrace = [
			"x" => t.X,
			"y" => t.Y,
			"type" => ptype,
			"mode" => mode,
			"name" => t.name
			]
		push!(data,ptrace)
	end
	return data
end

# send scatter plot data (multiple traces)
function plot_scatter_data(traces::Array{scatter_trace,1}, ptype::String, mode::String, filename::String, layout=None)
	data = get_scatter_data(traces,ptype,mode)
	if layout != None
		response = Plotly.plot([data], ["layout" => layout, "filename" => filename, "fileopt" => "overwrite"])
	else
		response = Plotly.plot([data], ["filename" => filename, "fileopt" => "overwrite"])
	end
	plot_url = response["url"]
end

# scatter trace type
type hist_trace
	X::Array{Float64,1}
	name::String
end

# generate plotting data from an array of histogram traces
function get_hist_data(traces::Array{hist_trace,1})
	data = Dict{String,Any}[]
	for t in traces
		ptrace = [
			"x" => t.X,
			"name" => t.name,
	  		"opacity" => 0.75,
	  		"histnorm" => "probability density",
	  		"type" => "histogram"
			]
		push!(data,ptrace)
	end
	return data
end

# plot 2 histograms
function plot_histograms(traces::Array{hist_trace,1}, filename::String, layout=None)
	data = get_hist_data(traces)

	#layout = ["title" => "Core Colinks Distribution",
	#  "xaxis" => ["title" => "# colinks"],
	#  "yaxis" => ["title" => "# vertices"],
	#  "barmode" => "overlay"]

	if layout != None
		response = Plotly.plot([data], ["layout" => layout, "filename" => filename, "fileopt" => "overwrite"])
	else
		response = Plotly.plot([data], ["filename" => filename, "fileopt" => "overwrite"])
	end
	plot_url = response["url"]
end
