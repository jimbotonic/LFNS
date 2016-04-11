include("../init.jl")
include("../graphs.jl")
include("../simulator.jl")
include("../solvers.jl")
include("../data.jl")
include("../metrics.jl")
include("../plotly.jl")

using Logging, ConfParser

@Logging.configure(level=DEBUG)

# loading config file
conf = ConfParse("../config.ini")
parse_conf!(conf)

# simple square lattice
n = 49
m = 49

g = generate_sq_lattice(n,m)

# initialize simulation parameters
sb = parse(Float64,retrieve(conf,"solvers","base_voltage"))
max_iter = round(Int,parse(Float64,retrieve(conf,"solvers","max_iter")))
epsilon = parse(Float64,retrieve(conf,"solvers","epsilon"))

o_args = Dict{Symbol,Any}()
o_args[:h] = parse(Float64,retrieve(conf,"rk","h"))
s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)

###
# injections initialization
###

alpha = 1.1
random = true
N = 1

if !random
	# initialize the injections with a uniform distribution
	U = init_unif_dist(n*m)
	# rescale distribution and set injections
	P_ref = init_P4(U)
else
	U = init_rand_dist(n*m)
	# rescale distribution and set injections
	#P_ref = init_P4(U)
	P_ref = init_P3(U)
end
P = P_ref*alpha
change_P(s.g,P)

###
# create vortices
###

if N == 1
	# vortex or antivortex in the middle
	# choose middle square
	i = 25
	j = 25
	T = create_vortex_on_sq_lattice2(n,m,i,j)
	#T = create_antivortex_on_sq_lattice2(n,m,i,j)
	@debug("angle ($i,$j): ", T[(i-1)*n+j])
	@debug("angle ($i,$(j+1)): ", T[(i-1)*n+j+1])
	@debug("angle ($(i+1),$j): ", T[i*n+j])
	@debug("angle ($(i+1),$(j+1)): ", T[i*n+j+1])
	
	trace1 = scatter_trace(Float64[],Float64[],"v1")
	traces = scatter_trace[]
	push!(traces,trace1)
	
	LP = Array{Array{Int,1}}(1)
	LP[1] = [i,j]
	
	push!(traces[1].X,i)
	push!(traces[1].Y,j)
elseif N == 2
	# create double vortex
	#i1 = 21; j1 = 21
	#i2 = 29; j2 = 29
	i1 = 25; j1 = 21
	i2 = 25; j2 = 29

	T1 = create_vortex_on_sq_lattice2(n,m,i1,j1)
	T2 = create_vortex_on_sq_lattice2(n,m,i2,j2)
	#T2 = create_antivortex_on_sq_lattice2(n,m,i2,j2)
	T = T1+T2

	trace1 = scatter_trace(Float64[],Float64[],"v1")
	trace2 = scatter_trace(Float64[],Float64[],"v2")

	traces = scatter_trace[]
	push!(traces,trace1)
	push!(traces,trace2)

	# array of the last vortices positions of the vortices
	LP = Array{Array{Int,1}}(2)
	LP[1] = [i1,j1]
	LP[2] = [i2,j2]

	push!(traces[1].X,i1)
	push!(traces[1].Y,j1)
	push!(traces[2].X,i2)
	push!(traces[2].Y,j2)
elseif N == 3
	# create double vortex
	i1 = 11; j1 = 11
	i2 = 11; j2 = 39
	i3 = 39; j3 = 25 

	T1 = create_vortex_on_sq_lattice2(n,m,i1,j1)
	T2 = create_vortex_on_sq_lattice2(n,m,i2,j2)
	T3 = create_vortex_on_sq_lattice2(n,m,i3,j3)
	T = T1+T2+T3

	trace1 = scatter_trace(Float64[],Float64[],"v1")
	trace2 = scatter_trace(Float64[],Float64[],"v2")
	trace3 = scatter_trace(Float64[],Float64[],"v3")

	traces = scatter_trace[]
	push!(traces,trace1)
	push!(traces,trace2)
	push!(traces,trace3)

	# array of the last vortices positions of the vortices
	LP = Array{Array{Int,1}}(3)
	LP[1] = [i1,j1]
	LP[2] = [i2,j2]
	LP[3] = [i3,j3]

	push!(traces[1].X,i1)
	push!(traces[1].Y,j1)
	push!(traces[2].X,i2)
	push!(traces[2].Y,j2)
	push!(traces[3].X,i3)
	push!(traces[3].X,j3)
elseif N == 4
	# create double vortex
	i1 = 11; j1 = 11
	i2 = 39; j2 = 39
	i3 = 11; j3 = 39
	i4 = 39; j4 = 11

	T1 = create_vortex_on_sq_lattice2(n,m,i1,j1)
	T2 = create_vortex_on_sq_lattice2(n,m,i2,j2)
	T3 = create_vortex_on_sq_lattice2(n,m,i3,j3)
	T4 = create_vortex_on_sq_lattice2(n,m,i4,j4)
	T = T1+T2+T3+T4

	trace1 = scatter_trace(Float64[],Float64[],"v1")
	trace2 = scatter_trace(Float64[],Float64[],"v2")
	trace3 = scatter_trace(Float64[],Float64[],"v3")
	trace4 = scatter_trace(Float64[],Float64[],"v4")

	traces = scatter_trace[]
	push!(traces,trace1)
	push!(traces,trace2)
	push!(traces,trace3)
	push!(traces,trace4)

	# array of the last vortices positions of the vortices
	LP = Array{Array{Int,1}}(4)
	LP[1] = [i1,j1]
	LP[2] = [i2,j2]
	LP[3] = [i3,j3]
	LP[4] = [i4,j4]

	push!(traces[1].X,i1)
	push!(traces[1].Y,j1)
	push!(traces[2].X,i2)
	push!(traces[2].Y,j2)
	push!(traces[3].X,i3)
	push!(traces[3].X,j3)
	push!(traces[4].Y,i4)
	push!(traces[4].Y,j4)
end

# get the contour cycle of the lattice
bcycle =  get_sq_lattice_contour_cycle(n,m)

# solver callback function
function callback_func(sp::SParams,n_iter::Int,error::Float64)
	@info("vorticity (iteration: $n_iter, error: $error): ", vorticity(sp.T,bcycle))
	A,B,V = find_vortices_in_sq_lattice(n,m,sp.T)
	for i in 1:length(LP)
		if LP[i] != [-1,-1]
			changed = false
			for j in 1:length(A)
				# we assume that vortices move by steps of distance 1 (infinite norm)
				d = maximum(abs(LP[i]-[A[j],B[j]]))
				if d <= 1
					# the vortex moved of 1 step
					if d == 1
						LP[i] = [A[j],B[j]]
						push!(traces[i].X,A[j])
						push!(traces[i].Y,B[j])
						@info("--- new pos of vortex $i: [", A[j], "," , B[j], "]")
					end
					changed = true
					break
				end
			end
			# the vortex has disappeared
			if !changed 
				@info("--- vortex $i has disappeared")
				LP[i] = [-1,-1]
			end
		end
	end
	return true
end

@info("start vorticity: ", vorticity(T,bcycle))
change_T(s.g,T)
state = simulation(s,callback_func)
@info("end vorticity: ", vorticity(state.T,bcycle))

# plot data
plot_scatter_data(traces,"scatter","markers+lines", "vortices_position", None)

