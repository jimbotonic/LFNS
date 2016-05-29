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

alpha = 0.6
random = true
N = 2

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

# dictionary of the last vortices positions of the vortices
global LP = Dict{Int,Array{Int,1}}()

if N == 1
	# vortex or antivortex in the middle
	# choose middle square
	i = 25
	j = 25
	T = create_vortex_on_sq_lattice(n,m,i,j)
	#T = create_antivortex_on_sq_lattice(n,m,i,j)
	#export_csv_data(T,"T_1-antivortex.csv")
	
	trace1 = scatter_trace(Float64[],Float64[],"v1")
	traces = scatter_trace[]
	push!(traces,trace1)
	
	LP[1] = [i,j]
	
	push!(traces[1].X,i)
	push!(traces[1].Y,j)
elseif N == 2
	# double vortex on the same row
	i1 = 25; j1 = 20
	i2 = 25; j2 = 31
	# double vortex on the same line
	#i1 = 20; j1 = 25
	#i2 = 31; j2 = 25
	# double vortex in the diagonal
	#i1 = 20; j1 = 20
	#i2 = 31; j2 = 31

	T1 = create_vortex_on_sq_lattice(n,m,i1,j1)
	T2 = create_vortex_on_sq_lattice(n,m,i2,j2)
	#T2 = create_antivortex_on_sq_lattice(n,m,i2,j2)
	T = T1+T2

	trace1 = scatter_trace(Float64[],Float64[],"v1")
	trace2 = scatter_trace(Float64[],Float64[],"v2")

	traces = scatter_trace[]
	push!(traces,trace1)
	push!(traces,trace2)

	LP[1] = [i1,j1]
	LP[2] = [i2,j2]

	push!(traces[1].X,i1)
	push!(traces[1].Y,j1)
	push!(traces[2].X,i2)
	push!(traces[2].Y,j2)
elseif N == 3
	# 3 vortices in diagonal
	#i1 = 23; j1 = 23
	#i2 = 25; j2 = 25
	#i3 = 27; j3 = 27 
	# 3 vortices on the line
	i1 = 25; j1 = 21
	i2 = 25; j2 = 25
	i3 = 25; j3 = 29 

	T1 = create_vortex_on_sq_lattice(n,m,i1,j1)
	#T2 = create_vortex_on_sq_lattice(n,m,i2,j2)
	T2 = create_antivortex_on_sq_lattice(n,m,i2,j2)
	T3 = create_vortex_on_sq_lattice(n,m,i3,j3)
	T = T1+T2+T3

	trace1 = scatter_trace(Float64[],Float64[],"v1")
	trace2 = scatter_trace(Float64[],Float64[],"v2")
	trace3 = scatter_trace(Float64[],Float64[],"v3")

	traces = scatter_trace[]
	push!(traces,trace1)
	push!(traces,trace2)
	push!(traces,trace3)

	LP[1] = [i1,j1]
	LP[2] = [i2,j2]
	LP[3] = [i3,j3]

	push!(traces[1].X,i1)
	push!(traces[1].Y,j1)
	push!(traces[2].X,i2)
	push!(traces[2].Y,j2)
	push!(traces[3].X,i3)
	push!(traces[3].Y,j3)
elseif N == 4
	# 4 vortices on a square -> gives vorticity=3
	#i1 = 21; j1 = 21
	#i2 = 29; j2 = 29
	#i3 = 21; j3 = 29
	#i4 = 29; j4 = 21
	# 4 vortices in the diagonal
	#i1 = 22; j1 = 22
	#i2 = 24; j2 = 24
	#i3 = 26; j3 = 26
	#i4 = 28; j4 = 28
	# 4 vortices center + triangle
	i1 = 25; j1 = 25
	i2 = 21; j2 = 25
	i3 = 29; j3 = 29
	i4 = 29; j4 = 21

	#T1 = create_vortex_on_sq_lattice(n,m,i1,j1)
	T1 = create_antivortex_on_sq_lattice(n,m,i1,j1)
	T2 = create_vortex_on_sq_lattice(n,m,i2,j2)
	T3 = create_vortex_on_sq_lattice(n,m,i3,j3)
	T4 = create_vortex_on_sq_lattice(n,m,i4,j4)
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

	LP[1] = [i1,j1]
	LP[2] = [i2,j2]
	LP[3] = [i3,j3]
	LP[4] = [i4,j4]

	push!(traces[1].X,i1)
	push!(traces[1].Y,j1)
	push!(traces[2].X,i2)
	push!(traces[2].Y,j2)
	push!(traces[3].X,i3)
	push!(traces[3].Y,j3)
	push!(traces[4].X,i4)
	push!(traces[4].Y,j4)
elseif N == 5
	# 4 vortices on a square -> gives vorticity=3
	i1 = 21; j1 = 21
	i2 = 29; j2 = 29
	i3 = 21; j3 = 29
	i4 = 29; j4 = 21
	i5 = 25; j5 = 25

	T1 = create_vortex_on_sq_lattice(n,m,i1,j1)
	T2 = create_vortex_on_sq_lattice(n,m,i2,j2)
	T3 = create_vortex_on_sq_lattice(n,m,i3,j3)
	T4 = create_vortex_on_sq_lattice(n,m,i4,j4)
	T5 = create_vortex_on_sq_lattice(n,m,i5,j5)
	T = T1+T2+T3+T4+T5

	trace1 = scatter_trace(Float64[],Float64[],"v1")
	trace2 = scatter_trace(Float64[],Float64[],"v2")
	trace3 = scatter_trace(Float64[],Float64[],"v3")
	trace4 = scatter_trace(Float64[],Float64[],"v4")
	trace5 = scatter_trace(Float64[],Float64[],"v5")

	traces = scatter_trace[]
	push!(traces,trace1)
	push!(traces,trace2)
	push!(traces,trace3)
	push!(traces,trace4)
	push!(traces,trace5)

	LP[1] = [i1,j1]
	LP[2] = [i2,j2]
	LP[3] = [i3,j3]
	LP[4] = [i4,j4]
	LP[5] = [i5,j5]

	push!(traces[1].X,i1)
	push!(traces[1].Y,j1)
	push!(traces[2].X,i2)
	push!(traces[2].Y,j2)
	push!(traces[3].X,i3)
	push!(traces[3].Y,j3)
	push!(traces[4].X,i4)
	push!(traces[4].Y,j4)
	push!(traces[5].X,i5)
	push!(traces[5].Y,j5)
end

# get the contour cycle of the lattice
bcycle =  get_sq_lattice_contour_cycle(n,m)

# record T vectors 
TS = Array{Array{Float64,1},1}()

# solver callback function
function callback_func(sp::SParams,n_iter::Int,error::Float64)
	@info("vorticity (iteration: $n_iter, error: $error): ", vorticity(sp.T,bcycle))

	global TS
	if n_iter < 10000 && n_iter % 100 == 0
		push!(TS,sp.T)
	end

	A,B,V = find_vortices_in_sq_lattice(n,m,sp.T)
	for k in keys(LP)
		changed = false
		for j in 1:length(A)
			# we assume that vortices move by steps of distance 1 (L1 norm)
			d = sum(abs(LP[k]-[A[j],B[j]]))
			if d == 1
				# the vortex moved of 1 step
				# update position
				LP[k] = [A[j],B[j]]
				push!(traces[k].X,A[j])
				push!(traces[k].Y,B[j])
				@info("--- new pos of vortex $k: [", A[j], "," , B[j], "]")
				changed = true
				break
			elseif d == 0
				# the vortex hasn't moved
				changed = true
				break
			end
		end
		# the vortex has disappeared
		if !changed 
			@info("--- vortex $k has disappeared")
			delete!(LP,k)
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

# export T vectors
#for i in 1:length(TS)
#	export_csv_data(TS[i],"./movie/$i.csv")	
#end
