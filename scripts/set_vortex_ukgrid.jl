using Logging, ConfParser

include("../init.jl")
include("../data.jl")
include("../graphs.jl")
include("../metrics.jl")
include("../simulator.jl")
include("../solvers.jl")

@Logging.configure(level=DEBUG)

# loading config file
conf = ConfParse("../config.ini")
parse_conf!(conf)

# load graph
g = load_serialized(ARGS[1])
n = length(vertices(g))

###
# computing PR
###

pr = ones(Float64,n)
pr = pr/n
pr = my_pagerank(g,pr,0.85,1e-11)

###
# getting contour cycle
###

# contour cycle starting from lowest node
bcycle = [1,5,4,6,7,38,39,40,41,37,35,117,49,50,51,52,54,68,69,77,78,80,81,89,95,97,96,100,101,102,104,106,108,109,107,94,92,91,90,87,86,83,82,75,74,60,59,31,29,28,22,20,16,15,14,13,12,10,119,11,3,2]

# correct lng of vertex 108 and 58
#vs = vertices(g)
#vs[108].lng = vs[107].lng
#vs[58].lng = vs[59].lng

# remove vertices (cut component 109 + branch 84-83) 
#ignore_vids = Set{Int}([84,110,111,112,113,114,115,116])

# get contour cycle
#bcycle = get_geolocalized_graph_contour_cycle(g,ignore_vids)

###
# injections initialization
###

alpha = 0.25
# :PR, :uniform, :random
init_type = :uniform

if init_type == :uniform
	# initialize the injections with a uniform distribution
	U = init_unif_dist(n)
	# rescale distribution and set injections
	P_ref = init_P4(U)
elseif init_type == :PR
	# rescale distribution and set injections
	P_ref = init_P4(pr)
else 
	U = init_rand_dist(n)
	# rescale distribution and set injections
	P_ref = init_P3(U)
end

P = P_ref*alpha

###
# initializing angles
###

X = Float64[]
Y = Float64[]
X2 = Float64[]
Y2 = Float64[]
nvids = Set{Int}()
on_vid = Dict{Int,Int}()
oid_v = Dict{Int,Bus}()
T = Float64[]

counter = 1
for v in vertices(g)
	push!(X,v.lat)
	push!(Y,v.lng)
	# ignore sinks
	if out_degree(v,g) > 1
		push!(nvids,v.id)
		on_vid[v.id] = counter
		oid_v[v.id] = v
		push!(X2,v.lat)
		push!(Y2,v.lng)
		counter += 1
	end
end

#xc,yc = mean(X2),mean(Y2)
xc,yc = mean(X),mean(Y)
@info("center coordinates ($xc,$yc)")

# compute angles 
for i in 1:length(X)
	push!(T,compute_edge_angle(xc, yc, X[i], Y[i]))
end

###
# launching simulation
###

# initialize simulation parameters
sb = parse(Float64,retrieve(conf,"solvers","base_voltage"))
max_iter = round(Int,parse(Float64,retrieve(conf,"solvers","max_iter")))
epsilon = parse(Float64,retrieve(conf,"solvers","epsilon"))

o_args = Dict{Symbol,Any}()
o_args[:h] = parse(Float64,retrieve(conf,"rk","h"))
s = Simulator(g,RK_solver1,o_args,sb,epsilon,max_iter)

@info("start vorticity: ", vorticity(T,bcycle))
change_P(s.g,P)
change_T(s.g,T)
state = simulation(s)
@info("end vorticity: ", vorticity(state.T,bcycle))
