include("../data.jl")
include("../plotly.jl")

pr = load_serialized("../data/eurogrid/dists/eurogrid_pc_pr.jld")
P = collect(load_csv_data("../data/eurogrid/eurogrid_pc_P.csv")[1])

println(pr[1:10])
println(P[1:10])

println("Pearson corr.: ", cor(pr,P))
println("Spearman corr.: ", corspearman(pr,P))

println("Pearson corr. (abs): ", cor(pr,abs(P)))
println("Spearman corr. (abs): ", corspearman(pr,abs(P)))

X = Float64[]
Y = Float64[]

for k in 1:length(pr)
	push!(X,pr[k])
	push!(Y,P[k])
end

# plot data
plot_scatter_data(X,Y,"scatter","markers", "PR_P-ref", None)
