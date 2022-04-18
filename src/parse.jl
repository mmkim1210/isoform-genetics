using CSV, DataFrames

@info "Parsing multivariate results"
files = readdir(joinpath(@__DIR__, "../results/multivariate"))
results = Vector{DataFrame}(undef, length(files))
for i in 1:length(files)
    results[i] = CSV.read(joinpath(@__DIR__, "../results/multivariate/$(files[i])"), delim = "\t", DataFrame)
end
results = reduce(vcat, results)
df = combine(groupby(results, :gene_name), nrow => :count)
sort(df, order(:count, rev = true))
df[df.count .== 21, :]
CSV.write(joinpath(@__DIR__, "../results/multivariate.tsv"), results, delim = "\t")

@info "Parsing pariwise bivariate results"
files = readdir(joinpath(@__DIR__, "../results/bivariate"))
results = Vector{DataFrame}(undef, length(files))
for i in 1:length(files)
    results[i] = CSV.read(joinpath(@__DIR__, "../results/bivariate/$(files[i])"), delim = "\t", DataFrame)
end
results = reduce(vcat, results)
CSV.write(joinpath(@__DIR__, "../results/bivariate.tsv"), results, delim = "\t")

@info "Parsing univariate results"
files = readdir(joinpath(@__DIR__, "../results/univariate"))
results = Vector{DataFrame}(undef, length(files))
for i in 1:length(files)
    results[i] = CSV.read(joinpath(@__DIR__, "../results/univariate/$(files[i])"), delim = "\t", DataFrame)
end
results = reduce(vcat, results)
CSV.write(joinpath(@__DIR__, "../results/univariate.tsv"), results, delim = "\t")