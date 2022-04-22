using CSV, DataFrames

@info "Parsing multivariate results"
files = readdir("results/multivariate")
mul = Vector{DataFrame}(undef, length(files))
for i in 1:length(files)
    rem(i, 100) == 0 ? println("parsed $(i) genes so far") : nothing
    mul[i] = CSV.read("results/multivariate/$(files[i])", delim = "\t", DataFrame)
end
mul = reduce(vcat, mul)
df = combine(groupby(mul, :gene_name), nrow => :count)
sort(df, order(:count, rev = true))
df[df.count .== 21, :]
CSV.write("results/multivariate.tsv", mul, delim = "\t")

@info "Parsing pairwise bivariate results"
files = readdir("results/bivariate")
bi = Vector{DataFrame}(undef, length(files))
for i in 1:length(files)
    rem(i, 100) == 0 ? println("parsed $(i) genes so far") : nothing
    bi[i] = CSV.read("results/bivariate/$(files[i])", delim = "\t", DataFrame)
end
bi = reduce(vcat, bi)
rename!(bi, :gene_id => :gene_name, :gene_name => :gene_id)
CSV.write("results/bivariate.tsv", bi, delim = "\t")

@info "Parsing univariate results"
files = readdir("results/univariate")
uni = Vector{DataFrame}(undef, length(files))
for i in 1:length(files)
    rem(i, 100) == 0 ? println("parsed $(i) genes so far") : nothing
    uni[i] = CSV.read("results/univariate/$(files[i])", delim = "\t", DataFrame)
end
uni = reduce(vcat, mul)
CSV.write("results/univariate.tsv", uni, delim = "\t")