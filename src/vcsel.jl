idx = parse(Int, ARGS[1])
include(joinpath(@__DIR__, "isoform-genetics.jl"))
using VCSEL, JLD

if !isfile(joinpath(@__DIR__, "../data/grms_by_chr.jld"))
    grms = []
    for chr in string.(collect(1:22))
        @info "Working on chromosome $(chr)"
        path = joinpath(@__DIR__, "../data/PE-chr")
        SnpArrays.filter(geno, trues(size(geno)[1]), geno.snp_info.chromosome .== chr; des = path)
        geno_subset = SnpData(path)
        grm_subset = grm(geno_subset.snparray, method = :GRM)
        grm_subset .*= 2
        push!(grms, grm_subset)
    end
    save(joinpath(@__DIR__, "../data/grms_by_chr.jld"), "grms", grms)
    for plink in ["bed", "bim", "fam"]
        rm(joinpath(@__DIR__, "../data/PE-chr.$(plink)"))
    end    
end

grms = load(joinpath(@__DIR__, "../data/grms_by_chr.jld"), "grms")
grms = convert(Vector{Matrix{Float64}}, grms)
push!(grms, Matrix{Float64}(I, size(grms[1])))

uni = CSV.read("results/univariate.tsv", DataFrame)
filter!(row -> row.p < 0.05, uni)

if uni.feature[idx] == "gene"
    ind = findfirst(isequal(uni.id[idx]), expr.pid)
    Yg = Array{Float64}(expr[ind, 7:end])
    vcm = VCModel(Yg, cov, grms)
else
    ind = findfirst(isequal(uni.id[idx]), expri.gid)
    Yi = Array{Float64}(expri[ind, 7:end])
    vcm = VCModel(Yi, cov, grms)
end
# @time vcselect!(vcm; penfun = L1Penalty(), λ = 10)
@time Σ̂path, β̂path, λpath, objpath, niterspath = 
    vcselectpath!(vcm; penfun = L1Penalty(), nλ = 50);

soln = rankvarcomps(Σ̂path, resvarcomp = true)
df = DataFrame(fill(100, 1, 22), "chr" .* string.(collect(1:22)))
for chr in 1:22
    if in(chr, soln[1])
        ind = findfirst(isequal(chr), soln[1])
        df[1, chr] = ind
    end
end
df.id .= uni.id[idx]

ispath(joinpath(@__DIR__, "../results/penalized")) || 
    mkpath(joinpath(@__DIR__, "../results/penalized"))

CSV.write("results/penalized/$(uni.id[idx]).tsv", df, delim = "\t")

# df = CSV.read("results/penalized/$(uni.id[idx]).tsv", DataFrame, delim = "\t")

# begin
#     colors = cgrad(ColorSchemes.tol_rainbow, 22, categorical = true)
#     r = rankvarcomps(Σ̂path, resvarcomp = true)[1]
#     f = Figure(resolution = (612, 306))
#     axs = [Axis(f[1, i]) for i in 1:2]
#     for  i in 1:(size(Σ̂path, 1) - 1)
#         lines!(axs[1], λpath,  Σ̂path[i, :], linewidth = 0.75, color = colors[i])
#     end
#     axs[1].ylabel = "σ̂²ᵢ"
#     axs[1].xlabel = "λ"
#     CairoMakie.xlims!(axs[1], minimum(λpath), maximum(λpath))
#     CairoMakie.ylims!(axs[1], 0, maximum(Σ̂path[1:(end - 1), :]) * 10 / 9)
#     for  i in 1:(size(Σ̂path, 1) - 1)
#         lines!(axs[2], λpath,  Σ̂path[i, :], linewidth = 0.75, color = colors[i])
#     end
#     lines!(axs[2], λpath,  Σ̂path[end, :], linewidth = 0.75, color = "black")
#     axs[2].ylabel = "σ̂²ᵢ"
#     axs[2].xlabel = "λ"
#     CairoMakie.xlims!(axs[2], minimum(λpath), maximum(λpath))
#     CairoMakie.ylims!(axs[2], 0, maximum(Σ̂path) * 10 / 9)
#     rowgap!(f.layout, 5)
#     colgap!(f.layout, 5)
#     # rowsize!(f.layout, 1, Aspect(1, 0.75))
#     hidedecorations!(axs[1], ticklabels = false, ticks = false, label = false)
#     hidedecorations!(axs[2], ticklabels = false, ticks = false, label = false)
#     Legend(f[1, 1], [LineElement(color = colors[i], linestyle = nothing) for i in r],
#         "chr " .* string.(r),
#         tellwidth = false, tellheight = false, rowgap = 0, 
#         framevisible = false, halign = :right, valign = :top,
#         patchsize = (5, 1.5), linewidth = 1, padding = (10, 3, 3, 3))
#     resize_to_layout!(f)
#     save("$(uni.id[idx])-vcsel.png", f, px_per_unit = 4)
#     f
# end