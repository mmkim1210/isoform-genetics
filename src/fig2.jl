include(joinpath(@__DIR__, "isoform-genetics.jl"))

@info "Loading regressed isoform expression data"
@time regressed = Arrow.Table(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-seqPCs-isoform.BED.arrow")) |> DataFrame
@assert (93_293, 861) == size(regressed)

results = DataFrame[]
for model in ["univariate", "bivariate", "multivariate"]
    push!(results, CSV.read("results/$(model).tsv", DataFrame))
end

resultsg = filter(row -> row.feature == "gene", results[1])
@assert size(resultsg, 1) == 22965
resultsi = filter(row -> row.feature == "isoform", results[1])
@assert size(resultsi, 1) == 89925

df1 = results[3][Not(isnan.(results[3].p_res)), :]
select!(df1, [3; 5:16; 4])
df2 = DataFrame()
genes_bionly = setdiff(unique(results[2].id), results[3].id)
for gene in genes_bionly
    ind = findfirst(results[2].id .== gene)
    storage = results[2][ind, [3, 6, 8, 10, 12, 14, 16, 23, 24, 25, 26, 27, 28, 4]]
    df2 = vcat(df2, DataFrame(storage))
end
rename!(df2, names(df1))
df_rg = vcat(df1, df2)
df_rg.rp .= 0.0
for i in 1:nrow(df_rg)
    rem(i, 100) == 0 ? println(i) : nothing
    ind1 = findfirst(isequal(df_rg.id[i]), regressed.gid)
    ind2 = findfirst(isequal(df_rg.pair[i]), regressed.gid)
    expr1 = Array(regressed[ind1, 7:end])
    expr2 = Array(regressed[ind2, 7:end])
    expr1 = vec(expr1)
    expr2 = vec(expr2)
    df_rg.rp[i] = cor(expr1, expr2)
end

df_rg_mul = results[3][Not(isnan.(results[3].p_res)), :]
filter!(row -> row.id != "ENST00000002125", df_rg_mul)
filter!(row -> row.pair != "ENST00000002125", df_rg_mul)
select!(df_rg_mul, 3:16)
df_rg_bi = DataFrame()
for i in 1:size(df_rg_mul, 1)
    ind = findfirst(results[2].id .== df_rg_mul.id[i] .&& results[2].pair .== df_rg_mul.pair[i])
    storage = results[2][ind, [3, 4, 6, 8, 10, 12, 14, 16, 23, 24, 25, 26, 27, 28]]
    df_rg_bi = vcat(df_rg_bi, DataFrame(storage))
end

df1 = DataFrame()
for gene in unique(results[3].id)
    ind = findfirst(results[3].id .== gene .&& results[3].pair .== gene)
    storage = results[3][ind, [3; 5:16]]
    df1 = vcat(df1, DataFrame(storage))
end
df2 = DataFrame()
for gene in genes_bionly
    ind = findfirst(results[2].id .== gene)
    storage = results[2][ind, [3, 5, 7, 9, 11, 13, 15, 17, 18, 19, 20, 21, 22]]
    df2 = vcat(df2, DataFrame(storage))
end
rename!(df2, names(df1))
df_h2 = vcat(df1, df2)
df_h2 = df_h2[Not(18), :]
df_uni = DataFrame()
for iso in df_h2.id
    ind = findfirst(isequal(iso), resultsi.id)
    storage = resultsi[ind, [3; 5:16]]
    df_uni = vcat(df_uni, DataFrame(storage))
end
@assert size(df_h2, 1) == 11338

hist(results[3].p_res[Not(isnan.(results[3].p_res))], bins = 300)
hist(results[3][Not(isnan.(results[3].p_res)), "h2/re_reml"], bins = 300)

let
    h²sg = [resultsg[:, "h2_cis_reml"], resultsg[:, "h2_trans_reml"], resultsg[:, "h2_reml"]]
    h²si = [resultsi[:, "h2_cis_reml"], resultsi[:, "h2_trans_reml"], resultsi[:, "h2_reml"]]
    colors = ["#B8E3FF", "#0583D2", "#16558F", "#61B0B7"]
    # https://colors.dopely.top/inside-colors/blue-color-palette-inspirations-with-names-and-hex-codes/
    f = Figure(resolution = (530, 792))
    g1 = f[1, 1] = GridLayout()
    g2 = f[2, 1] = GridLayout()
    g3 = f[3, 1] = GridLayout()
    g4 = f[1, 2] = GridLayout()
    g5 = f[2, 2] = GridLayout()
    g6 = f[3, 2] = GridLayout()
    colgap!(f.layout, 5)
    rowgap!(f.layout, 2.5)
    @info "Distribution of h² estimates for univariate model"
    axs1 = [Axis(g1[1, j]; xticks = -1:0.5:1) for j in 1:3]
    for (i, h²) in enumerate(h²sg)
        hist!(axs1[1], h², offset = -2i, scale_to = 1.5, bins = 100, color = colors[i])
        boxplot!(axs1[1], fill(-2i - 0.3, length(h²)), h², orientation = :horizontal,
            show_outliers = false, width = 0.4, color = colors[i])
    end
    for (i, h²) in enumerate(h²si)
        hist!(axs1[2], h², offset = -2i, scale_to = 1.5, bins = 100, color = colors[i])
        boxplot!(axs1[2], fill(-2i - 0.3, length(h²)), h², orientation = :horizontal,
            show_outliers = false, width = 0.4, color = colors[i])
    end
    Label(g1[0, 1:2], text = "Distribution of h²SNP estimates", textsize = 8, tellwidth = false)
    Label(g1[1, 3], "Dist. of rg estimates", textsize = 8, tellwidth = false)
    Label(g1[2, 1, Top()], "Genes", textsize = 6)
    Label(g1[2, 2, Top()], "Isoforms", textsize = 6)
    Label(g1[2, 3, Top()], "Isoforms", textsize = 6)
    @info "Distribution of r₉ estimates for multivariate model"
    for i in 1:3
        hist!(axs1[3], df_rg[:, 7 + i], offset = -2i, scale_to = 1.5, bins = 100, color = colors[i])
        boxplot!(axs1[3], fill(-2i - 0.3, size(df_rg, 1)), df_rg[:, 7 + i], orientation = :horizontal,
            show_outliers = false, width = 0.4, color = colors[i])
    end
    [axs1[j].yticks = ([-2, -4, -6], fill("", 3)) for j in 1:3]
    [hidexdecorations!(axs1[j], ticks = false, ticklabels = false) for j in 1:3]
    [hideydecorations!(axs1[j]) for j in 1:3]
    [xlims!(axs1[j], -0.05, 1.05) for j in 1:2]
    xlims!(axs1[3], -1.05, 1.05)
    rowsize!(g1, 2, Aspect(3, 1))
    Label(g1[1, 1:2, Left()], "a", textsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    Label(g1[1, 3, Left()], "d", textsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    colgap!(g1, 2.5)
    rowgap!(g1, 0)
    @info "Comparison of variance components estimates between univariate and multivariate models"
    axs2 = [Axis(g2[1, j]) for j in 1:3]
    [scatter!(axs2[j], df_uni[:, 1 + j], df_h2[:, 1 + j], 
        color = (colors[j], 0.5), markersize = 2) for j in 1:3]
    [abline!(axs2[j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    for j in 1:3
        data = DataFrame(x = df_uni[:, 1 + j ], y = df_h2[:, 1 + j])
        model = lm(@formula(y ~ x), data)
        abline!(axs2[j], coef(model)[1], coef(model)[2], color = :gold, linewidth = 0.5)
    end
    [hidedecorations!(axs2[j], ticks = false, ticklabels = false) for j in 1:3]
    maxs = [maximum(vcat(df_h2[:, 1 + j], df_uni[:, 1 + j])) for j in 1:3]
    [ylims!(axs2[j], 0, maxs[j] * 6 / 5) for j in 1:3]
    [xlims!(axs2[j], 0, maxs[j] * 6 / 5) for j in 1:3]
    rowsize!(g2, 1, Aspect(3, 1))
    Label(g2[0, 1:3], text = "Comparison of variance components estimates", textsize = 8, tellwidth = false)
    Label(g2[2, 1, Top()], "Cis-genetic variance", textsize = 6)
    Label(g2[2, 2, Top()], "Trans-genetic variance", textsize = 6)
    Label(g2[2, 3, Top()], "Residual variance", textsize = 6)
    Label(g2[1, 1:3, Left()], "b", textsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    Label(g2[3, 1:3], text = "Univariate estimates", textsize = 6)
    Label(g2[2, 0], text = "Multivariate estimates", textsize = 6, rotation = π / 2, tellheight = false)
    colgap!(g2, 2.5)
    rowgap!(g2, 0)
    @info "Comparison of h² estimates between univariate and multivariate models"
    axs3 = [Axis(g3[1, j]) for j in 1:3]
    [scatter!(axs3[j], df_uni[:, 7 + j], df_h2[:, 7 + j], 
        color = (colors[j], 0.5), markersize = 2) for j in 1:3]
    [abline!(axs3[j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    for j in 1:3
        data = DataFrame(x = df_uni[:, 7 + j ], y = df_h2[:, 7 + j])
        model = lm(@formula(y ~ x), data)
        abline!(axs3[j], coef(model)[1], coef(model)[2], color = :gold, linewidth = 0.5)
    end
    [hidedecorations!(axs3[j], ticks = false, ticklabels = false) for j in 1:3]
    [ylims!(axs3[j], -0.05, 1.05) for j in 1:3]
    [xlims!(axs3[j], -0.05, 1.05) for j in 1:3]
    rowsize!(g3, 1, Aspect(3, 1))
    Label(g3[0, 1:3], text = "Comparison of h²SNP estimates", textsize = 8, tellwidth = false)
    Label(g3[2, 1, Top()], "h2cis", textsize = 6)
    Label(g3[2, 2, Top()], "h2trans", textsize = 6)
    Label(g3[2, 3, Top()], "h2SNP", textsize = 6)
    Label(g3[1, 1:3, Left()], "c", textsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    Label(g3[3, 1:3], text = "Univariate estimates", textsize = 6)
    Label(g3[2, 0], text = "Multivariate estimates", textsize = 6, rotation = π / 2, tellheight = false)
    colgap!(g3, 2.5)
    rowgap!(g3, 0)
    @info "Comparison of variance components estimates between bivariate and multivariate models"
    axs4 = [Axis(g4[1, j]) for j in 1:3]
    [scatter!(axs4[j], df_rg_bi[:, 2 + j], df_rg_mul[:, 2 + j], 
        color = (colors[j], 0.5), markersize = 2) for j in 1:3]
    [abline!(axs4[j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    for j in 1:3
        data = DataFrame(x = df_rg_bi[:, 2 + j], y = df_rg_mul[:, 2 + j])
        model = lm(@formula(y ~ x), data)
        abline!(axs4[j], coef(model)[1], coef(model)[2], color = :gold, linewidth = 0.5)
    end
    [hidedecorations!(axs4[j], ticks = false, ticklabels = false) for j in 1:3]
    maxs = [maximum(vcat(df_rg_bi[:, 2 + j], df_rg_mul[:, 2 + j])) for j in 1:3]
    [ylims!(axs4[j], -maxs[j] * 6 / 5, maxs[j] * 6 / 5) for j in 1:3]
    [xlims!(axs4[j], -maxs[j] * 6 / 5, maxs[j] * 6 / 5) for j in 1:3]
    rowsize!(g4, 1, Aspect(3, 1))
    Label(g4[0, 1:3], text = "Comparison of variance components estimates", textsize = 8, tellwidth = false)
    Label(g4[2, 1, Top()], "Cis-genetic covariance", textsize = 6)
    Label(g4[2, 2, Top()], "Trans-genetic covariance", textsize = 6)
    Label(g4[2, 3, Top()], "Residual covariance", textsize = 6)
    Label(g4[1, 1:3, Left()], "e", textsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    Label(g4[3, 1:3], text = "Pairwise bivariate estimates", textsize = 6)
    Label(g4[2, 0], text = "Multivariate estimates", textsize = 6, rotation = π / 2, tellheight = false)
    colgap!(g4, 2.5)
    rowgap!(g4, 0)
    @info "Comparison of rg estimates between bivariate and multivariate models"
    axs5 = [Axis(g5[1, j]) for j in 1:3]
    [scatter!(axs5[j], df_rg_bi[:, 8 + j], df_rg_mul[:, 8 + j], 
        color = (colors[j], 0.5), markersize = 2) for j in 1:3]
    [abline!(axs5[j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    for j in 1:3
        data = DataFrame(x = df_rg_bi[:, 8 + j], y = df_rg_mul[:, 8 + j])
        model = lm(@formula(y ~ x), data)
        abline!(axs5[j], coef(model)[1], coef(model)[2], color = :gold, linewidth = 0.5)
    end
    [hidedecorations!(axs5[j], ticks = false, ticklabels = false) for j in 1:3]
    [ylims!(axs5[j], -1.05, 1.05) for j in 1:3]
    [xlims!(axs5[j], -1.05, 1.05) for j in 1:3]
    rowsize!(g5, 1, Aspect(3, 1))
    Label(g5[0, 1:3], text = "Comparison of rg estimates", textsize = 8, tellwidth = false)
    Label(g5[2, 1, Top()], "rg,cis", textsize = 6)
    Label(g5[2, 2, Top()], "rg,trans", textsize = 6)
    Label(g5[2, 3, Top()], "re", textsize = 6)
    Label(g5[1, 1:3, Left()], "f", textsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    Label(g5[3, 1:3], text = "Pairwise bivariate estimates", textsize = 6)
    Label(g5[2, 0], text = "Multivariate estimates", textsize = 6, rotation = π / 2, tellheight = false)
    colgap!(g5, 2.5)
    rowgap!(g5, 0)
    @info "Comparison of genetic and phenotypic correlation"
    axs6 = [Axis(g6[1, j]) for j in 1:3]
    [scatter!(axs6[j], df_rg.rp, df_rg[:, 7 + j], 
        color = (colors[j], 0.5), markersize = 2) for j in 1:3]
    [abline!(axs6[j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    for j in 1:3
        data = DataFrame(x = df_rg.rp, y = df_rg[:, 7 + j])
        model = lm(@formula(y ~ x), data)
        abline!(axs6[j], coef(model)[1], coef(model)[2], color = :gold, linewidth = 0.5)
    end
    [hidedecorations!(axs6[j], ticks = false, ticklabels = false) for j in 1:3]
    [ylims!(axs6[j], -1.05, 1.05) for j in 1:3]
    [xlims!(axs6[j], -1.05, 1.05) for j in 1:3]
    rowsize!(g6, 1, Aspect(3, 1))
    Label(g6[0, 1:3], text = "Comparison with phenotypic correlation", textsize = 8, tellwidth = false)
    Label(g6[2, 1, Top()], "rg,cis", textsize = 6)
    Label(g6[2, 2, Top()], "rg,trans", textsize = 6)
    Label(g6[2, 3, Top()], "re", textsize = 6)
    Label(g6[1, 1:3, Left()], "f", textsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    Label(g6[3, 1:3], text = "Phenotypic correlation", textsize = 6)
    Label(g6[2, 0], text = "Multivariate estimates", textsize = 6, rotation = π / 2, tellheight = false)
    colgap!(g6, 2.5)
    rowgap!(g6, 0)
    # resize_to_layout!(f)
    save("figure2.pdf", f, pt_per_unit = 1)
end

let
    @info "Comparing ML vs REML estimates"
    @info "Plotting genes"
    f = Figure(resolution = (306, 460))
    axs_title = [Axis(f[1, j]) for j in 1:3]
    [hidespines!(axs_title[j]) for j in 1:3]
    [hidedecorations!(axs_title[j]) for j in 1:3]
    Label(f[1, 1:3, Top()], "Comparison of variance components estimates", textsize = 8)
    rowsize!(f.layout, 1, 0.1)
    Label(f[1, 1, Bottom()], "Cis-genetic variance", textsize = 6)
    Label(f[1, 2, Bottom()], "Trans-genetic variance", textsize = 6)
    Label(f[1, 3, Bottom()], "Residual variance", textsize = 6)
    axs = [Axis(f[i, j]) for i in 2:3, j in 1:3]
    [errorbars!(axs[1, j], resultsg[:, 4 + j], resultsg[:, 16 + j], resultsg[:, 7 + j], 
        linewidth = 0.25, color = ("#9658B2", 0.25)) for j in 1:3]
    [errorbars!(axs[1, j], resultsg[:, 4 + j], resultsg[:, 16 + j], resultsg[:, 19 + j], 
        linewidth = 0.25, color = ("#9658B2", 0.25), direction = :x) for j in 1:3]
    [scatter!(axs[1, j], resultsg[:, 4 + j], resultsg[:, 16 + j], color = ("#4062D8", 0.75), markersize = 2) for j in 1:3]
    [abline!(axs[1, j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    rowsize!(f.layout, 2, Aspect(3, 1))
    maxsg = [maximum(vcat(resultsg[:, 4 + j], resultsg[:, 16 + j])) for j in 1:3]
    [ylims!(axs[1, j], 0, maxsg[j] * 5 / 4) for j in 1:3]
    [xlims!(axs[1, j], 0, maxsg[j] * 5 / 4) for j in 1:3]
    [hidedecorations!(axs[1, j], ticks = false, ticklabels = false) for j in 1:3]
    @info "Plotting isoforms"
    [errorbars!(axs[2, j], resultsi[:, 4 + j], resultsi[:, 16 + j], resultsi[:, 7 + j], 
        linewidth = 0.25, color = ("#9658B2", 0.25)) for j in 1:3]
    [errorbars!(axs[2, j], resultsi[:, 4 + j], resultsi[:, 16 + j], resultsi[:, 19 + j], 
        linewidth = 0.25, color = ("#9658B2", 0.25), direction = :x) for j in 1:3]
    [scatter!(axs[2, j], resultsi[:, 4 + j], resultsi[:, 16 + j], color = ("#4062D8", 0.75), markersize = 2) for j in 1:3]
    [abline!(axs[2, j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    rowsize!(f.layout, 3, Aspect(3, 1))
    maxsi = [maximum(vcat(resultsi[:, 4 + j], resultsi[:, 16 + j])) for j in 1:3]
    [ylims!(axs[2, j], 0, maxsi[j] * 5 / 4) for j in 1:3]
    [xlims!(axs[2, j], 0, maxsi[j] * 5 / 4) for j in 1:3]
    [hidedecorations!(axs[2, j], ticks = false, ticklabels = false) for j in 1:3]
    Label(f[4, 1:3], text = "REML estimates", textsize = 6)
    Label(f[2:3, 0], text = "ML estimates", textsize = 6, rotation = pi / 2, tellheight = false)
    Box(f[2, 5], color = :gray90)
    Label(f[2, 5], "Genes", tellheight = false, textsize = 6, rotation = -π / 2, padding = (0, 0, 3, 3))
    Box(f[3, 5], color = :gray90)
    Label(f[3, 5], "Isoforms", tellheight = false, textsize = 6, rotation = -π / 2, padding = (0, 0, 3, 3))
    colgap!(f.layout, 5)
    rowgap!(f.layout, 5)
    colgap!(f.layout, 4, 2)
    rowgap!(f.layout, 1, 1)
    resize_to_layout!(f)
    save("remlvsmle_vc.png", f, px_per_unit = 4)
end

let
    @info "Comparing ML vs REML estimates"
    @info "Plotting genes"
    f = Figure(resolution = (306, 460))
    axs_title = [Axis(f[1, j]) for j in 1:3]
    [hidespines!(axs_title[j]) for j in 1:3]
    [hidedecorations!(axs_title[j]) for j in 1:3]
    Label(f[1, 1:3, Top()], "Comparison of heritability estimates", textsize = 8)
    rowsize!(f.layout, 1, 0.1)
    Label(f[1, 1, Bottom()], "Cis-heritability", textsize = 6)
    Label(f[1, 2, Bottom()], "Trans-heritability", textsize = 6)
    Label(f[1, 3, Bottom()], "Total heritability", textsize = 6)
    axs = [Axis(f[i, j]) for i in 2:3, j in 1:3]
    [errorbars!(axs[1, j], resultsg[:, 10 + j], resultsg[:, 22 + j], resultsg[:, 13 + j], 
        linewidth = 0.25, color = ("#9658B2", 0.25)) for j in 1:3]
    [errorbars!(axs[1, j], resultsg[:, 10 + j], resultsg[:, 22 + j], resultsg[:, 25 + j], 
        linewidth = 0.25, color = ("#9658B2", 0.25), direction = :x) for j in 1:3]
    [scatter!(axs[1, j], resultsg[:, 10 + j], resultsg[:, 22 + j], color = ("#4062D8", 0.75), markersize = 2) for j in 1:3]
    [abline!(axs[1, j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    rowsize!(f.layout, 2, Aspect(3, 1))
    maxsg = [maximum(vcat(resultsg[:, 10 + j], resultsg[:, 22 + j])) for j in 1:3]
    [ylims!(axs[1, j], 0, maxsg[j] * 5 / 4) for j in 1:3]
    [xlims!(axs[1, j], 0, maxsg[j] * 5 / 4) for j in 1:3]
    [hidedecorations!(axs[1, j], ticks = false, ticklabels = false) for j in 1:3]
    @info "Plotting isoforms"
    [errorbars!(axs[2, j], resultsi[:, 10 + j], resultsi[:, 22 + j], resultsi[:, 13 + j], 
        linewidth = 0.25, color = ("#9658B2", 0.25)) for j in 1:3]
    [errorbars!(axs[2, j], resultsi[:, 10 + j], resultsi[:, 22 + j], resultsi[:, 25 + j], 
        linewidth = 0.25, color = ("#9658B2", 0.25), direction = :x) for j in 1:3]
    [scatter!(axs[2, j], resultsi[:, 10 + j], resultsi[:, 22 + j], color = ("#4062D8", 0.75), markersize = 2) for j in 1:3]
    [abline!(axs[2, j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    rowsize!(f.layout, 3, Aspect(3, 1))
    maxsi = [maximum(vcat(resultsi[:, 10 + j], resultsi[:, 22 + j])) for j in 1:3]
    [ylims!(axs[2, j], 0, maxsi[j] * 5 / 4) for j in 1:3]
    [xlims!(axs[2, j], 0, maxsi[j] * 5 / 4) for j in 1:3]
    [hidedecorations!(axs[2, j], ticks = false, ticklabels = false) for j in 1:3]
    Label(f[4, 1:3], text = "REML estimates", textsize = 6)
    Label(f[2:3, 0], text = "ML estimates", textsize = 6, rotation = pi / 2, tellheight = false)
    Box(f[2, 5], color = :gray90)
    Label(f[2, 5], "Genes", tellheight = false, textsize = 6, rotation = -π / 2, padding = (0, 0, 3, 3))
    Box(f[3, 5], color = :gray90)
    Label(f[3, 5], "Isoforms", tellheight = false, textsize = 6, rotation = -π / 2, padding = (0, 0, 3, 3))
    colgap!(f.layout, 5)
    rowgap!(f.layout, 5)
    colgap!(f.layout, 4, 2)
    rowgap!(f.layout, 1, 1)
    resize_to_layout!(f)
    save("remlvsmle_h2.png", f, px_per_unit = 4)
end

let
    @info "Comparison of heritability estimates for genes with a single isoform"
    gencode_single = filter(row -> row.feature == "transcript", gencode)
    gencode_single = combine(groupby(gencode_single, :gene_id), nrow => :count)
    filter!(row -> row.count == 1, gencode_single)
    genes_single = intersect(gencode_single.gene_id, resultsg.gene_id)
    genes_single = intersect(genes_single, resultsi.gene_id)
    @assert length(genes_single) == 7246
    resultsg_single = Matrix{Float64}(undef, length(genes_single), 3)
    resultsi_single = Matrix{Float64}(undef, length(genes_single), 3)
    for (i, gene) in enumerate(genes_single)
        cols = ["h2_cis_reml", "h2_trans_reml", "h2_reml"]
        resultsg_single[i, :] .= Array(resultsg[findfirst(isequal(gene), resultsg.gene_id), cols])
        resultsi_single[i, :] .= Array(resultsi[findfirst(isequal(gene), resultsi.gene_id), cols])
    end
    f = Figure(resolution = (306, 460))
    axs_title = [Axis(f[1, j]) for j in 1:3]
    [hidespines!(axs_title[j]) for j in 1:3]
    [hidedecorations!(axs_title[j]) for j in 1:3]
    Label(f[1, 1:3, Top()], "Comparison of heritability estimates for genes with a single isoform", textsize = 8)
    rowsize!(f.layout, 1, 0.1)
    Label(f[1, 1, Bottom()], "Cis-heritability", textsize = 6)
    Label(f[1, 2, Bottom()], "Trans-heritability", textsize = 6)
    Label(f[1, 3, Bottom()], "Total heritability", textsize = 6)
    axs = [Axis(f[2, j]) for j in 1:3]
    [scatter!(axs[j], resultsg_single[:, j], resultsi_single[:, j], 
        color = ("#4062D8", 0.75), markersize = 2) for j in 1:3]
    [abline!(axs[j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    [hidedecorations!(axs[j], ticks = false, ticklabels = false) for j in 1:3]
    rowsize!(f.layout, 2, Aspect(3, 1))
    [xlims!(axs[j], 0, 1) for j in 1:3]
    [ylims!(axs[j], 0, 1) for j in 1:3]
    Label(f[3, 1:3], text = "Gene-level estimates", textsize = 6)
    Label(f[2, 0], text = "Isoform-level estimates", textsize = 6, rotation = pi / 2, tellheight = false)
    colgap!(f.layout, 5)
    rowgap!(f.layout, 1)
    resize_to_layout!(f)
    save("single_isoform_genes.png", f, px_per_unit = 4)
end

let
    @info "Comparison of standard error of h² estimates between univariate and multivariate models"
    f = Figure(resolution = (306, 460))
    axs_title = [Axis(f[1, j]) for j in 1:3]
    [hidespines!(axs_title[j]) for j in 1:3]
    [hidedecorations!(axs_title[j]) for j in 1:3]
    Label(f[1, 1:3, Top()], "Comparison of standard error of heritability estimates", textsize = 8)
    rowsize!(f.layout, 1, 0.1)
    Label(f[1, 1, Bottom()], "Cis-heritability", textsize = 6)
    Label(f[1, 2, Bottom()], "Trans-heritability", textsize = 6)
    Label(f[1, 3, Bottom()], "Total heritability", textsize = 6)
    axs = [Axis(f[2, j]) for j in 1:3]
    [scatter!(axs[j], df_uni[:, 10 + j], df[:, 10 + j], 
        color = ("#4062D8", 0.75), markersize = 2) for j in 1:3]
    [abline!(axs[j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    [hidedecorations!(axs[j], ticks = false, ticklabels = false) for j in 1:3]
    rowsize!(f.layout, 2, Aspect(3, 1))
    maxs = [maximum(vcat(df[:, 10 + j], df_uni[:, 10 + j])) for j in 1:3]
    [xlims!(axs[j], 0, maxs[j] * 6 / 5) for j in 1:3]
    [ylims!(axs[j], 0, maxs[j] * 6 / 5) for j in 1:3]
    Label(f[3, 1:3], text = "Univariate model", textsize = 6)
    Label(f[2, 0], text = "Multivariate model", textsize = 6, rotation = pi / 2, tellheight = false)
    colgap!(f.layout, 5)
    rowgap!(f.layout, 1)
    resize_to_layout!(f)
    save("h2_se_univsmul.png", f, px_per_unit = 4)
end

begin
    @info "Comparison of standard error of variance components estimates between univariate and multivariate models"
    f = Figure(resolution = (306, 460))
    axs_title = [Axis(f[1, j]) for j in 1:3]
    [hidespines!(axs_title[j]) for j in 1:3]
    [hidedecorations!(axs_title[j]) for j in 1:3]
    Label(f[1, 1:3, Top()], "Comparison of standard error of variance components estimates", textsize = 8)
    rowsize!(f.layout, 1, 0.1)
    Label(f[1, 1, Bottom()], "Cis-genetic variance", textsize = 6)
    Label(f[1, 2, Bottom()], "Trans-genetic variance", textsize = 6)
    Label(f[1, 3, Bottom()], "Residual variance", textsize = 6)
    axs = [Axis(f[2, j]) for j in 1:3]
    [scatter!(axs[j], df_uni[:, 4 + j], df[:, 4 + j], 
        color = ("#4062D8", 0.75), markersize = 2) for j in 1:3]
    [abline!(axs[j], 0, 1, color = ("#CB3C33", 1), linewidth = 0.5) for j in 1:3]
    [hidedecorations!(axs[j], ticks = false, ticklabels = false) for j in 1:3]
    rowsize!(f.layout, 2, Aspect(3, 1))
    maxs = [maximum(vcat(df[:, 4 + j], df_uni[:, 4 + j])) for j in 1:3]
    [xlims!(axs[j], 0, maxs[j] * 6 / 5) for j in 1:3]
    [ylims!(axs[j], 0, maxs[j] * 6 / 5) for j in 1:3]
    Label(f[3, 1:3], text = "Univariate model", textsize = 6)
    Label(f[2, 0], text = "Multivariate model", textsize = 6, rotation = pi / 2, tellheight = false)
    colgap!(f.layout, 5)
    rowgap!(f.layout, 1)
    resize_to_layout!(f)
    save("vc_se_univsmul.png", f, px_per_unit = 4)
end

@info "Comparison of standard errors of rg between pairwise bivariate and multivariate models"
# cis, trans, residual separately for covariance and correlation

@info "Comparison between h² and eQTL result"
# overlap with venn diagram

@info "Comparison with only cis effects model"
# For genes and isoforms, h² and rg

@info "Polygenicity of cis effects"
# polygenicity (# of conditionally independent signals) + vcsel + index SNP R²

# https://juliadatascience.io/makie_layouts
