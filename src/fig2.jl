include(joinpath(@__DIR__, "isoform-genetics.jl"))

results = DataFrame[]
for model in ["univariate", "bivariate", "multivariate"]
    push!(results, CSV.read("results/$(model).tsv", DataFrame))
end

resultsg = filter(row -> row.feature == "gene", results[1])
resultsi = filter(row -> row.feature == "isoform", results[1])

begin
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
    rowgap!(f.layout, 1, 2)
    resize_to_layout!(f)
    save("remlvsmle_vc.png", f, px_per_unit = 4)
end

begin
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
    rowgap!(f.layout, 1, 2)
    resize_to_layout!(f)
    save("remlvsmle_h2.png", f, px_per_unit = 4)
end

@info "Comparison of heritability estimates for genes with a single isoform"
gencode_single = filter(row -> row.feature == "transcript", gencode)
gencode_single = combine(groupby(gencode_single, :gene_id), nrow => :count)
filter!(row -> row.count == 1, gencode_single)
genes_single = intersect(gencode_single.gene_id, resultsg.gene_id)
genes_single = intersect(genes_single, resultsi.gene_id) # 7,246 genes with a single transcript
resultsg_single = Matrix{Float64}(undef, length(genes_single), 3)
resultsi_single = Matrix{Float64}(undef, length(genes_single), 3)
for (i, gene) in enumerate(genes_single)
    cols = ["h2_cis_reml", "h2_trans_reml", "h2_reml"]
    resultsg_single[i, :] .= Array(resultsg[findfirst(isequal(gene), resultsg.gene_id), cols])
    resultsi_single[i, :] .= Array(resultsi[findfirst(isequal(gene), resultsi.gene_id), cols])
end

begin
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
    rowgap!(f.layout, 2)
    save("single_isoform_genes.png", f, px_per_unit = 4)
end

h²sg = [resultsg[:, "h2_cis_reml"], resultsg[:, "h2_trans_reml"], resultsg[:, "h2_reml"]]
h²si = [resultsi[:, "h2_cis_reml"], resultsi[:, "h2_trans_reml"], resultsi[:, "h2_reml"]]

begin
    f = Figure(resolution = (306, 460))
    @info "Distribution of h² estimates for univariate model"
    axs1 = [Axis(f[1, j]) for j in 1:2]
    for (i, h²) in enumerate(h²sg)
        hist!(axs1[1], h², offset = -i, scale_to = 0.8, strokecolor = (:black, 0.5), strokewidth = 0.25,
            bins = 100)
    end
    for (i, h²) in enumerate(h²si)
        hist!(axs1[2], h², offset = -i, scale_to = 0.8, strokecolor = (:black, 0.5), strokewidth = 0.25,
            bins = 100)
    end
    rowsize!(f.layout, 1, 60)
    # add boxplots below
    colgap!(f.layout, 5)
    save("figure2.png", f, px_per_unit = 4)
end