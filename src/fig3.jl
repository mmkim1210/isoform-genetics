include(joinpath(@__DIR__, "isoform-genetics.jl"))

gwas, kgp_raw = loadgwas()

uni = CSV.read("results/univariate.tsv", DataFrame)
bi = CSV.read("results/bivariate.tsv", DataFrame)
mul = CSV.read("results/multivariate.tsv", DataFrame)

df = combine(groupby(mul, :gene_name), nrow => :count)
sort!(df, order(:count, rev = true))
df = df[df.count .>= 28, :]
# genes = df.gene_name[Not(45)]
# genes = genes[6:end]
# findfirst(isequal("CLHC1"), genes)

function calculate_h2_uni(result, gene)
    storage = filter(row -> row.feature == "isoform" && row.gene_id == gene.gene_id && row.p < 0.05, result)
    storage.h2_cis_reml, storage.h2_trans_reml, storage.se_h2_cis_reml, storage.se_h2_trans_reml
end

function calculate_h2_bi(result, gene, heritable_isoforms)
    storage1 = filter(row -> row.gene_id == gene.gene_id, result)
    d = length(heritable_isoforms)
    h2cis = Matrix{Float64}(undef, d - 1, d)
    h2trans = Matrix{Float64}(undef, d - 1, d)
    h2cis_se = Matrix{Float64}(undef, d - 1, d)
    h2trans_se = Matrix{Float64}(undef, d - 1, d)
    for i in 1:d
        storage2 = filter(row -> row.id == heritable_isoforms[i], storage1)
        h2cis[:, i] = storage2.h2_cis_reml
        h2trans[:, i] = storage2.h2_trans_reml
        h2cis_se[:, i] = storage2.se_h2_cis_reml
        h2trans_se[:, i] = storage2.se_h2_trans_reml
    end
    h2cis, h2trans, h2cis_se, h2trans_se
end

function calculate_rg_bi(result, gene, heritable_isoforms)
    storage1 = filter(row -> row.gene_id == gene.gene_id, result)
    d = length(heritable_isoforms)
    rg = [Matrix{Float64}(undef, d, d) for _ in 1:3]
    rg_p = [Matrix{Float64}(undef, d, d) for _ in 1:3]
    for i in 1:d
        for j in (i + 1):d
            ind = findfirst(storage1.id .== heritable_isoforms[i] .&& 
                storage1.pair .== heritable_isoforms[j])
            rg[1][i, j] = rg[1][j, i] = storage1.rg_cis_reml[ind]
            rg[2][i, j] = rg[2][j, i] = storage1.rg_trans_reml[ind]
            rg[3][i, j] = rg[3][j, i] = storage1.re_reml[ind]
            rg_p[1][i, j] = rg_p[1][j, i] = storage1.p_cis[ind]
            rg_p[2][i, j] = rg_p[2][j, i] = storage1.p_trans[ind]
            rg_p[3][i, j] = rg_p[3][j, i] = storage1.p_res[ind]
        end
    end
    rg, rg_p
end

function parse_mul(result, gene, heritable_isoforms)
    storage1 = filter(row -> row.gene_id == gene.gene_id, result)
    d = length(heritable_isoforms)
    h2 = [Vector{Float64}(undef, d) for _ in 1:2]
    h2_se = [Vector{Float64}(undef, d) for _ in 1:2]
    rg = [Matrix{Float64}(undef, d, d) for _ in 1:3]
    rg_p = [Matrix{Float64}(undef, d, d) for _ in 1:3]
    for i in 1:nrow(storage1)
        ind1 = findfirst(isequal(storage1.id[i]), heritable_isoforms)
        ind2 = findfirst(isequal(storage1.pair[i]), heritable_isoforms)
        if ind1 != ind2
            rg[1][ind1, ind2] = rg[1][ind2, ind1] = storage1[i, "h2/rg_cis_reml"]
            rg[2][ind1, ind2] = rg[2][ind2, ind1] = storage1[i, "h2/rg_trans_reml"]
            rg[3][ind1, ind2] = rg[3][ind2, ind1] = storage1[i, "h2/re_reml"]
            rg_p[1][ind1, ind2] = rg_p[1][ind2, ind1] = storage1[i, "p_cis"]
            rg_p[2][ind1, ind2] = rg_p[2][ind2, ind1] = storage1[i, "p_trans"]
            rg_p[3][ind1, ind2] = rg_p[3][ind2, ind1] = storage1[i, "p_res"]
        else
            h2[1][ind1] = storage1[i, "h2/rg_cis_reml"]
            h2[2][ind1] = storage1[i, "h2/rg_trans_reml"]
            h2_se[1][ind1] = storage1[i, "se_h2/rg_cis_reml"]
            h2_se[2][ind1] = storage1[i, "se_h2/rg_trans_reml"]
        end
    end
    h2, h2_se, rg, rg_p
end

# for gene_name in genes
gene_name = "ATP9B"
@info "Working on $(gene_name)" # RBM23, COA1, FAM153A, WDR27, PCM1, AKAP13, SDHAP1, ZNF638, ZNF615
window = 1e6
@time gene = Gene(string.(gene_name), gencode, expr, expri, covariates, 1e6, geno, "cis")
range1 = gene.start - window
range2 = gene.stop + window
storage = filter(row -> row.feature == "gene" && row.gene_id == gene.gene_id, uni)
if !isempty(storage.p)
    storage.p[1] < 0.05 ? gene.gene_heritable = true : nothing
end
qtls, _ = runqtl(gene)
storage = filter(row -> row.feature == "isoform" && row.gene_id == gene.gene_id && row.p < 0.05, uni)
heritable_isoforms = unique(storage.id)
qtls_ind = []
for isoform in heritable_isoforms
    push!(qtls_ind, findfirst(isequal(isoform), gene.expressed_isoforms))
end
qtls_ind .+= 1
pushfirst!(qtls_ind, 1)
qtls_subset = qtls[qtls_ind]
@info "Subsetting GWAS results"
@time dfs = subsetgwas(gwas, gene.chr, range1, range2)
gwas_ind = []
for i in 1:length(dfs)
    if issig(dfs[i])
        push!(gwas_ind, i)
    end
end
@info "Subsetting 1000 Genomes"
@time kgp = subsetref(kgp_raw, gene.chr, range1, range2, joinpath.(@__DIR__, "../data/kgp.filtered"))

txid = [heritable_isoforms; setdiff(gene.expressed_isoforms, heritable_isoforms)]
n = length(heritable_isoforms)
m = length(gene.all_isoforms)
range1 < 0 ? range1 = 0 : nothing
range2 > GM.GRCh37_totlength[gene.chr] ? range2 = GM.GRCh37_totlength[gene.chr] : nothing
titles = [GM.gwas[key].title for key in keys(GM.gwas)]

function pvalues(n, p)
    pvalues = Matrix{Float64}(undef, n^2, 3)
    counter = 1
    for i in 1:n
        for j in 1:n
            j == i ? pvalues[counter, 1] = 1 : pvalues[counter, 1] = p[j, i]
            pvalues[counter, 2] = i - 0.5
            pvalues[counter, 3] = -j + 0.5
            counter += 1
        end
    end
    ind = findall(x -> x < 0.05, pvalues[:, 1]) 
    pvalues[ind, :]
end

begin
    f = Figure(size = (530, 792))
    g1 = f[1, 1] = GridLayout()
    g2 = f[1, 2] = GridLayout()
    @info "Plotting panel a"
    axs1 = [Axis(g1[1, 1:4]; alignmode = Outside()), Axis(g1[2, 1:4]), Axis(g1[3, 1:4])]
    rs, _, range1_iso, range2_iso = GM.plotisoforms!(axs1[1], gene_name, gencode; 
        orderby = heritable_isoforms, isoformcolor = "gray60", height = 0.1, text = :left,
        highlight = (txid, [fill("#CB3C33", n); fill("#4062D8", length(txid) - n)]))
    rowsize!(g1, 1, rs)
    GM.labelgenome(g1[1, 1:4, Bottom()], gene.chr, range1_iso, range2_iso)
    Label(g1[1, 1:4, Top()], "$(gene_name) isoforms", fontsize = 8)
    Label(g1[1, 1:4, TopLeft()], "a", fontsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    ax2 = Axis(g1[1, 1:4])
    ylims!(ax2, 0.875 - (m - 1) * 0.125, 1.05)
    ax2.yticks = ([0.95 - (j - 1) * 0.125 for j in 1:m], ["← Iso" .* string.(1:n); fill("", m - n)])
    hidexdecorations!(ax2)
    hideydecorations!(ax2, ticklabels = false)
    hidespines!(ax2)
    ax2.yticklabelsize = 4
    ax2.yaxisposition = :right
    Legend(g1[1, 2:4], [LineElement(color = "#CB3C33", linestyle = nothing), 
        LineElement(color = "#4062D8", linestyle = nothing),
        LineElement(color = "gray60", linestyle = nothing)],
        ["heritable", "non-heritable, brain-expressed", "non-brain-expressed"],
        tellwidth = false, tellheight = false, rowgap = 0, 
        framevisible = false, halign = :left, valign = :bottom,
        patchsize = (5, 1.5), linewidth = 1, padding = (10, 3, 3, 3))

    @info "Plotting panel b"
    h²cis_uni, h²trans_uni, h2cis_uni_se, h2trans_uni_se = calculate_h2_uni(uni, gene)
    h²cis_bi, h²trans_bi, h²cis_bi_se, h²trans_bi_se = calculate_h2_bi(bi, gene, heritable_isoforms)
    rg_bi, rg_bi_p = calculate_rg_bi(bi, gene, heritable_isoforms)
    h2_mul, h2_mul_se, rg_mul, rg_mul_p = parse_mul(mul, gene, heritable_isoforms)
    xs = range(1, 5, length = n + 1)
    ticks = Vector{Float64}(undef, n)
    for i in 1:n
        ys = [h²cis_uni[i]; h²cis_bi[:, i]; h2_mul[1][i]]
        barplot!(axs1[2], xs, ys, fillto = -1, color = 1:(n + 1), colormap = (:jpurple), gap = 0)
        errorbars!(axs1[2], xs, ys, [h2cis_uni_se[i]; h²cis_bi_se[:, i]; h2_mul_se[1][i]], 
            linewidth = 0.25, color = ("black", 0.8))
        ticks[i] = mean(extrema(xs))
        xs = xs .+ (4 + step(xs) * 4)
    end
    xs = range(1, 5, length = n + 1)
    for i in 1:n
        ys = [h²trans_uni[i]; h²trans_bi[:, i]; h2_mul[2][i]]
        barplot!(axs1[3], xs, ys, fillto = -1, color = 1:(n + 1), colormap = (:jgreen), gap = 0)
        errorbars!(axs1[3], xs, ys, [h2trans_uni_se[i]; h²trans_bi_se[:, i]; h2_mul_se[2][i]], 
            linewidth = 0.25, color = ("black", 0.85))
        xs = xs .+ (4 + step(xs) * 4)
    end
    xlims!(axs1[2], 1 - step(xs) * 4, minimum(xs))
    ylims!(axs1[2], 0, 0.14)
    hidexdecorations!(axs1[2])
    hideydecorations!(axs1[2], ticks = false, ticklabels = false)
    xlims!(axs1[3], 1 - step(xs) * 4, minimum(xs))
    ylims!(axs1[3], 0, 0.95)
    hidexdecorations!(axs1[3], ticklabels = false)
    hideydecorations!(axs1[3], ticks = false, ticklabels = false)
    Label(g1[2, 1:4, Top()], "Comparison of heritability estimates", fontsize = 8)
    axs1[2].yticks = 0:0.05:0.22
    axs1[2].xticks = (ticks, "Iso" .* string.(1:n))
    # ax3 = Axis(g1[2, 1], title = L"Comparison of $h^2_{SNP}$ estimates", titlesize = 8)
    axs1[3].yticks = 0:0.25:0.8
    axs1[3].xticks = (ticks, "Iso" .* string.(1:n))
    rowsize!(g1, 2, 45)
    rowsize!(g1, 3, 45)
    Label(g1[2, 1:4, TopLeft()], "b", fontsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    # colors = cgrad(ColorSchemes.jgreen, 4, categorical = true)
    colors = ["#9658B2", "#389826"]
    Legend(g1[2, 1:4], [PolyElement(color = colors[i], strokecolor = :transparent) for i in 1:2], ["cis", "trans"],
        tellwidth = false, tellheight = false, rowgap = 0, halign = :left, valign = :top,
        framevisible = false, patchsize = (3, 3), strokewidth = 0.1, padding = (10, 3, 3, 3))

    @info "Plotting panel c"
    axs2 = [Axis(g1[i, j]) for i in 4:7, j in 1:3]
    for i in 1:n
        for j in (i + 1):n
            rg_mul[1][i, j] = rg_bi[1][i, j]
            rg_mul[2][i, j] = rg_bi[2][i, j]
            rg_mul[3][i, j] = rg_bi[3][i, j]
            rg_mul_p[1][i, j] = rg_bi_p[1][i, j]
            rg_mul_p[2][i, j] = rg_bi_p[2][i, j]
            rg_mul_p[3][i, j] = rg_bi_p[3][i, j]
        end
    end
    GM.plotrg!(axs2[2, 1], rg_mul[1], diagonal = false, circle = false)
    GM.plotrg!(axs2[2, 2], rg_mul[2], diagonal = false, circle = false)
    GM.plotrg!(axs2[2, 3], rg_mul[3], diagonal = false, circle = false)
    ps = [pvalues(n, rg_mul_p[i]) for i in 1:3]
    [scatter!(axs2[2, i], ps[i][:, 2], ps[i][:, 3], marker = '✳', markersize = 7, color = :black) for i in 1:3]
    [axs2[2, i].xticks = (0.5:(n - 0.5), "Iso" .* string.(1:n)) for i in 1:3]
    [axs2[2, i].yticks = (-0.5:-1:-(n - 0.5), "Iso" .* string.(1:n)) for i in 1:3]
    [axs2[2, i].xticklabelsize = 6 for i in 1:3]
    [axs2[2, i].yticklabelsize = 6 for i in 1:3]
    [axs2[2, i].xticklabelpad = 0 for i in 1:3]
    [axs2[2, i].yticklabelpad = 1 for i in 1:3]
    [axs2[2, i].xaxisposition = :top for i in 1:3]
    [axs2[2, i].yaxisposition = :left for i in 1:3]
    [xlims!(axs2[2, i], 0, n) for i in 1:3]
    [ylims!(axs2[2, i], -n, 0) for i in 1:3]
    [hidedecorations!(axs2[2, i], ticklabels = false) for i in 1:3]
    [Label(g1[5, j, Right()], "Pairwise bivariate", fontsize = 6, rotation = -π / 2) for j in 1:3]
    [Label(g1[5, j, Bottom()], "Multivariate", fontsize = 6) for j in 1:3]
    for i in 1:3
        axs2[2, i].xaxisposition = :top
        if n >= 7
            axs2[2, i].xticklabelrotation = π / 4
        end
    end
    Colorbar(g1[5, 4], limits = (-1, 1), ticks = -1:1:1,
        colormap = :RdBu_10, label = "", ticksize = 0, tickwidth = 0,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
        labelsize = 6, width = 5, spinewidth = 0.5, tellheight = false)
    rowsize!(g1, 5, Aspect(3, 1))
    Label(g1[4, 1:4, TopLeft()], "c", fontsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    Label(g1[4, 1:4, Top()], "Comparison of genetic correlation estimates", fontsize = 8)
    Label(g1[4, 1, Bottom()], rich(rich("Cis", font = :italic), "-SNP effects"), fontsize = 6)
    Label(g1[4, 2, Bottom()], rich(rich("Trans", font = :italic), "-SNP effects"), fontsize = 6)
    Label(g1[4, 3, Bottom()], "Residual effects", fontsize = 6)
    rowsize!(g1, 4, 0.1)
    [hidespines!(axs2[1, j]) for j in 1:3]
    [hidedecorations!(axs2[1, j]) for j in 1:3]

    @info "Plotting panel d"
    rp = calculate_rp(gene, uni)
    for i in 1:n
        for j in (i + 1):n
            rg_mul[1][i, j] = rp[i, j]
            rg_mul[2][i, j] = rp[i, j]
            rg_mul[3][i, j] = rp[i, j]
            rg_mul_p[1][i, j] = 1
            rg_mul_p[2][i, j] = 1
            rg_mul_p[3][i, j] = 1
        end
    end
    GM.plotrg!(axs2[4, 1], rg_mul[1], diagonal = false, circle = false)
    GM.plotrg!(axs2[4, 2], rg_mul[2], diagonal = false, circle = false)
    GM.plotrg!(axs2[4, 3], rg_mul[3], diagonal = false, circle = false)
    ps = [pvalues(n, rg_mul_p[i]) for i in 1:3]
    [scatter!(axs2[4, i], ps[i][:, 2], ps[i][:, 3], marker = '✳', markersize = 7, color = :black) for i in 1:3]
    [axs2[4, i].xticks = (0.5:(n - 0.5), "Iso" .* string.(1:n)) for i in 1:3]
    [axs2[4, i].yticks = (-0.5:-1:-(n - 0.5), "Iso" .* string.(1:n)) for i in 1:3]
    [axs2[4, i].xticklabelsize = 6 for i in 1:3]
    [axs2[4, i].yticklabelsize = 6 for i in 1:3]
    [axs2[4, i].xticklabelpad = 0 for i in 1:3]
    [axs2[4, i].yticklabelpad = 1 for i in 1:3]
    [axs2[4, i].xaxisposition = :top for i in 1:3]
    [axs2[4, i].yaxisposition = :left for i in 1:3]
    [xlims!(axs2[4, i], 0, n) for i in 1:3]
    [ylims!(axs2[4, i], -n, 0) for i in 1:3]
    [hidedecorations!(axs2[4, i], ticklabels = false) for i in 1:3]
    [Label(g1[7, j, Right()], "Phenotypic correlation", fontsize = 6, rotation = -π / 2) for j in 1:3]
    [Label(g1[7, j, Bottom()], "Multivariate", fontsize = 6) for j in 1:3]
    for i in 1:3
        axs2[4, i].xaxisposition = :top
        if n >= 7
            axs2[4, i].xticklabelrotation = π / 4
        end
    end
    Colorbar(g1[7, 4], limits = (-1, 1), ticks = -1:1:1,
        colormap = :RdBu_10, label = "", ticksize = 0, tickwidth = 0,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
        labelsize = 6, width = 5, spinewidth = 0.5, tellheight = false)
    rowsize!(g1, 7, Aspect(3, 1))
    Label(g1[6, 1:4, TopLeft()], "d", fontsize = 12, font = "Arial bold", padding = (0, 5, 0, 0), halign = :right)
    Label(g1[6, 1:4, Top()], "Comparison with phenotypic correlation", fontsize = 8)
    Label(g1[6, 1, Bottom()], rich(rich("Cis", font = :italic), "-SNP effects"), fontsize = 6)
    Label(g1[6, 2, Bottom()], rich(rich("Trans", font = :italic), "-SNP effects"), fontsize = 6)
    Label(g1[6, 3, Bottom()], "Residual effects", fontsize = 6)
    rowsize!(g1, 6, 0.1)
    [hidespines!(axs2[3, j]) for j in 1:3]
    [hidedecorations!(axs2[3, j]) for j in 1:3]
    ax3 = Axis(g1[8, 1:4])
    hidespines!(ax3)
    hidedecorations!(ax3)
    rowgap!(g1, 5)
    colgap!(g1, 5)
    rowgap!(g1, 4, 1)
    rowgap!(g1, 6, 1)
    rowgap!(g1, 2, 2.5)
    @info "Plotting panel e"
    # l = length(gwas_ind) + n + 2
    l = n + 2
    axs3 = [Axis(g2[i, 1]) for i in 1:l]
    qtls_plot = subsetgwas(qtls_subset, gene.chr, gene.start - 0.5e6, gene.stop + 0.5e6)
    for i in 1:(n + 1)
        if issig(qtls_plot[i])
            GM.plotlocus!(axs3[i], gene.chr, gene.start - 0.5e6, gene.stop + 0.5e6, qtls_plot[i]; ld = kgp)
        else
            GM.plotlocus!(axs3[i], gene.chr, gene.start - 0.5e6, gene.stop + 0.5e6, qtls_plot[i])
        end
        rowsize!(g2, i, 30)
        if i == 1
            Label(g2[i, 1, Top()], "$(gene_name) (Gene)", fontsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
        else
            Label(g2[i, 1, Top()], "$(gene_name) (Iso$(i - 1))", fontsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
        end
    end
    # for i in 1:length(gwas_ind)
    #     if issig(dfs[gwas_ind[i]])
    #         GM.plotlocus!(axs3[i + n + 1], gene.chr, range1, range2, dfs[gwas_ind[i]]; ld = kgp)
    #         if dfs[gwas_ind[i]].BP[findmin(dfs[gwas_ind[i]].P)[2]] < (range1 + range2) / 2
    #             Label(g2[i + n + 1, 1, Top()], "$(titles[gwas_ind[i]])", fontsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
    #         else
    #             Label(g2[i + n + 1, 1, Top()], "$(titles[gwas_ind[i]])", fontsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
    #         end    
    #     else
    #         GM.plotlocus!(axs3[i + n + 1], gene.chr, range1, range2, dfs[gwas_ind[i]])
    #         Label(g2[i + n + 1, 1, Top()], "$(titles[gwas_ind[i]])", fontsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
    #     end
    #     rowsize!(g2, i + n + 1, 30)
    # end
    rs = GM.plotgenes!(axs3[end], gene.chr, gene.start - 0.5e6, gene.stop + 0.5e6, gencode; height = 0.1)
    rowsize!(g2, n + 2, rs)
    GM.labelgenome(g2[n + 2, 1, Bottom()], gene.chr, gene.start - 0.5e6, gene.stop + 0.5e6)
    Colorbar(g2[1:(n + 1), 2], limits = (0, 1), ticks = 0:1:1, height = 20,
        colormap = [:gray60, :red2], label = "LD", ticksize = 0, tickwidth = 0,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
        labelsize = 6, width = 5, spinewidth = 0.5)
    Label(g2[1:(n + 1), 0], text = rich("-log", subscript("10"), rich("P", font = :italic)), fontsize = 6, rotation = pi / 2)
    for i in 1:(n + 2)
        vlines!(axs3[i], gene.start, color = (:gold, 0.5), linewidth = 0.5)
        vlines!(axs3[i], gene.stop, color = (:gold, 0.5), linewidth = 0.5)
    end
    for i in 1:(n + 1)
        lines!(axs3[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
    end
    # rowsize!(g2, length(gwas_ind) + n + 2, rs)
    # GM.labelgenome(g2[length(gwas_ind) + n + 2, 1, Bottom()], gene.chr, range1, range2)
    # Colorbar(g2[1:(length(gwas_ind) + n + 1), 2], limits = (0, 1), ticks = 0:1:1, height = 20,
    #     colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
    #     tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
    #     labelsize = 6, width = 5, spinewidth = 0.5)
    # Label(g2[1:(length(gwas_ind) + n + 1), 0], text = "-log[p]", fontsize = 6, rotation = pi / 2)
    # for i in 1:(length(gwas_ind) + n + 2)
    #     vlines!(axs3[i], gene.start, color = (:gold, 0.5), linewidth = 0.5)
    #     vlines!(axs3[i], gene.stop, color = (:gold, 0.5), linewidth = 0.5)
    # end
    # for i in 1:(length(gwas_ind) + n + 1)
    #     lines!(axs3[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
    # end
    rowgap!(g2, 5)
    colgap!(g2, 5)
    colgap!(f.layout, 1, 0)
    resize_to_layout!(f)
    save("figures/figure3.pdf", f, pt_per_unit = 1)
    # save("figure3-$(gene.gene_name).png", f, px_per_unit = 4)
end

for plink in ["bed", "bim", "fam"]
    rm("data/$(gene.gene_name)-cis.$(plink)")
end