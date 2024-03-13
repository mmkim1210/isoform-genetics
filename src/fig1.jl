include(joinpath(@__DIR__, "isoform-genetics.jl"))

uni = CSV.read("results/univariate.tsv", DataFrame)
bi = CSV.read("results/bivariate.tsv", DataFrame)
mul = CSV.read("results/multivariate.tsv", DataFrame)

df = combine(groupby(mul, :gene_name), nrow => :count)
sort!(df, order(:count, rev = true))
df = df[df.count .>= 21, :]
findfirst(isequal("FAM200B"), df.gene_name) # potential candidates: KLHL24, SNX19, PPP2R3C, FIGNL1, PHLDB1, MRPL43

gene = Gene("KLHL24", gencode, expr, expri, covariates, 1e6, geno, "cis")
qtls, mega = runqtl(gene)
df = mul[mul.gene_name .== "KLHL24", :]
isoforms = unique(df.id)
gencode_subset = filter(x -> x.transcript_id in isoforms, gencode)
isoforms_sorted = unique(gencode_subset.transcript_id)
d = length(isoforms)
idx = [findfirst(isequal(isoforms_sorted[i]), gene.expressed_isoforms) for i in 1:d]
qtls_subset = GM.flipalleles(qtls[idx .+ 1])

function parse_mul(result, transcript_id)
    d = length(transcript_id)
    Σ̂ = [Matrix{Float64}(undef, d, d) for _ in 1:3]
    for i in 1:nrow(result)
        ind1 = findfirst(isequal(result.id[i]), transcript_id)
        ind2 = findfirst(isequal(result.pair[i]), transcript_id)
        if ind1 != ind2
            for j in 1:3
                Σ̂[1][ind1, ind2] = Σ̂[1][ind2, ind1] = 4 * result[i, "var/covar_cis_reml"]
                Σ̂[2][ind1, ind2] = Σ̂[2][ind2, ind1] = 4 * result[i, "var/covar_trans_reml"]
                Σ̂[3][ind1, ind2] = Σ̂[3][ind2, ind1] = 4 * result[i, "var/covar_res_reml"]
            end
        else
            for j in 1:3
                Σ̂[1][ind1, ind2] = result[i, "var/covar_cis_reml"]
                Σ̂[2][ind1, ind2] = result[i, "var/covar_trans_reml"]
                Σ̂[3][ind1, ind2] = result[i, "var/covar_res_reml"]
            end
        end
    end
    Σ̂
end
Σ̂ = parse_mul(df, isoforms_sorted)
Σ̂p = Σ̂[1] + Σ̂[2] + Σ̂[3]
Σ̂p[7, 7] = Σ̂p[6, 6] = 2.0
Σ̂ ./= maximum(vec(Σ̂p))
Σ̂p ./= maximum(vec(Σ̂p))

begin
    f = Figure(size = (460, 460))
    g1 = f[1, 1:2] = GridLayout()
    g2 = f[2, 1] = GridLayout()
    g3 = f[2, 2] = GridLayout()
    gs = [g1, g2, g3]
    ax1 = Axis(g1[1, 1])
    ax2 = Axis(g1[1, 2])
    ax3 = Axis(g1[1, 3])
    ax4 = Axis(g2[1, 1:4])
    ax5 = Axis(g2[2, 1])
    ax6 = Axis(g2[2, 2])
    ax7 = Axis(g2[2, 3])
    ax8 = Axis(g2[2, 4])
    ax9 = Axis(g2[3, 1:4])
    ax10 = Axis(g3[1, 1])
    ax11 = Axis(g3[2, 1])
    ax12 = Axis(g3[3, 1])
    ax13 = Axis(g3[4, 1])
    axs = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13]
    @info "Plotting isoforms"
    _, _, range1, range2 = GM.plotisoforms!(axs[2], gene.gene_name, gencode_subset; 
            isoformcolor = "black", height = 0.1, text = :left)
    xlims!(axs[2], range1 - 1e4, range2 + 1e4)
    ylims!(axs[2], 0.875 - (d + 1) * 0.125, 1.175)
    rowsize!(gs[1], 1, 18 * (d + 2) * 0.125 / 0.5)
    colsize!(gs[1], 2, 230)
    snps = filter(x -> range1 - 1e4 < x < range2 + 1e4, gene.snpinfo.position)
    barplot!(axs[2], snps,
        fill(0.875 - (d - 0.75) * 0.125, length(snps)), fillto = 0.875 - d * 0.125,
        color = "#389826")
    hidedecorations!(axs[2])
    Label(gs[1][1, 2, Top()], rich("1. Define ", rich("cis", font = :italic), "- and ", rich("trans", font = :italic), "-SNPs"), fontsize = 8)
    Label(gs[1][1, 2, Left()], "SNPs    Isoforms", fontsize = 6, rotation = π / 2)
    [hidespines!(axs[i]) for i in [1, 3]]
    [hidedecorations!(axs[i]) for i in [1, 3]]
    @info "Plotting variance components model"
    GM.plotrg!(axs[5], Σ̂p, diagonal = true, circle = false)
    GM.plotrg!(axs[6], Σ̂[1], diagonal = true, circle = false)
    GM.plotrg!(axs[7], Σ̂[2], diagonal = true, circle = false)
    GM.plotrg!(axs[8], Σ̂[3], diagonal = true, circle = false)
    [xlims!(axs[i], 0, d) for i in 5:8]
    [ylims!(axs[i], -d, 0) for i in 5:8]
    [hidedecorations!(axs[i]) for i in 5:8]
    Label(gs[2][2, 1, Left()], "Isoforms", fontsize = 6, rotation = π / 2)
    [Label(gs[2][2, i, Bottom()], "Isoforms", fontsize = 6) for i in 1:4]
    Label(gs[2][2, 1, Top()], "Phenotypic covariance", fontsize = 6)
    Label(gs[2][2, 2, Top()], rich(rich("Cis", font = :italic), "-genetic covariance"), fontsize = 6)
    Label(gs[2][2, 3, Top()], rich(rich("Trans", font = :italic), "-genetic covariance"), fontsize = 6)
    Label(gs[2][2, 4, Top()], "Residual covariance", fontsize = 6)
    Label(gs[2][1, 1:4, Top()], "2. Fit multivariate variance components model", fontsize = 8)
    rowgap!(g2, 3)
    colgap!(g2, 12.5)
    rowsize!(gs[2], 2, Aspect(4, 1))
    hidespines!(axs[4])
    hidedecorations!(axs[4])
    rowsize!(gs[2], 1, 0.1)
    hidespines!(axs[9])
    hidedecorations!(axs[9])
    @info "Plotting eQTL analysis"
    Label(gs[3][1, 1, Top()], rich("Run ", rich("cis", font = :italic), "-eQTL analysis"), fontsize = 8)
    rowsize!(gs[3], 1, 0.1)
    rowgap!(gs[3], 1, 3)
    hidespines!(axs[10])
    hidedecorations!(axs[10])
    qtls_subset2 = filter(x -> range1 - 1e5 < x.BP < range2 + 1e5, qtls_subset)
    x = range(0.5, 9.5, length = size(qtls_subset2, 1))
    y = range(0.5, d - 0.5, length = d)
    z = hcat(
        sign.(qtls_subset2.Z_9) .* -log10.(qtls_subset2.P_9),
        sign.(qtls_subset2.Z_8) .* -log10.(qtls_subset2.P_8),
        sign.(qtls_subset2.Z_7) .* -log10.(qtls_subset2.P_7),
        sign.(qtls_subset2.Z_6) .* -log10.(qtls_subset2.P_6),
        sign.(qtls_subset2.Z_5) .* -log10.(qtls_subset2.P_5),
        sign.(qtls_subset2.Z_4) .* -log10.(qtls_subset2.P_4),
        sign.(qtls_subset2.Z_3) .* -log10.(qtls_subset2.P_3),
        sign.(qtls_subset2.Z_2) .* -log10.(qtls_subset2.P_2),
        sign.(qtls_subset2.Z_1) .* -log10.(qtls_subset2.P_1), 
        sign.(qtls_subset2.Z) .* -log10.(qtls_subset2.P)
        )
    cmax = ceil(max(abs.(extrema(z))...)) * 1.25
    heatmap!(axs[11], x, y, z, colorrange = (-cmax, cmax), colormap = :RdBu_11)
    rowsize!(gs[3], 2, 52.5)
    Label(gs[3][2, 1, Top()], "SNP effects", fontsize = 6)
    Label(gs[3][2, 1, Left()], "Isoforms", fontsize = 6, rotation = π / 2)
    # xlims!(axs[11], 0.5 - step(x) / 2, 9.5 + step(x) / 2)
    hidedecorations!(axs[11])
    LD = cor(gene.snparray[:, findall(x -> range1 - 1e5 < x < range2 + 1e5, gene.snpinfo.position)])
    LD .= LD.^2
    GM.plotld!(axs[12], LD; colormap = cgrad(:Greens_9, 9, categorical = true))
    rowsize!(g3, 3, Aspect(1, 0.5))
    rowgap!(g3, 2, 0)
    ylims!(axs[12], 0, 5)
    xlims!(axs[12], 0, 10)
    hidespines!(axs[12])
    hidedecorations!(axs[12])
    hidespines!(axs[13])
    hidedecorations!(axs[13])
    rowgap!(f.layout, 10)
    colgap!(f.layout, 10)
    colsize!(f.layout, 1, 250)
    save("figures/figure1.pdf", f, pt_per_unit = 1)
end

for plink in ["bed", "bim", "fam"]
    rm("data/$(gene.gene_name)-cis.$(plink)")
end