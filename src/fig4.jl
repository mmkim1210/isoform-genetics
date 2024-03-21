include(joinpath(@__DIR__, "isoform-genetics.jl"))

@info "Load GWAS results."
delete!(GeneticsMakie.gwas, "mdd")
gwas = []
for key in keys(GeneticsMakie.gwas)
    push!(gwas, DataFrame(Arrow.Table(joinpath.("data/gwas/processed/", key * ".tsv.arrow"))))
end

@info "Loading 1000 Genomes"
@time kgp_raw = SnpData(joinpath(@__DIR__, "../data/1kg/kgp.EUR.maf0.05.geno"))

genes_to_focus = ["XRN2", "TBL1XR1", "SYNE1", "SYT1"]
isoforms_to_focus = [
    ["ENST00000430571"], 
    ["ENST00000422066", "ENST00000474363"], 
    ["ENST00000472563", "ENST00000461872", "ENST00000466159"], 
    ["ENST00000549454", "ENST00000457153", "ENST00000552624"]
]
gwas_to_focus = [["asd"], ["scz"], ["bd"], ["ea"]]

genes = []
qtls_all = []
for gene_name in genes_to_focus
    @info "Working on $(gene_name)"
    window = 1e6
    @time gene = Gene(gene_name, gencode, expr, expri, covariates, 1e6, geno, "cis")
    push!(genes, gene)
    @time qtls, _ = runqtl(gene)
    push!(qtls_all, qtls)
end

begin
    f = Figure(resolution = (530, 792))
    g1 = f[1, 1] = GridLayout()
    g2 = f[1, 2] = GridLayout()
    ga = g1[1, 1] = GridLayout()
    gb = g2[1, 1] = GridLayout()
    gc = g1[2, 1] = GridLayout()
    gd = g2[2, 1] = GridLayout()
    gs = [ga, gb, gc, gd]

    window = 0.5e6
    for (i, gene_name) in enumerate(genes_to_focus)
        @info "Working on $(gene_name)"
        range1 = genes[i].start - window
        range2 = genes[i].stop + window
        @info "Subsetting GWAS results"
        @time dfs = subsetgwas(gwas, genes[i].chr, range1, range2)
        @info "Subsetting 1000 Genomes"
        @time kgp = subsetref(kgp_raw, genes[i].chr, range1, range2, joinpath.(@__DIR__, "../data/kgp.filtered"))
        @info "Plotting LocusZoom"
        n = length(isoforms_to_focus[i])
        m = length(gwas_to_focus[i])
        range1 < 0 ? range1 = 0 : nothing
        range2 > GM.GRCh37_totlength[genes[i].chr] ? range2 = GM.GRCh37_totlength[genes[i].chr] : nothing
        titles = [GM.gwas[key].title for key in gwas_to_focus[i]]
        axs = [Axis(gs[i][j, 1]) for j in 1:(n + 3 + m)]
        @info "Plotting QTL"
        ind = []
        for isoform in isoforms_to_focus[i]
            push!(ind, findfirst(isequal(isoform), genes[i].expressed_isoforms))
        end        
        for j in 1:(n + 1)
            if j == 1
                GM.plotlocus!(axs[j], genes[i].chr, range1, range2, qtls_all[i][1])
                Label(gs[i][j, 1, Top()], "$(genes[i].gene_id)", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            else
                GM.plotlocus!(axs[j], genes[i].chr, range1, range2, qtls_all[i][ind[j - 1] + 1]; ld = kgp)
                Label(gs[i][j, 1, Top()], "$(isoforms_to_focus[i][j - 1])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end
            rowsize!(gs[i], j, 30)
        end
        @info "Plotting GWAS"
        ind = findall(in(gwas_to_focus[i]), collect(keys(GeneticsMakie.gwas)))
        for j in 1:m
            GM.plotlocus!(axs[j + n + 1], genes[i].chr, range1, range2, dfs[ind[j]]; ld = kgp)
            Label(gs[i][j + n + 1, 1, Top()], "$(titles[j])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            rowsize!(gs[i], j + n + 1, 30)
        end
        @info "Plotting genes"
        rs = GM.plotgenes!(axs[n + m + 2], genes[i].chr, range1, range2, gencode; height = 0.1)
        rowsize!(gs[i], n + m + 2, rs)
        GM.labelgenome(gs[i][n + m + 2, 1, Bottom()], genes[i].chr, range1, range2)
        @info "Plotting isoforms"
        gencode_subset = filter(x -> x.gene_name == gene_name, gencode)
        gencode_subset = gencode_subset[findall(in(genes[i].expressed_isoforms), gencode_subset.transcript_id), :]
        rs, _, range1_iso, range2_iso = GM.plotisoforms!(axs[n + m + 3], gene_name, gencode_subset; 
            orderby = isoforms_to_focus[i], isoformcolor = "#4062D8", height = 0.1, text = :left,
            highlight = (isoforms_to_focus[i], fill("#CB3C33", length(isoforms_to_focus[i]))))
        rowsize!(gs[i], n + m + 3, rs)
        GM.labelgenome(gs[i][n + m + 3, 1, Bottom()], genes[i].chr, range1_iso, range2_iso)
        for j in 1:(n + 2 + m)
            vlines!(axs[j], genes[i].start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs[j], genes[i].stop, color = (:gold, 0.5), linewidth = 0.5)
        end
        for j in 1:(n + 1 + m)
            lines!(axs[j], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        end
        Colorbar(gs[i][1:(n + m + 1), 2], limits = (0, 1), ticks = 0:1:1, height = 20,
            colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
            tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
            labelsize = 6, width = 5, spinewidth = 0.5)
        Label(gs[i][1:(n + m + 1), 1, Left()], text = "-log[p]", textsize = 6, rotation = π / 2)
        rowgap!(gs[i], 5)
        colgap!(gs[i], 5)
    end
    rowgap!(f.layout, 2)
    colgap!(f.layout, 0)
    # resize_to_layout!(f)
    # save("figure4.svg", f, pt_per_unit = 1)
    save("figure4.pdf", f, pt_per_unit = 1)
end