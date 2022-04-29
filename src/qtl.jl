include(joinpath(@__DIR__, "isoform-genetics.jl"))

gwas, kgp_raw = loadgwas()

# genes = ["TBL1XR1", "SYT1", "SF3B1", "NFIB", "PTCH1", "GNAI1", "GATM", "ANK3", "WAC", "USP7", "CSDE1", "CAMTA1", "ZBTB20"]
# genes = ["TMEM106B", "SOX5", "PAX6", "MYT1L", "NRXN3", "RSRC1", "FAM104A", "HERC1", "KDM3B", "FBXO11", "LARP4B", "ANKS1B", "ASXL3", "CTTNBP2", "ZBTB7A"]
# genes = ["PTEN", "STARD9", "GIGYF2", "AFF3", "TRPM3", "CTNNA2", "UIMC1", "TSHZ3"]

genes = ["CRYBA1"]

for gene_name in genes
    !isfile("figs/$(gene_name)-qtl.png") || continue
    @info "Working on $(gene_name)"
    window = 1e6
    @time gene = Gene(gene_name, gencode, expr, expri, cov, 1e6, geno, nothing)
    range1 = gene.start - window
    range2 = gene.stop + window
    # @time qtls, _ = runqtl(gene)
    @info "Subsetting GWAS results"
    @time dfs = subsetgwas(gwas, gene.chr, range1, range2)
    @info "Subsetting 1000 Genomes"
    @time kgp = subsetref(kgp_raw, gene.chr, range1, range2, joinpath.(@__DIR__, "../data/kgp.filtered"))
    @info "Plotting LocusZoom"
    @time begin
        # n = length(gene.expressed_isoforms)
        range1 < 0 ? range1 = 0 : nothing
        range2 > GM.GRCh37_totlength[gene.chr] ? range2 = GM.GRCh37_totlength[gene.chr] : nothing
        titles = [GM.gwas[key].title for key in keys(GM.gwas)]

        f = Figure(resolution = (306, 792))
        gwas_ind = []
        for i in 1:length(gwas)
            if issig(dfs[i])
                push!(gwas_ind, i)
            end
        end
        axs = [Axis(f[i, 1]) for i in 1:(1 + length(gwas_ind))]
        # @info "Plotting QTL"
        # for i in 1:(n + 1)
        #     if issig(qtls[i])
        #         GM.plotlocus!(axs[i + 1], gene.chr, range1, range2, qtls[i]; ld = kgp)
        #     else
        #         GM.plotlocus!(axs[i + 1], gene.chr, range1, range2, qtls[i])
        #     end
        #     rowsize!(f.layout, i + 1, 30)
        #     if i == 1
        #         Label(f[i + 1, 1, Top()], "$(gene_name) (Gene)", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
        #     else
        #         Label(f[i + 1, 1, Top()], "$(gene_name) ($(gene.expressed_isoforms[i - 1]))", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
        #     end
        # end
        @info "Plotting GWAS"
        for (i, ind) in enumerate(gwas_ind)    
            GM.plotlocus!(axs[i], gene.chr, range1, range2, dfs[ind]; ld = kgp)
            if dfs[ind].BP[argmin(dfs[ind].P)] < (range1 + range2) / 2
                Label(f[i, 1, Top()], "$(titles[ind])", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
            else
                Label(f[i, 1, Top()], "$(titles[ind])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end    
            rowsize!(f.layout, i , 30)
        end
        rs = GM.plotgenes!(axs[length(gwas_ind) + 1], gene.chr, range1, range2, gencode; height = 0.1)
        rowsize!(f.layout, length(gwas_ind) + 1, rs)
        GM.labelgenome(f[1 + length(gwas_ind), 1, Bottom()], gene.chr, range1, range2)
        Colorbar(f[1:length(gwas_ind), 2], limits = (0, 1), ticks = 0:1:1, height = 20,
            colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
            tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
            labelsize = 6, width = 5, spinewidth = 0.5)
        Label(f[1:length(gwas_ind), 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
        for i in 1:(1 + length(gwas_ind))
            vlines!(axs[i], gene.start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs[i], gene.stop, color = (:gold, 0.5), linewidth = 0.5)
        end
        for i in 1:length(gwas_ind)
            lines!(axs[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        end
        rowgap!(f.layout, 5)
        colgap!(f.layout, 5)
        colgap!(f.layout, 1, 5)
        resize_to_layout!(f)
        save("figs/$(gene_name)-qtl.png", f, px_per_unit = 4)
    end
end

# SLC12A8 distinct loci + COL28A1 + DDX6 + PAK1, STAU1 (DTU) DNAJA3 + ELAC2
# "XRN2", "SYNE1", "TBL1XR1", "SYT1" + WAC CTNNB1 GIGYF2
# df = qtls[findfirst(isequal("ENST00000549454"), gene.expressed_isoforms) + 1]
# df = innerjoin(df, dfs[gwas_ind[3]]; on = :BP, makeunique = true)
# scatter(-log10.(df.P), -log10.(df.P_1))