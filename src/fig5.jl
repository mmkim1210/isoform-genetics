using CSV, DataFrames, GeneticsMakie, CairoMakie, ColorSchemes, SnpArrays, GLM, Arrow
using Distributions, CategoricalArrays, GLM, MultivariateStats
using GeneticsMakie, CairoMakie
using LinearAlgebra, Statistics
const GM = GeneticsMakie

CairoMakie.activate!(type = "png")
set_theme!(;
    fontsize = 6, 
    font = "Arial",
    Axis = (
        xticksize = 3,
        xticklabelsize = 6,
        yticksize = 3,
        yticklabelsize = 6,
        xlabelsize = 6,
        ylabelsize = 6,
        xtickalign = 0,
        ytickalign = 0,
        xtickwidth = 0.75,
        ytickwidth = 0.75,
        xgridwidth = 0.25,
        ygridwidth = 0.25,
        spinewidth = 0.75,
        xticklabelpad = 0,
        yticklabelpad = 2
        )
    )

if !isfile(joinpath(@__DIR__, "../data/gencode.v40.annotation.parsed.gtf.arrow"))
    run(`curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz --output ./data/gencode.v40.annotation.gtf.gz`)
    @time gencode = CSV.read(joinpath(@__DIR__, "../data/gencode.v40.annotation.gtf.gz"), DataFrame,
        delim = "\t", skipto = 6, header = ["seqnames", "source", "feature", "start", 
        "end", "score", "strand", "phase", "info"])
    GeneticsMakie.parsegtf!(gencode)
    select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_id, :gene_name, :gene_type, :transcript_id)
    @time Arrow.write(joinpath(@__DIR__, "../data/gencode.v40.annotation.parsed.gtf.arrow"), gencode)
    run(`rm ./data/gencode.v40.annotation.gtf.gz`)
end 
@time gencode = Arrow.Table(joinpath(@__DIR__, "../data/gencode.v40.annotation.parsed.gtf.arrow"))|> DataFrame

if !isfile(joinpath(@__DIR__, "../results/long-read/gffcmp.filtered.gencode.arrow"))
    @time gtf = CSV.read(joinpath(@__DIR__, "../results/long-read/gffcmp.filtered.gencode.gtf.gz"), DataFrame,
        delim = "\t", skipto = 4, header = ["seqnames", "source", "feature", "start", 
        "end", "score", "strand", "phase", "info"]) # hg38 coordinates
    gtf.seqnames = replace.(gtf.seqnames, "chr" => "")
    for col in ["transcript_id", "gene_id", "gene_name", "cmp_ref"]
        storage = Vector(undef, size(gtf, 1))
        for (n, m) in enumerate(match.(Regex("$(col) \"(.*?)\";"), gtf.info))
            if isnothing(m)
                storage[n] = missing
            else
                storage[n] = getindex(m.captures, 1)
            end
        end
        gtf[:, col] = storage
    end
    for col in ["num_samples"]
        storage = Vector(undef, size(gtf, 1))
        for (n, m) in enumerate(match.(Regex("$(col) (.*?);"), gtf.info))
            if isnothing(m)
                storage[n] = missing
            else
                storage[n] = getindex(m.captures, 1)
            end
        end
        gtf[:, col] = storage
    end
    for col in ["cmp_ref_gene"]
        storage = Vector(undef, size(gtf, 1))
        for (n, m) in enumerate(match.(Regex("$(col) \"(.*?)\""), gtf.info))
            if isnothing(m)
                storage[n] = missing
            else
                storage[n] = getindex(m.captures, 1)
            end
        end
        gtf[:, col] = storage
    end
    for i in 1:nrow(gtf)
        if ismissing(gtf.gene_name[i])
            gtf.gene_name[i] = gtf.gene_name[i - 1]
            gtf.cmp_ref[i] = gtf.cmp_ref[i - 1]
            gtf.cmp_ref_gene[i] = gtf.cmp_ref_gene[i - 1]
        end
    end
    # gtf.cmp_ref_gene = [getindex(i, 1) for i in split.(gtf.cmp_ref_gene, ".")]
    # gtf.cmp_ref = [getindex(i, 1) for i in split.(gtf.cmp_ref, ".")]
    @time Arrow.write(joinpath(@__DIR__, "../results/long-read/gffcmp.filtered.gencode.arrow"), gtf)
end
gtf = Arrow.Table(joinpath(@__DIR__, "../results/long-read/gffcmp.filtered.gencode.arrow")) |> DataFrame

xrn2 = gtf[findall(isequal("XRN2"), gtf.gene_name), :]
xrn2[findall(isequal("TCONS_00525094"), xrn2.transcript_id), :] # gencode transcript
xrn2[findall(isequal("TCONS_00530677"), xrn2.transcript_id), :] # isoform of interest

# storage = xrn2[findall(isequal("TCONS_00530677"), xrn2.transcript_id), :]
# filter!(row -> row.feature == "exon", storage)
# storage2 = xrn2[findall(isequal("TCONS_00525094"), xrn2.transcript_id), :]
# filter!(row -> row.feature == "exon", storage2)
# for i in 1:nrow(storage)
#     ind = findfirst(isequal(storage.start[i]), storage2.start)
#     if storage.end[i] == storage2.end[ind]
#         println("exon $i has a match")
#     else
#         println("exon $i does not have a match")
#     end
# end

function subsetref(ref::SnpData, chr::AbstractString, range1::Real, range2::Real, path::AbstractString)
    SnpArrays.filter(ref, trues(size(ref)[1]), GeneticsMakie.findlocus(ref, chr, range1, range2); des = path)
    SnpData(path)
end

issig(P::AbstractVector; p = 5e-8) = any(P .< p)
issig(df::DataFrame; p = 5e-8) = issig(df.P; p = p)

geno = SnpData(joinpath(@__DIR__, "../data/fetal/concat.all.eur.filtered"))
chr, start, stop = "20", minimum(xrn2.start), maximum(xrn2.end)
geno = subsetref(geno, chr, start - 1e6, stop + 1e6, "data/fetal/obrien.filtered")

expr = CSV.read("data/fetal/obrien-expr-strict.tsv", DataFrame, delim = "\t")
meta = CSV.read("data/fetal/obrien-metadata.tsv", DataFrame, delim = "\t")
row_meta = [findfirst(isequal(parse(Int, id)), meta.Subject) for id in names(expr)[2:end]]
meta = meta[row_meta, :]

seq = meta[:, [4; 5; 7:12; 14:88]]
storage = Matrix(seq)
# storage = (storage .- mean(storage, dims = 1))
storage = (storage .- mean(storage, dims = 1)) ./ std(storage, dims = 1)
p = fit(PCA, storage', maxoutdim = 50)
pca = MultivariateStats.transform(p, storage')
principalvars(p) / sum(principalvars(p))

covariate = hcat(
    DataFrame(permutedims(pca), :auto),
    select(meta, [:Age, :Sex, :gPC1, :gPC2, :gPC3, :gPC4, :RIN]), 
    )
[sum(ismissing.(col)) for col in eachcol(covariate)]

row_snp = [findfirst(isequal(id), geno.person_info.iid) for id in names(expr)[2:end]]
snparray, snpinfo = (convert(Matrix{Float64}, geno.snparray[row_snp, :]), geno.snp_info)
isoforms = ["TCONS_00530677", "TCONS_00525094", "TCONS_00530678", "TCONS_00530699", "TCONS_00530705",
    "TCONS_00530665", "TCONS_00530679", "TCONS_00530701", "TCONS_00530706"]

ind = Int[]
for i in 1:length(isoforms)
    push!(ind, findfirst(isequal(isoforms[i]), expr.transcript_id))
end
g = permutedims(sum(Array{Float64}(expr[ind, 2:end]), dims = 1))

qtls = Vector{DataFrame}(undef, length(isoforms) + 1)
for i in 1:(length(isoforms) + 1)
    if i == 1
        @info "Working on gene expression"
        covariate.e .= g
    else
        @info "Working on $(isoforms[i - 1])"
        covariate.e = Array{Float64}(expr[findfirst(isequal(isoforms[i - 1]), expr.transcript_id), 2:end])
        # covariate.e = log2.(covariate.e .+ 1)
    end
    qtl = DataFrame(
        SNP = String[], BETA = Float64[], SE = Float64[], Z = Float64[], 
        P = Float64[], BP = Int64[], A1 = String[], A2 = String[], ID = String[]
        )
    for j in 1:size(snparray, 2)
        any(isnan.(snparray[:, j])) ? continue : nothing
        covariate.g =  snparray[:, j]
        model = lm(term(:e) ~ sum(term.(Symbol.(names(covariate, Not(:e))))), covariate)
        i == 1 ? id = "Gene-level" : id = isoforms[i - 1]
        push!(qtl,
            [snpinfo.snpid[j], coef(model)[end], coeftable(model).cols[2][end], 
            coeftable(model).cols[3][end], coeftable(model).cols[4][end],
            snpinfo.position[j], snpinfo.allele1[j],
            snpinfo.allele2[j], id]
            )
    end
    qtl.CHR .= snpinfo.chromosome[1]
    qtls[i] = qtl
end

num_samples = Int[]
for isoform in unique(xrn2.transcript_id)
    push!(num_samples, parse(Int, xrn2[findfirst(isequal(isoform), xrn2.transcript_id), :num_samples]))
end
colors = ["#4062D8", "#389826", "#9658B2", "#CB3C33"]
cs = []
for i in eachindex(num_samples)
    push!(cs, colors[findfirst(isequal(num_samples[i]), unique(num_samples))])
end

begin
    f = Figure(resolution = (530, 792))
    @info "Plotting isoforms"
    g1 = f[1, 1] = GridLayout()
    g2 = f[1, 2] = GridLayout()
    ax = Axis(g1[1, 1])
    rs, chr, range1, range2 = GM.plotisoforms!(ax, "XRN2", xrn2; height = 0.1,
        orderby = ["TCONS_00530677", "TCONS_00525094", "TCONS_00530678", "TCONS_00530699", "TCONS_00530705",
            "TCONS_00530665", "TCONS_00530679", "TCONS_00530701", "TCONS_00530706"],
        highlight = (unique(xrn2.transcript_id), cs), text = :left)
    # vspan!(ax, start, stop; color = (:gray, 0.2))
    GM.labelgenome(g1[1, 1, Bottom()], chr, range1, range2)
    Legend(g1[1, 1], [PolyElement(color = colors[i], strokecolor = :transparent) for i in 1:4], string.(unique(num_samples)),
        tellwidth = false, tellheight = false, rowgap = 0, halign = :left, valign = :bottom,
        framevisible = false, patchsize = (3, 3), strokewidth = 0.1, padding = (10, 3, 3, 3))
    Label(g1[1, 1, Top()], "XRN2 isoforms from long-read sequencing", textsize = 8)
    m = length(unique(xrn2.transcript_id))
    ax2 = Axis(g1[1, 1])
    ylims!(ax2, 0.875 - (m - 1) * 0.125, 1.05)
    ax2.yticks = ([0.95 - (j - 1) * 0.125 for j in 1:m], [["← Isoform of interest", "← In GENCODE v40"]; fill("", m - 2)])
    hidexdecorations!(ax2)
    hideydecorations!(ax2, ticklabels = false)
    hidespines!(ax2)
    ax2.yticklabelsize = 4
    ax2.yaxisposition = :right
    rowsize!(g1, 1, rs)
    ax3 = Axis(g1[2, 1])
    hidedecorations!(ax3)
    hidespines!(ax3)
    rowgap!(g1, 1.5)
    @info "Plotting isoQTL results"
    d = length(qtls)
    axs = [Axis(g2[i, 1]) for i in 1:(d + 1)]
    range1, range2 = start - 0.5e6, stop + 0.5e6
    for i in 1:d
        if issig(qtls[i]; p = 1e-7)
            GeneticsMakie.plotlocus!(axs[i], chr, range1, range2, filter(row -> !isnan(row.P), qtls[i]); ld = geno)
        else
            GeneticsMakie.plotlocus!(axs[i], chr, range1, range2, filter(row -> !isnan(row.P), qtls[i]))
        end
        rowsize!(g2, i, 30)
        Label(g2[i, 1, Top()], "$(qtls[i].ID[1])", textsize = 6, halign = :left, padding = (7.5, 7.5, -7.5, 0))
    end
    rs = GM.plotgenes!(axs[end], chr, range1, range2, gencode; height = 0.1)
    rowsize!(g2, d + 1, rs)
    GM.labelgenome(g2[d + 1, 1, Bottom()], chr, range1, range2)
    Colorbar(g2[1:d, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
        colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
        labelsize = 6, width = 5, spinewidth = 0.5)
    Label(g2[1:d, 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
    rowgap!(g2, 5)
    colgap!(g2, 5)
    colgap!(f.layout, 0)
    for i in 1:(d + 1)
        vlines!(axs[i], start, color = (:gold, 0.5), linewidth = 0.5)
        vlines!(axs[i], stop, color = (:gold, 0.5), linewidth = 0.5)
    end
    for i in 1:d
        lines!(axs[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
    end
    resize_to_layout!(f)
    save("figure5.pdf", f, pt_per_unit = 1)
end

# liftover ASD gwas results
# repeat for all the other genes nearby
# add Sashimi plot

# q1: GENCODE/gencode.v40.annotation.gtf
# q2: Patowary/Isoform_annotations_pseudoBulk.GRCh38.no_gene_names.gtf
# q3: Chun/ALL_combined_isoSeqccs_demux_flnc.5merge.collapsed.rep_classification.filtered_lite.gtf
# q4: GTEx/flair_filter_transcripts.gtf
# q5: “Mill/Filtered - Primary dataset/HumanCTX_sqantitamafiltered.final.classification.gtf”
# q6: PacBio/Alzheimer_IsoSeq2019.postFilter.gtf
# q7: PacBio/UHR_IsoSeq.gff

# RBFOX1, CACNA1C, GRIN2A, SP4, FAM120A, STAG1, CSMD1, DRD2, KMT2E