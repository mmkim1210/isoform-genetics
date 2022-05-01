using CSV, DataFrames, GeneticsMakie, CairoMakie
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

@time gencode = CSV.read(joinpath(@__DIR__, "../results/long-read/gencode.v40.XRN2.gtf"), DataFrame)

begin
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
    # gtf.cmp_ref_gene = [getindex(i, 1) for i in split.(gtf.cmp_ref_gene, ".")]
    # gtf.cmp_ref = [getindex(i, 1) for i in split.(gtf.cmp_ref, ".")]
end

gene_id = gtf.cmp_ref_gene[findfirst(isequal("XRN2"), gtf.gene_name)]
xrn2 = gtf[findall(isequal(gene_id), gtf.cmp_ref_gene), :]
xrn2 = gtf[findall(in(xrn2.transcript_id), gtf.transcript_id), :]
xrn2.gene_name .= "XRN2"

filter(row -> row.feature == "exon", gencode)
xrn2[findall(isequal("TCONS_00525094"), xrn2.transcript_id), :] # gencode transcript
xrn2[findall(isequal("TCONS_00530677"), xrn2.transcript_id), :] # isoform of interest

storage = xrn2[findall(isequal("TCONS_00530677"), xrn2.transcript_id), :]
filter!(row -> row.feature == "exon", storage)
storage2 = xrn2[findall(isequal("TCONS_00525094"), xrn2.transcript_id), :]
filter!(row -> row.feature == "exon", storage2)
for i in 1:nrow(storage)
    ind = findfirst(isequal(storage.start[i]), storage2.start)
    if storage.end[i] == storage2.end[ind]
        println("exon $i has a match")
    else
        println("exon $i does not have a match")
    end
end

nsamples = Int[]
for isoform in unique(xrn2.transcript_id)
    push!(nsamples, parse(Int, xrn2[findfirst(isequal(isoform), xrn2.transcript_id), :num_samples]))
end
colors = ["#4062D8", "#389826", "#9658B2", "#CB3C33"]
cs = []
for i in eachindex(nsamples)
    push!(cs, colors[findfirst(isequal(nsamples[i]), unique(nsamples))])
end
@time let
    f = Figure(resolution = (306, 792))
    ax = Axis(f[1, 1])
    rs, chr, range1, range2 = GM.plotisoforms!(ax, "XRN2", xrn2; height = 0.1,
        orderby = ["TCONS_00530677", "TCONS_00525094"],
        highlight = (unique(xrn2.transcript_id), cs), text = :left)
    GM.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
    Legend(f[1, 1], [PolyElement(color = colors[i], strokecolor = :transparent) for i in 1:4], string.(unique(nsamples)),
        tellwidth = false, tellheight = false, rowgap = 0, halign = :left, valign = :bottom,
        framevisible = false, patchsize = (3, 3), strokewidth = 0.1, padding = (10, 3, 3, 3))
    rowsize!(f.layout, 1, rs)
    resize_to_layout!(f)
    Label(f[1, 1, Top()], "XRN2 isoforms from long-read sequencing", textsize = 8)
    m = length(unique(xrn2.transcript_id))
    ax2 = Axis(f[1, 1])
    ylims!(ax2, 0.875 - (m - 1) * 0.125, 1.05)
    ax2.yticks = ([0.95 - (j - 1) * 0.125 for j in 1:m], [["← Isoform of interest", "← In GENCODE"]; fill("", m - 2)])
    hidexdecorations!(ax2)
    hideydecorations!(ax2, ticklabels = false)
    hidespines!(ax2)
    ax2.yticklabelsize = 4
    ax2.yaxisposition = :right
    save("figs/$(gene_id)-long-read.png", f, px_per_unit = 4)
end

# q1: GENCODE/gencode.v40.annotation.gtf
# q2: Patowary/Isoform_annotations_pseudoBulk.GRCh38.no_gene_names.gtf
# q3: Chun/ALL_combined_isoSeqccs_demux_flnc.5merge.collapsed.rep_classification.filtered_lite.gtf
# q4: GTEx/flair_filter_transcripts.gtf
# q5: “Mill/Filtered - Primary dataset/HumanCTX_sqantitamafiltered.final.classification.gtf”
# q6: PacBio/Alzheimer_IsoSeq2019.postFilter.gtf
# q7: PacBio/UHR_IsoSeq.gff

# RBFOX1, CACNA1C, GRIN2A, SP4, FAM120A, STAG1, CSMD1, DRD2, KMT2E