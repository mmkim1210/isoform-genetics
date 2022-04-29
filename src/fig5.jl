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

@time gencode = CSV.read(joinpath(@__DIR__, "../results/long-read/gffcmp.filtered.gencode.gtf.gz"), DataFrame,
    delim = "\t", skipto = 4, header = ["seqnames", "source", "feature", "start", 
    "end", "score", "strand", "phase", "info"]) # hg38 coordinates
gencode.seqnames = replace.(gencode.seqnames, "chr" => "")
for col in ["transcript_id", "gene_id", "gene_name", "cmp_ref"]
    storage = Vector(undef, size(gencode, 1))
    for (n, m) in enumerate(match.(Regex("$(col) \"(.*?)\";"), gencode.info))
        if isnothing(m)
            storage[n] = missing
        else
            storage[n] = getindex(m.captures, 1)
        end
    end
    gencode[:, col] = storage
end
for col in ["num_samples"]
    storage = Vector(undef, size(gencode, 1))
    for (n, m) in enumerate(match.(Regex("$(col) (.*?);"), gencode.info))
        if isnothing(m)
            storage[n] = missing
        else
            storage[n] = getindex(m.captures, 1)
        end
    end
    gencode[:, col] = storage
end
for col in ["cmp_ref_gene"]
    storage = Vector(undef, size(gencode, 1))
    for (n, m) in enumerate(match.(Regex("$(col) \"(.*?)\""), gencode.info))
        if isnothing(m)
            storage[n] = missing
        else
            storage[n] = getindex(m.captures, 1)
        end
    end
    gencode[:, col] = storage
end
# gencode.cmp_ref_gene = [getindex(i, 1) for i in split.(gencode.cmp_ref_gene, ".")]
# gencode.cmp_ref = [getindex(i, 1) for i in split.(gencode.cmp_ref, ".")]

xrn2 = gencode[findall(isequal("ENSG00000088930.8"), gencode.cmp_ref_gene), :]
xrn2 = gencode[findall(in(xrn2.transcript_id), gencode.transcript_id), :]
xrn2.gene_name .= "XRN2"
xrn2[findall(isequal("TCONS_00530677"), xrn2.transcript_id), :]
nsamples = Int[]
for isoform in unique(xrn2.transcript_id)
    push!(nsamples, parse(Int, xrn2[findfirst(isequal(isoform), xrn2.transcript_id), :num_samples]))
end
@time let
    f = Figure(resolution = (306, 792))
    ax = Axis(f[1, 1])
    rs, chr, range1, range2 = GM.plotisoforms!(ax, "XRN2", xrn2; height = 0.1,
        highlight = (["TCONS_00530677"], ["#CB3C33"]))
    GM.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
    rowsize!(f.layout, 1, rs)
    resize_to_layout!(f)
    save("figs/XRN2-long-read.png", f, px_per_unit = 4)
end

# q1: GENCODE/gencode.v40.annotation.gtf
# q2: Patowary/Isoform_annotations_pseudoBulk.GRCh38.no_gene_names.gtf
# q3: Chun/ALL_combined_isoSeqccs_demux_flnc.5merge.collapsed.rep_classification.filtered_lite.gtf
# q4: GTEx/flair_filter_transcripts.gtf
# q5: “Mill/Filtered - Primary dataset/HumanCTX_sqantitamafiltered.final.classification.gtf”
# q6: PacBio/Alzheimer_IsoSeq2019.postFilter.gtf
# q7: PacBio/UHR_IsoSeq.gff