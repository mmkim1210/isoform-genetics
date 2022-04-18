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


@time gencode = CSV.read(joinpath(@__DIR__, "../results/XRN2/gffcmp.XRN2.filtered.gtf"), DataFrame,
    delim = "\t", skipto = 5, header = ["seqnames", "source", "feature", "start", 
    "end", "score", "strand", "phase", "info"])
gencode.gene_id .= "ENSG00000088930"
gencode.gene_name .= "XRN2"
storage = Vector(undef, size(gencode, 1))
for (n, m) in enumerate(match.(Regex("transcript_id \"(.*?)\";"), gencode.info))
    if isnothing(m)
        storage[n] = gencode.gene_id[n]
    else
        storage[n] = getindex(m.captures, 1)
    end
end
gencode[:, :transcript_id] = storage
gencode.seqnames = replace.(gencode.seqnames, "chr" => "")

@time let
    f = Figure(resolution = (306, 792))
    ax = Axis(f[1, 1])
    rs, chr, range1, range2 = GM.plotisoforms!(ax, "XRN2", gencode; height = 0.1,
        highlight = (["TCONS_00000039"], ["#CB3C33"]))
    GM.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
    rowsize!(f.layout, 1, rs)
    resize_to_layout!(f)
    save("figs/XRN2-long-read.png", f, px_per_unit = 4)
end