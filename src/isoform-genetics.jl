using CSV, DataFrames, CairoMakie, ColorSchemes
using SnpArrays, GLM, Arrow, Distributions
using GeneticsMakie, MultiResponseVarianceComponentModels
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

mutable struct Gene
    gene_name          :: String
    gene_id            :: String
    chr                :: String
    start              :: Real
    stop               :: Real
    all_isoforms       :: Vector{String}
    expressed_isoforms :: Vector{String}
    heritable_isoforms :: Union{Nothing, Vector{Int}}
    gene_heritable     :: Bool
    Yg                 :: Vector{Float64}
    Yi                 :: Union{Vector{Float64}, Matrix{Float64}}
    X                  :: Matrix{Float64} # covariates
    snparray           :: Union{Nothing, Matrix{Float64}} # cis SNPs only
    snpinfo            :: Union{Nothing, DataFrame}
    grmcis             :: Union{Nothing, Matrix{Float64}}
    grmtrans           :: Union{Nothing, Matrix{Float64}}
    ncisnps            :: Union{Nothing, Int}
end

function Gene(
    gene     :: String,
    gencode  :: DataFrame,
    expr     :: DataFrame,
    expri    :: DataFrame,
    X        :: Union{Nothing, Matrix{Float64}},
    window   :: Real,
    geno     :: Union{Nothing, SnpData},
    whichgrm :: Union{Nothing, String},
    )

    if startswith(gene, "ENSG")
        gene_id = gene
        gene_name = gencode.gene_name[findfirst(isequal(gene), gencode.gene_id)]
    else
        gene_id = gencode.gene_id[findfirst(isequal(gene), gencode.gene_name)]
        gene_name = gene
    end
    chr, start, stop = GM.findgene(gene_name, gencode)
    df = filter(x -> (x.gene_id == gene_id) && (x.feature .== "transcript"), gencode)
    all_isoforms = df.transcript_id
    df = filter(x -> x.pid == gene_id, expr)
    Yg = permutedims(Matrix{Float64}(df[:, 7:end]))
    Yg = vec(Yg)
    df = filter(x -> x.pid == gene_id, expri)
    Yi = permutedims(Matrix{Float64}(df[:, 7:end]))
    size(Yi, 2) == 1 ? Yi = vec(Yi) : nothing
    expressed_isoforms = df.gid
    heritable_isoforms = nothing
    gene_heritable = false
    # isnothing(X) ? X = ones(size(Yg, 1), 1) : nothing
    if !isnothing(geno) && isfile(joinpath(@__DIR__, "../data/$(gene_name)-cis.bed"))
        genocis = SnpData(joinpath(@__DIR__, "../data/$(gene_name)-cis"))
        snparray, snpinfo = (convert(Matrix{Float64}, genocis.snparray), genocis.snp_info)
        ncisnps = size(genocis.snparray, 2)
    elseif !isnothing(geno)
        SnpArrays.filter(geno; des = joinpath(@__DIR__, "../data/$(gene_name)-cis"), 
            f_snp = x -> (x[:chromosome] == chr && x[:position] > (start - window) && x[:position] < (stop + window)))
        genocis = SnpData(joinpath(@__DIR__, "../data/$(gene_name)-cis"))
        snparray, snpinfo = (convert(Matrix{Float64}, genocis.snparray), genocis.snp_info)
        ncisnps = size(genocis.snparray, 2)
    else
        snparray, snpinfo, ncisnps = nothing, nothing, nothing
    end
    if whichgrm == "cis"
        grmcis = grm(genocis.snparray, method = :GRM)
        grmcis .*= 2
        grmtrans = nothing
    elseif whichgrm == "both"
        grmcis = grm(genocis.snparray, method = :GRM)
        grmcis .*= 2
        if !isfile(joinpath(@__DIR__, "../data/$(gene_name)-trans.bed"))
            SnpArrays.filter(geno; des = joinpath(@__DIR__, "../data/$(gene_name)-trans"), 
                f_snp = x -> !(x[:chromosome] == chr && x[:position] > (start - window) && x[:position] < (stop + window)))
        end
        genotrans = SnpData(joinpath(@__DIR__, "../data/$(gene_name)-trans"))
        grmtrans = grm(genotrans.snparray, method = :GRM)
        grmtrans .*= 2
    else
        grmcis = grmtrans = nothing
    end
    Gene(
        gene_name, gene_id, chr, start, stop, 
        all_isoforms, expressed_isoforms, heritable_isoforms,
        gene_heritable, Yg, Yi, X, snparray, snpinfo, grmcis, grmtrans, ncisnps
        )
end

@info "Loading gene expression data"
if !isfile(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-gene.BED.arrow"))
    @time expr = CSV.read(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-gene.BED.gz"), delim = "\t", DataFrame)
    select!(expr, Not([:CMC_MSSM_307, :CMC_PITT_077, :Br1977, :Br1876, :Br2173]))
    @time Arrow.write(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-gene.BED.arrow"), expr)
end
@time expr = Arrow.Table(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-gene.BED.arrow")) |> DataFrame
@assert (24_905, 861) == size(expr)

@info "Loading isoform expression data"
if !isfile(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-isoform.BED.arrow"))
    @time expri = CSV.read(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-isoform.BED.gz"), delim = "\t", DataFrame)
    select!(expri, Not([:CMC_MSSM_307, :CMC_PITT_077, :Br1977, :Br1876, :Br2173]))
    @time Arrow.write(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-isoform.BED.arrow"), expri)
end
@time expri = Arrow.Table(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-isoform.BED.arrow")) |> DataFrame
@assert (93_293, 861) == size(expri)

@info "Loading covariates for expression"
@time covariates = CSV.read(joinpath(@__DIR__, "../data/expression/PsychENCODE-EUR-covariates.tsv"), delim = "\t", DataFrame)
covariates = covariates[Not([257, 439, 665, 757, 831]), :]
covariates = Matrix(covariates)
@assert (855, 45) == size(covariates)

@info "Loading GENCODE annotation"
if !isfile(joinpath(@__DIR__, "../data/gencode.v19.annotation.parsed.gtf.arrow"))
    run(`curl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz --output ./data/gencode.v19.annotation.gtf.gz`)
    @time gencode = CSV.read(joinpath(@__DIR__, "../data/gencode.v19.annotation.gtf.gz"), DataFrame,
        delim = "\t", skipto = 6, header = ["seqnames", "source", "feature", "start", 
        "end", "score", "strand", "phase", "info"])
    GeneticsMakie.parsegtf!(gencode)
    select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_id, :gene_name, :gene_type, :transcript_id)
    @time Arrow.write(joinpath(@__DIR__, "../data/gencode.v19.annotation.parsed.gtf.arrow"), gencode)
    run(`rm ./data/gencode.v19.annotation.gtf.gz`)
end 
@time gencode = Arrow.Table(joinpath(@__DIR__, "../data/gencode.v19.annotation.parsed.gtf.arrow"))|> DataFrame
@assert (2_619_444, 9) == size(gencode)

@info "Loading genotype data"
@time geno = SnpData(joinpath(@__DIR__, "../data/genotype/Capstone4.HRC.European.unique.frontal.nochr.filter.unrelated"))
@assert (855, 4_685_674) == size(geno.snparray)
@assert names(expri)[7:end] == geno.person_info.iid

function runqtl(gene::Gene)
    @info "Running cis-eQTL analysis"
    d = size(gene.Yi, 2) + 1
    qtls = Vector{DataFrame}(undef, d)
    data = DataFrame(g = gene.snparray[:, 1], e = gene.Yg)
    data = hcat(data, DataFrame(gene.X[:, 2:end], :auto))
    for i in 1:d
        @info "Fitting $i / $d phenotype"
        df = DataFrame(
            SNP = String[], BETA = Float64[], SE = Float64[], Z = Float64[], 
            P = Float64[], BP = Int64[], A1 = String[], A2 = String[], ID = String[]
            )
        for j in 1:size(gene.snparray, 2)
            any(isnan.(gene.snparray[:, j])) ? continue : nothing
            if i > 1
                data.e = gene.Yi[:, i - 1]
            end
            data.g =  gene.snparray[:, j]
            model = lm(term(:e) ~ sum(term.(Symbol.(names(data, Not(:e))))), data)
            i == 1 ? id = gene.gene_id : id = gene.expressed_isoforms[i - 1]
            push!(df,
                [gene.snpinfo.snpid[j], coef(model)[2], coeftable(model).cols[2][2], 
                coeftable(model).cols[3][2], coeftable(model).cols[4][2],
                gene.snpinfo.position[j], gene.snpinfo.allele1[j],
                gene.snpinfo.allele2[j], id]
                )
        end
        storage = Vector{String}(undef, length(df.SNP))
        for i in 1:nrow(df)
            snp = split(df.SNP[i], ";")
            length(snp) > 1 ? storage[i] = snp[2] : storage[i] = snp[1]
        end
        df.SNP = storage
        df.CHR .= gene.chr
        qtls[i] = df
    end
    qtls, reduce(vcat, qtls)
end

function calculate_rp(gene::Gene, univariate)
    storage = filter(row -> row.feature == "isoform" && row.gene_id == gene.gene_id && row.p < 0.05, univariate)
    ind = [findfirst(isequal(iso), gene.expressed_isoforms) for iso in storage.id]
    d = length(ind)
    regressed = Matrix{Float64}(undef, size(gene.Yi, 1), d)
    data = DataFrame(e = gene.Yg)
    data = hcat(data, DataFrame(gene.X[:, 2:end], :auto))
    for i in 1:d
        data.e = gene.Yi[:, ind[i]]
        model = lm(term(:e) ~ sum(term.(Symbol.(names(data, Not(:e))))), data)
        regressed[:, i] = data.e - Matrix(data[:, 2:end]) * coef(model)[2:end]
    end
    cor(regressed)
end

function loadgwas()
    @info "Load GWAS results."
    delete!(GM.gwas, "mdd")
    gwas = []
    for key in keys(GM.gwas)
        push!(gwas, DataFrame(Arrow.Table(joinpath.("data/gwas/processed/", key * ".tsv.arrow"))))
    end
    @info "Loading 1000 Genomes"
    @time kgp_raw = SnpData(joinpath(@__DIR__, "../data/1kg/kgp.EUR.maf0.05.geno"))
    gwas, kgp_raw
end

function subsetref(ref::SnpData, chr::AbstractString, range1::Real, range2::Real, path::AbstractString)
    SnpArrays.filter(ref, trues(size(ref)[1]), GeneticsMakie.findlocus(ref, chr, range1, range2); des = path)
    SnpData(path)
end

function subsetgwas(gwas, chr, range1, range2)
    dfs = Vector{DataFrame}(undef, length(gwas))
    for i in 1:length(gwas)
        dfs[i] = gwas[i][GeneticsMakie.findlocus(gwas[i], chr, range1, range2), :]
    end
    dfs
end

issig(P::AbstractVector; p = 5e-8) = any(P .< p)
issig(df::DataFrame; p = 5e-8) = issig(df.P; p = p)