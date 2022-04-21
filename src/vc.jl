idx = parse(Int, ARGS[1])
include(joinpath(@__DIR__, "isoform-genetics.jl"))

function fit_uni(Y, X, V)
    model = MRVC(Y, X, V)
    MRVCs.fit!(model, maxiter = 3000, reml = true)
    h², h²se = calculate_h2(model)
    model_mle = MRVC(Y, X, V)
    MRVCs.fit!(model_mle, maxiter = 3000)
    h²_mle, h²se_mle = calculate_h2(model_mle)
    model0 = MRVC(Y, X, V[end])
    MRVCs.fit!(model0, maxiter = 3000)
    p = calculate_p(model_mle, model0)
    vcat(
        [model.Σ[k][1] for k in 1:3], 
        [sqrt(model.Σcov[k, k]) for k in 1:3],
        vec(h²), vec(h²se),
        [model_mle.Σ[k][1] for k in 1:3], 
        [sqrt(model_mle.Σcov[k, k]) for k in 1:3],
        vec(h²_mle), vec(h²se_mle), p
        )
end

function runvc_uni!(gene::Gene)
    @info "Fitting univariate model"
    V = [gene.grmcis, gene.grmtrans, Matrix{Float64}(I, size(gene.grmcis))]
    @info "Fitting gene"
    success_gene = Int[]
    try
        result = permutedims(fit_uni(gene.Yg, gene.X, V))
        push!(success_gene, 1)
    catch e
        println(e)
        result = Matrix{Float64}(undef, 0, 25)
    end
    d = size(gene.Yi, 2)
    success_isoform = Int[]
    for i in 1:d
        @info "Fitting $i / $d isoform"
        try
            result = vcat(result, permutedims(fit_uni(gene.Yi[:, i], gene.X, V)))
            push!(success_isoform, i)
        catch e
            println(e)
        end
    end
    result = DataFrame(result, [
        :var_cis_reml, :var_trans_reml, :var_res_reml, :se_var_cis_reml,
        :se_var_trans_reml, :se_var_res_reml, :h2_cis_reml, :h2_trans_reml,
        :h2_reml, :se_h2_cis_reml, :se_h2_trans_reml, :se_h2_reml, 
        :var_cis_mle, :var_trans_mle, :var_res_mle, :se_var_cis_mle, 
        :se_var_trans_mle, :se_var_res_mle, :h2_cis_mle, :h2_trans_mle,
        :h2_mle, :se_h2_cis_mle, :se_h2_trans_mle, :se_h2_mle, :p
        ])
    result.gene_name .= gene.gene_name
    result.gene_id .= gene.gene_id
    if isempty(success_gene)
        result.id = gene.expressed_isoforms[success_isoform]
        result.feature = fill("isoform", length(success_isoform))
        ind = findall(result.p .< 0.05)
        !isempty(ind) ? gene.heritable_isoforms = success_isoform[ind] : nothing
    else
        result.id = vcat(gene.gene_id, gene.expressed_isoforms[success_isoform])
        result.feature = vcat("gene", fill("isoform", length(success_isoform)))
        result.p[1] < 0.05 ? gene.gene_heritable = true : nothing
        ind = findall(result.p[2:end] .< 0.05)
        !isempty(ind) ? gene.heritable_isoforms = success_isoform[ind] : nothing
    end
    result.ncisnps .= gene.ncisnps
    select!(result, :gene_name, :gene_id, :id, :feature, :)
    CSV.write(joinpath(@__DIR__, "../results/univariate/$(gene.gene_id).tsv"), result, delim = "\t")
    if !isnothing(gene.heritable_isoforms) && length(gene.heritable_isoforms) > 1
        remove = Int[]
        d = length(gene.heritable_isoforms)
        for i in 1:d
            if i in remove
                continue
            end    
            for j in (i + 1):d
                if cor(gene.Yi[:, gene.heritable_isoforms[i]], gene.Yi[:, gene.heritable_isoforms[j]]) == 1
                    push!(remove, j)
                    continue
                end        
            end
        end
        gene.heritable_isoforms = gene.heritable_isoforms[Not(remove)]
    end
end

function runvc_bi(gene::Gene)
    if isnothing(gene.heritable_isoforms) || length(gene.heritable_isoforms) == 1
        return
    end
    result = DataFrame(
        gene_name           = String[],
        gene_id             = String[],
        id                  = String[],
        pair                = String[],
        var_cis_reml        = Float64[],
        covar_cis_reml      = Float64[],
        var_trans_reml      = Float64[],
        covar_trans_reml    = Float64[],
        var_res_reml        = Float64[],
        covar_res_reml      = Float64[],
        se_var_cis_reml     = Float64[],
        se_covar_cis_reml   = Float64[],
        se_var_trans_reml   = Float64[],
        se_covar_trans_reml = Float64[],
        se_var_res_reml     = Float64[],
        se_covar_res_reml   = Float64[],
        h2_cis_reml         = Float64[],
        h2_trans_reml       = Float64[],
        h2_reml             = Float64[],
        se_h2_cis_reml      = Float64[],
        se_h2_trans_reml    = Float64[],
        se_h2_reml          = Float64[],
        rg_cis_reml         = Float64[],
        rg_trans_reml       = Float64[],
        re_reml             = Float64[],
        se_rg_cis_reml      = Float64[],
        se_rg_trans_reml    = Float64[],
        se_re_reml          = Float64[],
        p_cis               = Float64[],
        p_trans             = Float64[],
        p_res               = Float64[],
        ncisnps             = Int[]
        )
    V = [gene.grmcis, gene.grmtrans, Matrix{Float64}(I, size(gene.grmcis))]
    d = length(gene.heritable_isoforms)
    @info "Fitting bivariate model for $d isoforms"
    for i in 1:d
        for j in (i + 1):d
            @info "Fitting pair ($i, $j)"
            model = MRVC(gene.Yi[:, [gene.heritable_isoforms[i], gene.heritable_isoforms[j]]], gene.X, V)
            MRVCs.fit!(model, maxiter = 3000, reml = true) 
            h², h²se = calculate_h2(model)
            r₉, r₉se = calculate_rg(model)
            r₉_parsed = [r₉[k][1, 2] for k in 1:3]
            r₉se_parsed = [r₉se[k][1, 2] for k in 1:3]
            zs = (r₉_parsed ./ r₉se_parsed).^2
            ps = ccdf.(Chisq(1), zs)
            push!(result, vcat(
                gene.gene_id,
                gene.gene_name,
                gene.expressed_isoforms[gene.heritable_isoforms[i]],
                gene.expressed_isoforms[gene.heritable_isoforms[j]],
                reduce(vcat, [model.Σ[k][1:2] for k in 1:3]),
                sqrt.(diag(model.Σcov)[[1, 2, 4, 5, 7, 8]]),
                h²[:, 1], h²se[:, 1],
                r₉_parsed,
                r₉se_parsed,
                ps, 
                gene.ncisnps
                ))
            push!(result, vcat(
                gene.gene_id,
                gene.gene_name,
                gene.expressed_isoforms[gene.heritable_isoforms[j]],
                gene.expressed_isoforms[gene.heritable_isoforms[i]],
                reduce(vcat, [model.Σ[k][4:-1:3] for k in 1:3]),
                sqrt.(diag(model.Σcov)[[3, 2, 6, 5, 9, 8]]),
                h²[:, 2], h²se[:, 2],
                r₉_parsed,
                r₉se_parsed,
                ps, 
                gene.ncisnps
                ))
        end
    end
    CSV.write(joinpath(@__DIR__, "../results/bivariate/$(gene.gene_id).tsv"), result, delim = "\t")
end

function runvc_mul(gene::Gene)
    if isnothing(gene.heritable_isoforms) || length(gene.heritable_isoforms) < 3
        return
    end
    d = length(gene.heritable_isoforms)
    np = (d * (d + 1)) >> 1
    ind = ones(Int, d)
    for i in 2:length(ind)
        ind[i] = ind[i - 1] + (d - i + 2) 
    end
    V = [gene.grmcis, gene.grmtrans, Matrix{Float64}(I, size(gene.grmcis))]
    d = length(gene.heritable_isoforms)
    @info "Fitting multivariate model for $d isoforms"
    model = MRVC(gene.Yi[:, gene.heritable_isoforms], gene.X, V)
    MRVCs.fit!(model, maxiter = 3000, reml = true, verbose = true)
    h², h²se = calculate_h2(model)
    r₉, r₉se = calculate_rg(model)
    h2rg = reduce(hcat, [MRVCs.vech(r₉[i]) for i in 1:3])
    h2rgse = reduce(hcat, [MRVCs.vech(r₉se[i]) for i in 1:3])
    for i in 1:3
        h2rg[ind, i] = h²[i, :]
        h2rgse[ind, i] = h²se[i, :]
    end
    zs = (h2rg ./ h2rgse).^2
    ps = ccdf.(Chisq(1), zs)
    ps[ind, :] .= NaN
    result = hcat(
        reduce(hcat, [MRVCs.vech(model.Σ[i]) for i in 1:3]),
        reduce(hcat, [diag(model.Σcov)[((i - 1) * np + 1):(i * np)] for i in 1:3]),
        h2rg, h2rgse, ps
        )
    result = DataFrame(result, [
        "var/covar_cis_reml", "var/covar_trans_reml", "var/covar_res_reml", 
        "se_var/covar_cis_reml", "se_var/covar_trans_reml", "se_var/covar_res_reml",
        "h2/rg_cis_reml", "h2/rg_trans_reml", "h2/re_reml",
        "se_h2/rg_cis_reml", "se_h2/rg_trans_reml", "se_h2/re_reml", 
        "p_cis", "p_trans", "p_res"
        ])
    ids = gene.expressed_isoforms[gene.heritable_isoforms]
    txind = MRVCs.vech([[j, i] for i in 1:d, j in 1:d])
    ids = reduce(vcat, [permutedims(ids[txind[i]]) for i in 1:length(txind)])    
    result.gene_name .= gene.gene_name
    result.gene_id .= gene.gene_id
    result.id = ids[:, 1]
    result.pair = ids[:, 2]
    result.ncisnps .= gene.ncisnps
    select!(result, :gene_name, :gene_id, :id, :pair, :)
    CSV.write(joinpath(@__DIR__, "../results/multivariate/$(gene.gene_id).tsv"), result, delim = "\t")
end

for folder in ["univariate", "bivariate", "multivariate"]
    ispath(joinpath(@__DIR__, "../results/$folder")) || 
        mkpath(joinpath(@__DIR__, "../results/$folder"))
end

function main(gencode, expr, expri, cov, geno, idx)
    for gene_id in union(expr.pid, expri.pid)[idx:(idx + 99)]
        gene_name = gencode.gene_name[findfirst(isequal(gene_id), gencode.gene_id)]
        (isfile(joinpath(@__DIR__, "../results/univariate/$(gene_id).tsv")) &&
	        !iseile(joinpath(@__DIR__, "../data/$(gene_name)-trans.bed"))) && continue
        @info "Working on $(gene_name) / $(gene_id)"
        @time gene = Gene(gene_id, gencode, expr, expri, cov, 1e6, geno, "both")
        @time runvc_uni!(gene)
        @time runvc_bi(gene)
        @time runvc_mul(gene)
        for plink in ["bed", "bim", "fam"]
            rm("data/$(gene.gene_name)-cis.$(plink)")
            rm("data/$(gene.gene_name)-trans.$(plink)")
        end
    end
end

main(gencode, expr, expri, cov, geno, idx)