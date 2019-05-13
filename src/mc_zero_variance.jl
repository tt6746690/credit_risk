

" Computes μ for the approximating standard normal of loss probability"
function approx_μ(pnc, weights, A)
    @. A = pnc * weights
    sum(A)
end


" Computes σ for the approximating standard normal of loss probability"
function approx_σ(pnc, lgc, ead, A, B, C)
    n, c = size(lgc)
    i = 1

    lgcs = Vector{SubArray}(undef, c)
    pncs = Vector{SubArray}(undef, c)
    for idx = 1:c
        lgcs[idx] = @view lgc[:,idx]
        pncs[idx] = @view pnc[:,idx]
    end

    for a = 1:c
        for b = 1:(a-1)
            @. A[:,i] = (lgcs[a] - lgcs[b])^2 * pncs[a] * pncs[b]
            i += 1
        end
    end

    sum!(B, A)
    @. C = ead^2 * B
    sqrt(sum(C) / sum(ead)^2)
end


"""
One level Monte Carlo simulation using
    conditional portfolio loss approaches in distribution to standard normal distribution
    (L - μ) / σ  ⟶ N(0, 1)
"""
function onelvl_mc(parameter::Parameter, nz::Int64, io::Union{IO, Nothing}=nothing)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

    C₀, C₁ = 0.7164, 31.395
    w = zeros(N)
    w = ead ./ sum(ead)
    M = zeros(N)

    Z = zeros(S)
    βZ = zeros(N)
    Zdist = MvNormal(S, 1)
    phi0 = zeros(N, C+1)
    phi  = @view phi0[:,2:end]
    pnc = zeros(N, C)
    approx_μ_A = similar(weights)
    approx_σ_A = zeros(N, convert(Int, (C-1)*C/2))
    approx_σ_B = zeros(N)
    approx_σ_C = zeros(N)

    estimates = zeros(nz)
    for i = 1:nz
        rand!(Zdist, Z)
        mul!(βZ, β, Z)
        @. phi = normcdf((H - βZ) / denom)
        diff!(pnc, phi0; dims=2)

        μ = approx_μ(pnc, weights, approx_μ_A)
        σ = approx_σ(pnc, lgc, ead, approx_σ_A, approx_σ_B, approx_σ_C)
        p = 1 - normcdf((l-μ)/σ)

        estimates[i] = p

        (i) % 1 == 0 && begin
            if io != nothing
                println(io, string(mean(estimates[1:i])))
                flush(io)
            end
        end
    end
    return mean(estimates)
end



"""
Approximating the zero variance importance sampler
"""
function onelvlISCLT_mc(parameter::Parameter, nz::Int64, io::Union{IO, Nothing}=nothing)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

    Z = zeros(S)
    βZ = zeros(N)
    phi0 = zeros(N, C+1)
    phi  = @view phi0[:,2:end]
    pnc = zeros(N, C)
    approx_μ_A = similar(weights)
    approx_σ_A = zeros(N, convert(Int, (C-1)*C/2))
    approx_σ_B = zeros(N)
    approx_σ_C = zeros(N)
    Zdist = MvNormal(S, 1)

    # sample from π(z) with slice sampler

    function π(Z)
        mul!(βZ, β, Z)
        @. phi = normcdf((H - βZ) / denom)
        diff!(pnc, phi0; dims=2)
        μ = approx_μ(pnc, weights, approx_μ_A)
        σ = approx_σ(pnc, lgc, ead, approx_σ_A, approx_σ_B, approx_σ_C)
        p = (1 - normcdf((l-μ)/σ)) * pdf(Zdist, Z)
        return max(p, 1e-100)
    end

    burn_ratio = 0.1
    initial_point = zeros(S)
    width = (S == 1) ? 15. : 0.5
    stepin = (S == 1) ? false : true

    Zs = slicesample(initial_point, π, nz;
        step_limit=20,
        width=width,
        burn=Int(burn_ratio * nz),
        thin=3,
        stepin=stepin)

    # TODO: train zero variance IS sampler with Gaussian Mixture

    estimates = zeros(nz)
    for i = 1:nz
        Z = Zs[:,i]
        mul!(βZ, β, Z)
        @. phi = normcdf((H - βZ) / denom)
        diff!(pnc, phi0; dims=2)

        μ = approx_μ(pnc, weights, approx_μ_A)
        σ = approx_σ(pnc, lgc, ead, approx_σ_A, approx_σ_B, approx_σ_C)
        p = (1 - normcdf((l-μ)/σ)) * pdf(Zdist, Z) / π(Z)
        display("p:$p I:$((1 - normcdf((l-μ)/σ))) * pdf: $(pdf(Zdist, Z)) / pipdf: $(π(Z))")

        estimates[i] = p

        (i) % 1 == 0 && begin
            if io != nothing
                println(io, string(mean(estimates[1:i])))
                flush(io)
            end
        end
    end
    return mean(estimates)
end
