include("utils.jl")

" Parameters related to the credit risk problem "
struct Parameter
    # Number of creditors
    N::Int64
    # Number of credit states
    C::Int64
    # Dimension for Systematic Risk Factor Z, i.e. Z ∼  N(0, I_S)
    S::Int64
    # Tail for computing P(L ⩾ l)
    l::Float64


    # Credit Migration Matrix,
    # CMM[n, c] has probability of n-th creditor moving to credit state c
    cmm::Array{Float64, 2}
    # Exposure at Default
    # Value lost if n-th creditor is in credit state c=1, i.e. defaults
    ead::Array{Float64, 1}
    # Percentage loss/gain
    # lgc[n, c] represent loss from creditor n move to credit state c
    lgc::Array{Float64, 2}

    # Initial credit state for each creditor
    # In case of binary credit states, initially at c=2, i.e. non-default state
    cn::Array{Int64, 1}

    # β[n, c] indicates n-th creditor's sensitivity to systematic risk factor Z
    β::Array{Float64, 2}

    # Threshold for creditor n migrating from current credit state `cn` to any `c`
    # H[n, c] represent threshold for migrating to c from current state
    H::Array{Float64, 2}


    # Computed values from above fields
    # √(1 - βᵗβ)
    denom::Array{Float64, 1}
    # weights
    # weights[n, c] := w_n * lgc[n, c] where w_n := ead[n] ./ sum(ead)
    weights::Array{Float64, 2}


    function Parameter(N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights)
        @checksize (N, C)   cmm
        @checksize (N,)     ead
        @checksize (N, C)   lgc
        @checksize (N, S)   β
        @checksize (N,)     cn
        @checksize (N, C)   H
        @checksize (N,)     denom
        @checksize (N, C)   weights
        new(N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights)
    end
end


function Parameter(N, C, S, l)
    N >= 0 || throw(ArgumentError("Invalid number of creditors: $N"))
    C >= 0 || throw(ArgumentError("Invalid number of credit states: $C"))
    S >= 0 || throw(ArgumentError("Invalid Dimension for systematic risk factor: $S"))
    # 0 <= l <= 1 || throw(ArgumentError("Invalid tail probability $l"))

    p = @. 0.01*(1 + sin(16π*(1:N)/N))
    cmm = zeros(N, C)
    cmm[:,1], cmm[:,2] = p, 1 .- p  # binary case ...
    ead = 0.5 .+ rand(N)
    lgc = zeros(N, C)
    lgc[:,1] = @. floor(5(1:N)/N)^2
    cn = fill(2, N)

    β = rand(N, S)
    # 1/sqrt(S) * Unif(-1,1)
    β = @. (2β-1)/sqrt(S)
    # β = fill(1/sqrt(S), N, S)   # No idiosyncratic risk factor
    # β = zeros(N, S)             # No systematic risk factor

    # H[n, c] = inverse_unit_Gaussian(∑ᵧ cmm[c(n), γ])
    cum_cmm = cumsum(cmm, dims=2)
    H = invnormcdf.(cum_cmm)

    # Computed
    denom = @. sqrt(1 - $sum((β).^2, dims=2))
    denom = reshape(denom, :)
    weights = ead ./ sum(ead)
    weights = weights .* lgc
    # Idea is the larger weights/exposure, there is a larger increase in
    # default probability during the inner level twisting pnc → qnc

    return Parameter(N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights)
end


function unpack(p::Parameter)
    return (
        p.N,
        p.C,
        p.S,
        p.l,
        p.cmm,
        p.ead,
        p.lgc,
        p.cn,
        p.β,
        p.H,
        p.denom,
        p.weights
    )
end
