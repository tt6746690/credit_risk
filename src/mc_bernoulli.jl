"""
    bernoulli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::Union{IO, Nothing}=nothing)

Monte Carlo Sampling for computing default probability: P(L ⩾ l), where
inner level sampling of `ε` is replaced with sampling of bernoulli random variables

W[n] := weights[n],    p[n]
     := 0         ,    1 - p[n]

where p[n] is a function of outer level sample `Z`

    `sample_size` represents `(nZ, nE)`, number of samples for
        systematic risk factor `Z` and idiosyncratic risk factor `ϵ`

For n=2500, (nz,ne) = (1000,1000), ~10s

"""
function bernoulli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::Union{IO, Nothing}=nothing)
    nz, ne = sample_size
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

    Z = zeros(S)
    phi0 = zeros(N, C+1)
    phi  = @view phi0[:,2:end]
    pnc = zeros(N, C)
    pn1 = @view pnc[:,1]
    u = zeros(N)
    W = zeros(N)
    w1 = @view weights[:,1]
    losses = zeros(N)

    Zdist = MvNormal(S, 1)

    estimate = 0
    estimates = zeros(ne)
    for i = 1:nz
        rand!(Zdist, Z)
        @. phi = normcdf((H - $(β*Z)) / denom)
        diff!(pnc, phi0; dims=2)
        for j = 1:ne
            rand!(u)
            @. W = (pn1 >= u)
            @. losses = w1 * W
            estimates[j] = (sum(losses) >= l)

            k = (i-1)*ne+j
            if k % 10 == 0
                record_current(io, i, j, k, estimate, estimates)
            end
        end
        estimate = (1-1/i)*estimate + (1/i)*mean(estimates)
    end
    return estimate
end
