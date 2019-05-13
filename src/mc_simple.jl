"""
Naive Monte Carlo Sampling for computing default probability: P(L ⩾ l).
Returns the Monte Carlo estimate of default probability

    `sample_size` represents `(nZ, nE)`, number of samples for
        systematic risk factor `Z` and idiosyncratic risk factor `ϵ`
"""
function simple_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::Union{IO, Nothing}=nothing)
    nz, ne = sample_size
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

    Z = zeros(S)
    E = zeros(N)
    Y = zeros(N)
    ind = zeros(N)
    w1 = @view weights[:,1]
    losses = zeros(N)

    Zdist = MvNormal(S, 1)
    Edist = MvNormal(N, 1)

    estimate = 0
    estimates = zeros(Int8, ne)
    for i = 1:nz
        rand!(Zdist, Z)
        for j = 1:ne
            rand!(Edist, E)
            @. Y = $(β*Z) + denom*E
            @. ind = Y <= H[:,1]
            @. losses = w1 * ind
            estimates[j] = (sum(losses) >= l)

            k = (i-1)*ne+j
            if k % 500 == 0
                record_current(io, i, j, k, estimate, estimates)
            end
        end
        estimate = (1-1/i)*estimate + (1/i)*mean(estimates)
    end
    return estimate
end
