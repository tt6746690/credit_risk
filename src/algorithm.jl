
import Random: rand!
import Distributions: MvNormal
import Statistics: mean

include("utils.jl")
include("parameter.jl")

"""
    simple_mc(parameter::Parameter, sample_size::Tuple{Number, Number})

Naive Monte Carlo Sampling for computing default probability: P(L ⩾ l).
Returns the Monte Carlo estimate of default probability

    `sample_size` represents `(nZ, nE)`, number of samples for
        systematic risk factor `Z` and idiosyncratic risk factor `ϵ`
"""
function simple_mc(parameter::Parameter, sample_size::Tuple{Number, Number})
    nz, ne = sample_size
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

    Z = zeros(S)
    E = zeros(N)
    estimates = []

    for i = 1:nz
        # ∼𝓝(0, I_S)
        rand!(MvNormal(S, 1), Z)
        for j = 1:ne
            # ∼𝓝(0, I_N)
            rand!(MvNormal(N, 1), E)
            Y = β*Z + denom.*E
            ind = Y .<= H[:,1]
            L = sum(weights[:,1] .* ind)
            push!(estimates, (L >= l))
        end
    end
    return mean(estimates)
end
