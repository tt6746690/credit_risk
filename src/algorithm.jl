import Distributions: MvNormal

include("parameter.jl")

"""
    simple_mc(parameter::Parameter, sample_size::Tuple{Number, Number})

Naive Monte Carlo Sampling for computing default probability: P(L ⩾ l).
Returns a generator outputing the Monte Carlo estimate of default probability
for each sampleing iteration.

    `sample_size` represents `(nZ, nE)`, number of samples for
        systematic risk factor `Z` and idiosyncratic risk factor `ϵ`
"""
function simple_mc(parameter::Parameter, sample_size::Tuple{Number, Number})

    (N, C, S, l, cmm, ead, lgc, cn, β, H) = unpack(parameter)

    # Outer-level sampling ∼𝓝(0, I_S)
    sampleZ = rand(MvNormal(S, 1))

    # Computing 𝟙_{nth creditor in credit state c}
    denom = @. sqrt(1 - β'β)

end
