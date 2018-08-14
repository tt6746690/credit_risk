
import Random: rand!
import Distributions: cdf, Normal, MvNormal
import Statistics: mean

include("utils.jl")
include("parameter.jl")

"""
    simple_mc(parameter::Parameter, sample_size::Tuple{Number, Number})

Naive Monte Carlo Sampling for computing default probability: P(L ⩾ l).
Returns the Monte Carlo estimate of default probability

    `sample_size` represents `(nZ, nE)`, number of samples for
        systematic risk factor `Z` and idiosyncratic risk factor `ϵ`

~20s, 18GB for (nz, ne) == (500, 500)
"""
function simple_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::IO)
    nz, ne = sample_size
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

    Z = zeros(S)
    E = zeros(N)
    Y = zeros(N)
    ind = zeros(N)

    Zdist = MvNormal(S, 1)
    Edist = MvNormal(N, 1)

    estimates = BitArray{1}()
    for i = 1:nz
        rand!(Zdist, Z)
        for j = 1:ne
            rand!(Edist, E)
            @. Y = $(β*Z) + denom*E
            @. ind = Y <= H[:,1]
            L = sum(weights[:,1] .* ind)
            push!(estimates, (L >= l))
            if (i*ne + j) % 10 == 0
                println(io, string(mean(estimates)))
                flush(io)
            end
        end
    end
    return mean(estimates)
end


"""
    bernoulli_mc(parameter::Parameter, sample_size::Tuple{Number, Number})

Monte Carlo Sampling for computing default probability: P(L ⩾ l), where
inner level sampling of `ε` is replaced with sampling of bernoulli random variables

W[n] := weights[n],    p[n]
     := 0         ,    1 - p[n]

where p[n] is a function of outer level sample `Z`

    `sample_size` represents `(nZ, nE)`, number of samples for
        systematic risk factor `Z` and idiosyncratic risk factor `ϵ`
"""
function bernoulli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::IO)
    nz, ne = sample_size
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

    Z = zeros(S)
    E = zeros(N)
    pinv1 = zeros(N)
    pn1z = zeros(N)
    u = zeros(N)
    ind = zeros(N)

    Φ = Normal()
    Zdist = MvNormal(S, 1)

    estimates = BitArray{1}()
    for i = 1:nz
        rand!(Zdist, Z)
        @. pinv1 = (H[:,1] - $(β*Z)) / denom
        for i = eachindex(pinv1)
            pn1z[i] = cdf(Φ, pinv1[i])
        end
        for j = 1:ne
            rand!(u)
            @. ind = (pn1z >= u)
            L = sum(weights[:,1] .* ind)
            push!(estimates, (L >= l))
            if (i*ne + j) % 500 == 0
                println(io, string(mean(estimates)))
                flush(io)
            end
        end
    end
    return mean(estimates)
end
