
import Base.MathConstants: e
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
"""
function simple_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::Union{IO, Nothing}=nothing)
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
            if io != nothing && (i*ne + j) % 500 == 0
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
function bernoulli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::Union{IO, Nothing}=nothing)
    nz, ne = sample_size
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

    Z = zeros(S)
	phi0 = zeros(N, C+1)
	phi  = @view phi0[:,2:end]
    pncz = zeros(N, C)
	pn1z = @view pncz[:,1]
    u = zeros(N)
    W = zeros(N)

    Φ = Normal()
	normcdf(x) = cdf(Φ, x)
    Zdist = MvNormal(S, 1)

    estimates = BitArray{1}()
    for i = 1:nz
        rand!(Zdist, Z)
        @. phi = normcdf(H - $(β*Z) / denom)
		pncz .= diff(phi0; dims=2)
        display(pn1z)
        for j = 1:ne
            rand!(u)
            @. W = (pn1z >= u)
            L = sum(weights[:,1] .* W)
            push!(estimates, (L >= l))
            if io != nothing && (i*ne + j) % 500 == 0
                println(io, string(mean(estimates)))
                flush(io)
            end
        end
    end
    return mean(estimates)
end


# " Computes Ψ = Σ_n ln( pn1z ⋅ e^{θ⋅w_n} ) "
# Ψ(θ, p, w) = sum(@. log(p*e^(θ*w)))
#
# " Minimizatinon objective function for inner level twisting "
# innerlevel_objective(p, θ, w, l) = Ψ(θ, p, w) - θ*l
#
# function outerlevel_twisting!(μ)
# 	fill!(μ, 0)
# end
#
# " Computes θ = argmin_θ { -θl + Ψ(θ, Z) } "
# function innerlevel_twisting!(θ)
# 	fill!(θ, 0)
# end
#
#
# " Given twisting parameter θ, compute twisted bernoulli probability q from p"
# function twisted_bernoulli!(q, θ, p)
#
# end
#
# """
#     glassermanli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::IO=Union{IO, Nothing}=nothing)
#
#     Monte Carlo Sampling for computing default probability: P(L ⩾ l), with importance sampling
#         of Z and weights w
#        ⋅ Z ∼ N(μ, I) where μ is shifted mean that minimizes variance
#        ⋅ W ∼ q       where q is shifted bernoulli probability that minimizes upper bound on the variance
# """
# function glassermanli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::IO=Union{IO, Nothing}=nothing)
#     nz, ne = sample_size
#     (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)
#
#     μ = zeros(S)
#     Z = zeros(S)
#     W = zeros(N)
#
#     Φ = Normal()
#
#     estimates = BitArray{1}()
#     for i = 1:nz
#         outerlevel_twisting!(μ)
#         rand!(MvNormal(μ, 1), Z)
#         @. pinv1 = (H[:,1] - $(β*Z)) / denom
#         for i = eachindex(pinv1)
#             pn1z[i] = cdf(Φ, pinv1[i])
#         end
#         for j = 1:ne
#             innerlevel_twisting!(θ)
#
#
#             rand!(u)
#             @. ind = (pn1z >= u)
#             L = sum(weights[:,1] .* ind)
#             push!(estimates, (L >= l))
#             if (i*ne + j) % 500 == 0
#                 println(io, string(mean(estimates)))
#                 flush(io)
#             end
#         end
#     end
#     return mean(estimates)
# end
