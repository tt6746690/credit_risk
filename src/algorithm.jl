import Optim
import Optim: optimize, ConjugateGradient

import LinearAlgebra: mul!, ⋅, dot
import Base.MathConstants: e
import Random: rand!
import Distributions: cdf, Normal, MvNormal
import Statistics: mean

include("utils.jl")
include("parameter.jl")


function record_current(io, i, j, estimate, estimates)
    if io != nothing
        prop = j/((i-1)*ne+j)
        est = (1-prop)*estimates + (prop)*mean(estimates[1:j])
        println(io, string(mean(est)))
        flush(io)
    end
end

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
            (i*ne + j) % 500 == 0 &&
                record_current(io, i, j, estimate, estimates)
        end
        estimate = (1-1/i)*estimate + (1/i)*mean(estimates)
    end
    return estimate
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
            (i*ne + j) % 500 == 0 &&
                record_current(io, i, j, estimate, estimates)
        end
        estimate = (1-1/i)*estimate + (1/i)*mean(estimates)
    end
    return estimate
end

function outerlevel_twisting!(μ)
    fill!(μ, 0)
end


function init_Ψ()
    x = nothing
    y = nothing
    " Computes Ψ = Σ_n ln( y ) where x = Σ_c p⋅e^{θ⋅w_n}"
    function Ψ(θ, p, w)
        if x == nothing && y == nothing
            x = similar(p)
            y = similar(p, axes(p, 1))
        end
        @. x = p*ℯ^(θ*w)
        sum!(y, x)
        mapreduce(log, +, y; init=0)
    end
    return Ψ
end

struct InnerLevelTwisting
    N::Int64
    C::Int64
    wp::Array{Float64, 2}
    # optimization related
    # Result of twist, wrapped into an array
    θ::Array{Float64, 1}
    Ψ::Any
    initialguess::Array{Float64, 1}

    function InnerLevelTwisting(N, C, wp, θ, Ψ, initialguess)
        @checksize (N, C)   wp
        @checksize (1,)     θ
        @checksize (1,)     initialguess
        new(N, C, wp, θ, Ψ, initialguess)
    end
end

function InnerLevelTwisting(N, C)
    wp = zeros(N, C)
    θ = [0]
    Ψ = init_Ψ()
    initialguess = [0]
    InnerLevelTwisting(N, C, wp, θ, Ψ, initialguess)
end

get_result(t::InnerLevelTwisting) = t.θ[1]
set_result!(t::InnerLevelTwisting, θ) = (t.θ[1] = θ)

" Computes inner level twisting parameter θ = argmin_θ { -θl + Ψ(θ, Z) } "
function twist!(t::InnerLevelTwisting, p, w, l)
    @. t.wp = w * p
    if l > sum(t.wp)
        function objective(θ)
            θ = θ[1]
            t.Ψ(θ, p, w) - θ*l
        end
        results = optimize(objective, t.initialguess, ConjugateGradient())
        minimizer = Optim.converged(results) ? Optim.minimizer(results) : 0
        t.initialguess[:] = minimizer
        set_result!(t, minimizer[1])
    else
        set_result!(t, 0)
    end
end


struct OuterLevelTwisting end

function twist!(t::OuterLevelTwisting, parameter::Parameter, optimizer)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)
    phi0 = zeros(N, C+1)
    phi  = @view phi0[:,2:end]
    pnc = zeros(N, C)

    Ψ = init_Ψ()
    innerlevel = InnerLevelTwisting(N, C)

    function objective(z)
        @checksize (S,) z
        @. phi = normcdf((H - $(β*z)) / denom)
        diff!(pnc, phi0; dims=2)
        twist!(innerlevel, pnc, weights, l)
        θ = get_result(innerlevel)
        θ*l - Ψ(θ, pnc, weights) + 0.5z'z
    end

    results = optimize(objective, zeros(S), optimizer)
    if Optim.converged(results)
        display("Converged!")
        return Optim.minimizer(results)
    else
        return zeros(S)
    end
end


"""
    glassermanli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::IO=Union{IO, Nothing}=nothing)

    Monte Carlo Sampling for computing default probability: P(L ⩾ l), with importance sampling
        of Z and weights w
       ⋅ Z ∼ N(μ, I) where μ is shifted mean that minimizes variance
       ⋅ W ∼ q       where q is shifted bernoulli probability that minimizes upper bound on the variance
"""
function glassermanli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, extra_params=(nothing, nothing), io::Union{IO, Nothing}=nothing)
    nz, ne = sample_size
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)
    (mu, theta) = extra_params

    μ = zeros(S)
    Z = zeros(S)
    phi0 = zeros(N, C+1)
    phi  = @view phi0[:,2:end]
    pnc = zeros(N, C)
    twist = zeros(N, C)
    mgf = zeros(N)
    qnc = zeros(N, C)
    qn1 = @view qnc[:,1]
    u = zeros(N)
    W = zeros(N)
    w1 = @view weights[:,1]
    losses = zeros(N)

    innerlevel = InnerLevelTwisting(N, C)
    Ψ = init_Ψ()

    if mu == nothing
        outerlevel_twisting!(μ)
    else
        μ = mu
    end
    Zdist = MvNormal(μ, 1)

    estimate = 0
    estimates = zeros(ne)
    for i = 1:nz
        rand!(Zdist, Z)
        @. phi = normcdf((H - $(β*Z)) / denom)
        diff!(pnc, phi0; dims=2)
        # Twist bernoulli distribution pnc → qnc
        if theta == nothing
            twist!(innerlevel, pnc, weights, l)
            θ = get_result(innerlevel)
        else
            θ = theta
        end
        if θ != 0
            display(θ)
        end
        @. twist = pnc*ℯ^(θ*weights)
        sum!(mgf, twist)
        @. qnc = twist / mgf

        for j = 1:ne
            # Sample weights from bernoulli distribution
            rand!(u)
            @. W = (qn1 >= u)
            L = w1 ⋅ W
            lr = e^(-μ'Z + 0.5μ'μ -θ*L + Ψ(θ, pnc, weights))
            estimates[j] = (L >= l)*lr
            (i*ne + j) % 500 == 0 &&
                record_current(io, i, j, estimate, estimates)
        end
        estimate = (1-1/i)*estimate + (1/i)*mean(estimates)
    end
    return estimate
end
