include("utils.jl")
include("parameter.jl")

function record_current(io, i, j, estimate, estimates)
    if io != nothing
        prop = j/((i-1)*ne+j)
        est = (1-prop)*estimate + (prop)*mean(estimates[1:j])
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
    θ::Array{Float64, 1}
    Ψ::Any
    initialguess::Array{Float64, 1}
    n_restarts::Int64

    function InnerLevelTwisting(N, C, wp, θ, Ψ, initialguess, n_restarts)
        @checksize (N, C)   wp
        @checksize (1,)     θ
        @checksize (1,)     initialguess
        new(N, C, wp, θ, Ψ, initialguess, n_restarts)
    end
end

function InnerLevelTwisting(N, C)
    wp = zeros(N, C)
    θ = [0]
    Ψ = init_Ψ()
    initialguess = [0]
    n_restarts = 5
    InnerLevelTwisting(N, C, wp, θ, Ψ, initialguess, n_restarts)
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
        optimizer = ConjugateGradient()
        for i = 1:t.n_restarts
            results = optimize(objective, t.initialguess, optimizer)
            if Optim.converged(results)
                minimizer = Optim.minimizer(results)
                set_result!(t, minimizer[1])
                t.initialguess[:] .= minimizer
                return
            end
            optimizer = BFGS()
        end

        display("Inner level optimization failed to converged with $(t.n_restarts) restarts.")
        set_result!(t, 0)
    else
        set_result!(t, 0)
    end
end


struct OuterLevelTwisting
    N::Int64
    C::Int64
    S::Int64
    μ::Array{Float64, 1}
    innerlevel::InnerLevelTwisting
    Ψ::Any
    n_restarts::Int64

    function OuterLevelTwisting(N, C, S, μ, innerlevel, Ψ, n_restarts)
        @checksize (S,)     μ
        new(N, C, S, μ, innerlevel, Ψ, n_restarts)
    end
end

function OuterLevelTwisting(N, C, S)
    μ = zeros(S)
    innerlevel = InnerLevelTwisting(N, C)
    Ψ = init_Ψ()
    #initialguess = zeros(S)
    #n_restarts = 20
    n_restarts = 2
    OuterLevelTwisting(N, C, S, μ, innerlevel, Ψ, n_restarts)
end


get_result(t::OuterLevelTwisting) = t.μ
set_result!(t::OuterLevelTwisting, μ) = (t.μ[:] = μ)


function twist_local!(t::OuterLevelTwisting, parameter::Parameter)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)
    phi0 = zeros(N, C+1)
    phi  = @view phi0[:,2:end]
    pnc = zeros(N, C)

    function objective(z)
        @checksize (S,) z
        @. phi = normcdf((H - $(β*z)) / denom)
        diff!(pnc, phi0; dims=2)
        twist!(t.innerlevel, pnc, weights, l)
        θ = get_result(t.innerlevel)
        θ*l - t.Ψ(θ, pnc, weights) + 0.5z'z
    end

    optimizers = Dict(
        "BFGS" => BFGS(),
        "LBFGS" => LBFGS(),
        "MomentumGradientDescent" => MomentumGradientDescent()
    )
    minimizers = []
    for i = 1:t.n_restarts
        initialguess = rand(S)*2.0 .- 1.0
        op = rand(optimizers)
        optimizer = op[2]
        println("Start $i optimization with $(op[1])")

        results = optimize(objective, initialguess, optimizer)
        if Optim.converged(results)
            minimizer = Optim.minimizer(results)
            #println("mu: $minimizer; f(mu): $(objective(minimizer))")
            push!(minimizers, minimizer)
        end
    end
    if size(minimizers) == 0
        set_result!(t, zeros(S))
    else
        fmins = [objective(x) for x in minimizers]
        # find index i wheere fmin is smallest
        println("minimizers: $fmins")
        #i = argmin(fmins)
        set_result!(t, minimizers[1])
    end
end


function twist_global!(t::OuterLevelTwisting, parameter::Parameter)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)
    phi0 = zeros(N, C+1)
    phi  = @view phi0[:,2:end]
    pnc = zeros(N, C)

    function objective(z)
        @checksize (S,) z
        @. phi = normcdf((H - $(β*z)) / denom)
        diff!(pnc, phi0; dims=2)
        twist!(t.innerlevel, pnc, weights, l)
        θ = get_result(t.innerlevel)
        θ*l - t.Ψ(θ, pnc, weights) + 0.5z'z
    end

    optimizers = Dict(
        "BFGS" => BFGS(),
        "LBFGS" => LBFGS(),
        "MomentumGradientDescent" => MomentumGradientDescent()
    )
    minimizers = []
    for i = 1:t.n_restarts
        initialguess = rand(S)*2.0 .- 1.0
        op = rand(optimizers)
        optimizer = op[2]
        println("Start $i optimization with $(op[1])")

        results = optimize(objective, initialguess, optimizer)
        if Optim.converged(results)
            minimizer = Optim.minimizer(results)
            #println("mu: $minimizer; f(mu): $(objective(minimizer))")
            push!(minimizers, minimizer)
        end
    end
    if size(minimizers) == 0
        set_result!(t, zeros(S))
    else
        fmins = [objective(x) for x in minimizers]
        # find index i wheere fmin is smallest
        println("minimizers: $fmins")
        i = argmin(fmins)
        set_result!(t, minimizers[i])
    end
end


"""
    glassermanli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::IO=Union{IO, Nothing}=nothing)

    Monte Carlo Sampling for computing default probability: P(L ⩾ l), with importance sampling
        of Z and weights w
       ⋅ Z ∼ N(μ, I) where μ is shifted mean that minimizes variance
       ⋅ W ∼ q       where q is shifted bernoulli probability that minimizes upper bound on the variance
"""
function glassermanli_mc_local(parameter::Parameter, sample_size::Tuple{Int64, Int64}, filename:: String, extra_params=(nothing, nothing))
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

    outerlevel = OuterLevelTwisting(N, C, S)
    innerlevel = InnerLevelTwisting(N, C)
    Ψ = init_Ψ()

    if mu == nothing
        twist_local!(outerlevel, parameter)
        μ = get_result(outerlevel)
    else
        μ = mu
    end
    Zdist = MvNormal(μ, 1)

    estimate = 0
    estimates = zeros(ne)
    open(filename, "w") do io
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
                (i*ne + j) % 5 == 0 && begin
                    prop = j/((i-1)*ne+j)
                    est = (1-prop)*estimate + (prop)*mean(estimates[1:j])
                    println(io, string(mean(est)))
                    flush(io)
                end
            end
            estimate = (1-1/i)*estimate + (1/i)*mean(estimates)
        end
    end
    return estimate
end


"""
    glassermanli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::IO=Union{IO, Nothing}=nothing)

    Monte Carlo Sampling for computing default probability: P(L ⩾ l), with importance sampling
        of Z and weights w
       ⋅ Z ∼ N(μ, I) where μ is shifted mean that minimizes variance
       ⋅ W ∼ q       where q is shifted bernoulli probability that minimizes upper bound on the variance
"""
function glassermanli_mc_global(parameter::Parameter, sample_size::Tuple{Int64, Int64}, filename:: String, extra_params=(nothing, nothing))
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

    outerlevel = OuterLevelTwisting(N, C, S)
    innerlevel = InnerLevelTwisting(N, C)
    Ψ = init_Ψ()

    if mu == nothing
        twist_global!(outerlevel, parameter)
        μ = get_result(outerlevel)
    else
        μ = mu
    end
    Zdist = MvNormal(μ, 1)

    estimate = 0
    estimates = zeros(ne)
    open(filename, "w") do io
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
                (i*ne + j) % 5 == 0 && begin
                    prop = j/((i-1)*ne+j)
                    est = (1-prop)*estimate + (prop)*mean(estimates[1:j])
                    println(io, string(mean(est)))
                    flush(io)
                end
            end
            estimate = (1-1/i)*estimate + (1/i)*mean(estimates)
        end
    end
    return estimate
end





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
    onelvl_mc(parameter::Parameter, nz::Int64, io::IO=Union{IO, Nothing}=nothing)

    One level Monte Carlo simulation using
        conditional portfolio loss approaches in distribution to standard normal distribution
        (L - μ) / σ  ⟶ N(0, 1)
"""
function onelvl_mc(parameter::Parameter, nz::Int64, io::Union{IO, Nothing}=nothing)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)

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

        (i) % 500 == 0 && begin
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

    Zs = slicesample(initial_point, π, nz;
        step_limit=20,
        width=0.5,
        burn=Int(burn_ratio * nz),
        thin=3)

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

        (i) % 500 == 0 && begin
            if io != nothing
                println(io, string(mean(estimates[1:i])))
                flush(io)
            end
        end
    end
    return mean(estimates)
end
