
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
    # Number of restarts before giving up
    n_restarts::Int64
    # Number of initial guesses to start optimization
    n_initial::Int64

    function OuterLevelTwisting(N, C, S, μ, innerlevel, Ψ, n_restarts, n_initial)
        @checksize (S,)     μ
        new(N, C, S, μ, innerlevel, Ψ, n_restarts, n_initial)
    end
end

""" n_initial=1 -> local;  n_initial>1 -> global """
function OuterLevelTwisting(N, C, S, n_initial)
    μ = zeros(S)
    innerlevel = InnerLevelTwisting(N, C)
    Ψ = init_Ψ()
    n_restarts = 5
    OuterLevelTwisting(N, C, S, μ, innerlevel, Ψ, n_restarts, n_initial)
end


get_result(t::OuterLevelTwisting) = t.μ
set_result!(t::OuterLevelTwisting, μ) = t.μ[:] = isa(μ, Array) ? μ : [μ]


function twist!(t::OuterLevelTwisting, parameter::Parameter)
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
    while length(minimizers) < t.n_initial
        for i = 1:t.n_restarts
            initialguess = rand(S)*2.0 .- 1.0
            op = rand(optimizers)
            optimizer = op[2]
            println("Start $i optimization with $(op[1])")

            results = optimize(objective, initialguess, optimizer)
            if Optim.converged(results)
                minimizer = Optim.minimizer(results)
                push!(minimizers, minimizer)
                break
            elseif i == t.n_restarts
                push!(minimizers, zeros(S))
            end
        end
    end

    fmins = [objective(x) for x in minimizers]
    println("minimizers: $minimizers")
    println("obj(minimizers): $fmins")
    set_result!(t, minimizers[argmin(fmins)])
end


"""
    glassermanli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, io::IO=Union{IO, Nothing}=nothing)

    Monte Carlo Sampling for computing default probability: P(L ⩾ l), with importance sampling
        of Z and weights w
       ⋅ Z ∼ N(μ, I) where μ is shifted mean that minimizes variance
       ⋅ W ∼ q       where q is shifted bernoulli probability that minimizes upper bound on the variance

    Access intermediate estimates

        # 1. create IOBuffer
        io = IOBuffer()

        # 2. feed io to the argument
        estimate = CreditRisk.glassermanli_mc(param, (100,100), (zeros(s), 0), io)

        # 3. get estimate as usual
        println(estimate)

        # 4. the intermediate estimates also available
        estimates = String(take!(io))
        println("---")
        println(estimates[end-100:end])

"""
function glassermanli_mc(parameter::Parameter, sample_size::Tuple{Int64, Int64}, extra_params=(nothing, nothing, nothing), io::Union{IO, Nothing}=nothing)
    nz, ne = sample_size
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)
    (mu, theta, n_initial) = extra_params

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

    outerlevel = OuterLevelTwisting(N, C, S, n_initial)
    innerlevel = InnerLevelTwisting(N, C)
    Ψ = init_Ψ()

    if mu == nothing
        twist!(outerlevel, parameter)
        μ = get_result(outerlevel)
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

            k = (i-1)*ne+j
            if k % 500 == 0
                record_current(io, i, j, k, estimate, estimates)
            end

        end
        estimate = (1-1/i)*estimate + (1/i)*mean(estimates)
    end
    return estimate
end
