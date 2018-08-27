include("parameter.jl")
include("algorithm.jl")

import BenchmarkTools: @benchmark
import Optim: optimize
import Optim: BFGS, LBFGS, ConjugateGradient, GradientDescent
import Optim: MomentumGradientDescent, AcceleratedGradientDescent
import PyPlot: plot, surf, savefig

" Inner level optimization on a univariate convex function
    Should almost always reach global min regardless of algorithm used
    So just find one that runs fastest
"
function innerlevel_optimizer()
    n = 2500
    c = 4
    s = 5
    l = 0.5
    nz = 1000
    ne = 1000
    param = Parameter(n,c,s,l)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(param)

    p = zeros(n,c)
    p[:,1] .= abs.(0.01 .* randn(n))
    p[:,2] .= 1 - p[:,1]

    function objective(θ)
        θ = θ[1]
        Ψ(θ, pnc, weights) - θ*l
    end

    methods = [
        BFGS()
        LBFGS()
        ConjugateGradient()  # fasted, least allocation
        GradientDescent()
        MomentumGradientDescent()
        AcceleratedGradientDescent()
    ]

    for method in methods
        @benchmark optimize(objective, [0.0], method = $method)
    end

end



" Visualize Inner level's objective function with many iterations of Z sampling
"
function innerlevel_objective()
    n = 2500
    c = 4
    s = 5
    l = 0.2
    param = Parameter(n,c,s,l)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(param)

    Ψ = init_Ψ()

    Zdist = MvNormal(S, 1)
    Z = zeros(S)
    pnc = zeros(N, C)
    phi0 = zeros(N, C+1)
    phi  = @view phi0[:,2:end]
    twist = zeros(N, C)
    mgf = zeros(N)
    qnc = zeros(N, C)

    for i in 1:1
        rand!(Zdist, Z)
        @. phi = normcdf((H - $(β*Z)) / denom)
        diff!(pnc, phi0; dims=2)
        println("pnc")
        display(pnc[1:20,1:2])

        function objective(θ)
            θ = θ[1]
            a = Ψ(θ, pnc, weights) - θ*l
            # b = sum(log.(sum(@. pnc*ℯ^(θ*weights); dims=2)))
            # @assert a == b
            return a
        end

        # plot
        xs = -100:1000
        ys = [objective([x]) for x in xs]
        display(Z)
        plot(xs, ys)
        savefig("$i.pdf", format="pdf")

        # optimization
        println("Criteria: $(sum(weights .* pnc) > l)")
        results = optimize(objective, [0.0], method = ConjugateGradient())

        if Optim.converged(results)
            θ = Optim.minimizer(results)[1]
            println("theta: $θ")
        else
            θ = 0
        end

        # θ = 300

        # twisted distribution
        # @. twist = pnc*ℯ^(θ*weights)
        # display(twist[1:20,1:2])
        # sum!(mgf, twist)
        # display(mgf[1:20])
        # @. qnc = twist / mgf
        #
        #
        # println("qnc")
        # display(qnc[1:20,1:2])
    end
end


innerlevel_objective()
