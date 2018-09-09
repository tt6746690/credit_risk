include("parameter.jl")
include("algorithm.jl")

import BenchmarkTools: @benchmark
# import ForwardDiff: gradient, hessian
import Optim: optimize
import Optim: BFGS, LBFGS, ConjugateGradient, GradientDescent
import Optim: MomentumGradientDescent, AcceleratedGradientDescent
import PyPlot: surf, savefig
import Plots: contour, pdf

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

    Ψ = init_Ψ()
    function objective(θ)
        θ = θ[1]
        Ψ(θ, pnc, weights) - θ*l
    end

    function g(θ, w, l)
        θ = θ[1]
        sum(@. w * p * ℯ^(θ * w)) - l
    end

    x = [20]
    # println(gradient(objective, x))
    println(g(x, weights, l))

    # methods = [
    #     BFGS()
    #     LBFGS()
    #     ConjugateGradient()  # fasted, least allocation
    #     GradientDescent()
    #     MomentumGradientDescent()
    #     AcceleratedGradientDescent()
    # ]
    #
    # # println("Gradient free")
    # # for method in methods
    # #     @benchmark optimize(objective, [0.0], method = $method)
    # # end
    #
    # println("With g!")
    # for method in methods
    #     @benchmark optimize(objective, g!, [0.0], method = $method)
    # end

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

    for i in 1:3
        rand!(Zdist, Z)
        @. phi = normcdf((H - $(β*Z)) / denom)
        diff!(pnc, phi0; dims=2)
        println("pnc")
        display(pnc[1:20,1:2])

        function objective(θ)
            θ = θ[1]
            Ψ(θ, pnc, weights) - θ*l
        end

        # plot
        xs = -100:1000
        ys = [objective([x]) for x in xs]
        plot(xs, ys)
        savefig("$i.pdf", format="pdf")

        # optimization
        results = optimize(objective, [0.0], method = ConjugateGradient())
        if Optim.converged(results)
            θ = Optim.minimizer(results)[1]
            println("theta: $θ")
        else
            println("Did not converge")
        end
    end
end


" Visualize Outer level's objective function with many iterations of Z sampling
"
function outerlevel_objective()
    n = 2500
    c = 4
    s = 2  # 2D for easy visualization
    l = 0.2
    step = 0.01

    # for l in 0:0.01:0.2
    for l in [0.2]
        param = Parameter(n,c,s,l)
        (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(param)

        phi0 = zeros(N, C+1)
        phi  = @view phi0[:,2:end]
        pnc = zeros(N, C)

        Ψ = init_Ψ()
        innerlevel = InnerLevelTwisting(N, C)

        function objective(z)
            @. phi = normcdf((H - $(β*z)) / denom)
            diff!(pnc, phi0; dims=2)
            twist!(innerlevel, pnc, weights, l)
            θ = get_result(innerlevel)
            θ*l - Ψ(θ, pnc, weights) + 0.5z'z
        end

        xs = -2:step:2
        ys = -2:step:2
        zs = [objective([x, y]) for x in xs, y in ys]
        #
        obj_wrapped(x, y) = objective([x, y])
        c = contour(xs, ys, obj_wrapped, levels=100)
        pdf(c, "contour_$l.pdf")

        # display(surf(zs))
        # savefig("surf_$l.pdf", format="pdf")

    end
end


" Try different values of μ equidistance from the zero. See if a certain direction
yield no bias, lower variance, faster convergence "
function outerlevel_different_mu()
    n = 2500
    c = 4
    s = 2
    l = 0.2

    nz = 1000
    ne = 1000

    param = Parameter(n,c,s,l)
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(param)

    # with S=2, sample points from the circle with radius 2
    r = 2.5
    for theta = 0:π/8:2π
        x = r * cos(theta)
        y = r * sin(theta)
        mu = [x, y]


        phi0 = zeros(N, C+1)
        phi  = @view phi0[:,2:end]
        pnc = zeros(N, C)

        Ψ = init_Ψ()
        innerlevel = InnerLevelTwisting(N, C)

        function objective(z)
            @. phi = normcdf((H - $(β*z)) / denom)
            diff!(pnc, phi0; dims=2)
            twist!(innerlevel, pnc, weights, l)
            θ = get_result(innerlevel)
            θ*l - Ψ(θ, pnc, weights) + 0.5z'z
        end

        @time p = glassermanli_mc(Parameter(n,c,s,l), (nz,ne), (mu, nothing))
        println("mu: $mu   \t f(mu)= $(objective(mu))    \t p:$p")
    end

end
