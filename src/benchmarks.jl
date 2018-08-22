include("parameter.jl")
include("algorithm.jl")

import BenchmarkTools: @benchmark
import Optim: optimize
import Optim: BFGS, LBFGS, ConjugateGradient, GradientDescent
import Optim: MomentumGradientDescent, AcceleratedGradientDescent

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
    (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, w) = unpack(param)

    p = zeros(n,c)
    p[:,1] .= abs.(0.01 .* randn(n))
    p[:,2] .= 1 - p[:,1]

    objective(θ) = let θ = θ[1], (p,w,l) = (p,w,l)
        innerlevel_objective(θ, p, w, l)
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
