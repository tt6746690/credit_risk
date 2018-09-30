using Test
using CreditRisk

import Distributions: MvNormal, Normal, pdf
import Statistics: mean, std


@testset "slice sampler" begin

    @testset "arguments" begin
        for d in [1, 2]
            x0 = zeros(d)
            g(x) = pdf(MvNormal(d, 1), x)

            for n in convert(Array{Int64}, [floor(x) for x in 10*rand(5)])
                m = 10
                samples = slicesample(x0, g, n; step_limit=m, burn=2, thin=1)
                @test size(samples) == (d, n)
            end

            for m in convert(Array{Int64}, [floor(x) for x in 10*rand(5)])
                n = 2
                samples = slicesample(x0, g, n; step_limit=m, burn=2, thin=1)
                @test size(samples) == (d, n)
            end
        end
    end

    @testset "correctness" begin
        n = 5000

        x0 = zeros(1)
        g(x) = pdf(Normal(), x[1])
        samples = slicesample(x0, g, n; step_limit=10)
        @test mean(samples) ≈ 0. atol = 0.1
        @test std(samples)  ≈ 1. atol = 0.1

        x0 = zeros(2)
        g(x) = pdf(MvNormal(2, 1), x)
        samples = slicesample(x0, g, n; step_limit=10)
        @test mean(samples) ≈ 0. atol = 0.1
        @test std(samples) ≈ 1.  atol = 0.1
    end
end

@testset "Computing Ψ" begin
    n, c, s, l = 10, 2, 2, 1
    param = Parameter(n,c,s,l)
    w = param.weights
    p = fill(0.02, n, c)
    p[:,2] = 1 .- p[:,1]

    memefficientΨ = CreditRisk.init_Ψ()
    expectedΨ(θ, p, w) = sum(log.(sum(@. p*ℯ^(θ*w); dims=2)))

    for θ in 0:2:10
        @test expectedΨ(θ, p, w) ≈ memefficientΨ(θ, p, w)
    end
end

@testset "Twisting a bernoulli distribution" begin
    n, c, s, l = 10, 2, 2, 1
    param = Parameter(n,c,s,l)
    w = param.weights
    p = fill(0.02, n, c)
    p[:,2] = 1 .- p[:,1]

    binary_twisted(p, θ, w) = begin
        out = @. p*ℯ^(θ*w) / (1 + p*(ℯ^(θ*w) - 1))
        out[:,2] = 1 .- out[:,1]
        out
    end
    multistate_twisted(p, θ, w) = @. p*ℯ^(θ*w) / $sum(@. p*ℯ^(θ*w); dims=2)

    for θ in 0:2:10
        @test binary_twisted(p, θ, w) ≈ multistate_twisted(p, θ, w)
    end

    multistate_twisted!(twist, mgf, q, p, θ, w) = begin
        @. twist = p*ℯ^(θ*w)
        sum!(mgf, twist)
        @. q = twist / mgf
    end

    twist = zeros(n,c)
    mgf = zeros(n)
    q = zeros(n,c)

    for θ in 0:2:10
        multistate_twisted!(twist, mgf, q, p, θ, w)
        @test multistate_twisted(p, θ, w) ≈ q
    end
end
