using Test
using CreditRisk


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
        display(multistate_twisted(p, θ, w))
        @test multistate_twisted(p, θ, w) ≈ q
    end
end
