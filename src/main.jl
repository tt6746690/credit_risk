include("CreditRisk.jl")

module tst


import Juno
import Serialization: serialize, deserialize
import Random: seed!
import Profile
import Profile: @profile
import BenchmarkTools: @btime, @benchmark

import Optim: ConjugateGradient

import Main.CreditRisk: twist!, OuterLevelTwisting, get_result
import Main.CreditRisk: init_Ψ, InnerLevelTwisting, normcdf, unpack, diff!
using Main.CreditRisk

# seed!(0)

n = 2500
c = 4
s = 2
l = 1
nz = 1000
ne = 1000
param = Parameter(n,c,s,l)

# nrep = 10
# ls = range(0; stop=0.2, length=11)
# make_replications((ls, nrep), "gl1.txt")
# ls = range(0.22; stop=0.4, length=10)
# make_replications((ls, nrep), "gl2.txt")
# ls = range(0.42; stop=0.6, length=10)
# make_replications((ls, nrep), "gl3.txt")

# @time p = bernoulli_mc(Parameter(n,c,s,l), (nz,ne))
# display(p)  # 0.005
  # 4.781180 seconds (526.12 k allocations: 44.665 MiB, 0.61% gc time)

@time p = glassermanli_mc(Parameter(n,c,s,l), (nz,ne))
display(p)

# result = twist!(OuterLevelTwisting(), Parameter(n,c,s,l), ConjugateGradient())
# display(result)

# import Optim
# import Optim: optimize, ConjugateGradient
# using PyPlot
#
# parameter = Parameter(n,c,s,l)
# (N, C, S, l, cmm, ead, lgc, cn, β, H, denom, weights) = unpack(parameter)
# phi0 = zeros(N, C+1)
# phi  = @view phi0[:,2:end]
# pnc = zeros(N, C)
#
# Ψ = init_Ψ()
# innerlevel = InnerLevelTwisting(N, C)
#
# function objective(z)
#     @. phi = normcdf((H - $(β*z)) / denom)
#     diff!(pnc, phi0; dims=2)
#     twist!(innerlevel, pnc, weights, l)
#     θ = get_result(innerlevel)
#     θ*l - Ψ(θ, pnc, weights) + 0.5z'z
# end
#
# xs = -1:0.1:1
# ys = -1:0.1:1
# zs = [objective([x, y]) for x in xs, y in ys]
# surf(zs)
# savefig("outerlevelobj.pdf",format="pdf",dpi=1000);


# results = optimize(objective, ones(S), ConjugateGradient())
# display(results)

# if Optim.converged(results)
#     display("Converged!")
#     return Optim.minimizer(results)
# else
#     return zeros(S)
# end

end
