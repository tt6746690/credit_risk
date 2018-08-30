include("CreditRisk.jl")

module tst


import Juno
import Serialization: serialize, deserialize
import Random: seed!
import Profile
import Profile: @profile
import BenchmarkTools: @btime, @benchmark
import Optim: ConjugateGradient

using Main.CreditRisk

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

# @time p = glassermanli_mc(Parameter(n,c,s,l), (nz,ne))
# display(p)

# result = twist!(OuterLevelTwisting(), Parameter(n,c,s,l), ConjugateGradient())
# display(result)

outerlevel_objective()

end
