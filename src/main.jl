import Pkg
Pkg.activate(".")

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
s = 5
l = 0.2
nz = 1000
ne = 1000
param = Parameter(n,c,s,l)

@time p = onelvl_mc(param, 100000)
display(p)


# nrep = 10
# ls = range(0; stop=0.2, length=11)
# make_replications((ls, nrep), "glassermanli1.txt")
# ls = range(0.22; stop=0.4, length=10)
# make_replications((ls, nrep), "glassermanli2.txt")
# ls = range(0.42; stop=0.6, length=10)
# make_replications((ls, nrep), "glassermanli3.txt")
# ls = range(0.62; stop=0.8, length=10)
# make_replications((ls, nrep), "glassermanli4.txt")


# display("bernoulli")
# for n = 1:10
#     p = bernoulli_mc(Parameter(n,c,s,l), (nz,ne))
#     display(p)
# end
#
# display("one level")
# for n = 1:10
#     p = onelvl_mc(param, ne*nz)
#     display(p)
# end

# 0.009933000000000006
# 0.01987399999999998
# 0.02910700000000006
# 0.03111500000000009
# 0.047101000000000066
# 0.03835700000000001
# 0.0319319999999999
# 0.04774800000000004
# 0.07832600000000016
# 0.07708899999999999


# bernoulli
# 0.009560999999999991
# 0.020446999999999976
# 0.03022600000000003
# 0.038724999999999996
# 0.045941
# 0.03898599999999998
# 0.031371000000000045
# 0.04926800000000004
# 0.0776170000000002
# 0.07477800000000003



# @time p = glassermanli_mc(Parameter(n,c,s,l), (nz,ne))
# display(p)


end
