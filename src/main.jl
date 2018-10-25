import Pkg
Pkg.activate(".")

include("CreditRisk.jl")

module tst

using Revise

import CreditRisk
import Serialization: serialize
import Profile: @profile
import BenchmarkTools: @btime, @benchmark


n = 2500
c = 4
s = 5
l = 0.2
nz = 10000
ne = 10000
param = CreditRisk.Parameter(n,c,s,l)

open("gl_long.txt", "w") do io
    @time estimates = CreditRisk.glassermanli_mc(param, (nz, ne))
    serialize(io, estimates)
end


# @time p = onelvl_mc(param, 10000)
# display(p)

# nrep = 10
# ls = range(0; stop=0.2, length=11)
# make_replications((ls, nrep), "glassermanli1.txt")
# ls = range(0.22; stop=0.4, length=10)
# make_replications((ls, nrep), "glassermanli2.txt")
# ls = range(0.42; stop=0.6, length=10)
# make_replications((ls, nrep), "glassermanli3.txt")
# ls = range(0.62; stop=0.8, length=10)
# make_replications((ls, nrep), "glassermanli4.txt")

end
