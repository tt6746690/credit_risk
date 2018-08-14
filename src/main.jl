include("CreditRisk.jl")

module tst
import Profile
import Profile: @profile
import BenchmarkTools: @btime, @benchmark

import Main.CreditRisk: Parameter, simple_mc

n = 2500
c = 4
s = 5
l = 0.2

nz = 500
ne = 500

@time p = simple_mc(Parameter(n,c,s,l), (nz, ne))
print(p)

end
