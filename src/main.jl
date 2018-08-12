import Profile: @profile
import BenchmarkTools: @btime, @benchmark

include("CreditRisk.jl")

module tst
import Main.CreditRisk: Parameter, simple_mc

n = 2500
c = 4
s = 5
l = 0.2

nz = 1000
ne = 1000

@time p = simple_mc(Parameter(n,c,s,l), (nz, ne))
print(p)

end
