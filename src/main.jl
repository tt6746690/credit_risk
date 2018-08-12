import Profile: @profile
import BenchmarkTools: @btime, @benchmark

include("CreditRisk.jl")

module tst
import Main.CreditRisk: Parameter, simple_mc

n = 20
c = 4
s = 10
l = 0.2

nz = 100
ne = 100

p = simple_mc(p, (nz, ne))
print(p)

end
