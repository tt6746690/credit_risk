import Profile: @profile
import BenchmarkTools: @btime, @benchmark

include("CreditRisk.jl")

module tst
import Main.CreditRisk: Parameter, simple_mc

n = 20
c = 4
s = 10
l = 0.2

p = Parameter(n,c,s,l)

nz = 10
ne = 10

simple_mc(p, (nz, ne))

end
