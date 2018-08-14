include("CreditRisk.jl")

module tst
import Random: seed!
import Profile: @profile
import BenchmarkTools: @btime, @benchmark

import Main.CreditRisk: Parameter, simple_mc, bernoulli_mc

n = 2500
c = 4
s = 5
l = 0.2

nz = 10000
ne = 10000

seed!(0)

# open("simple_mc", "w") do io
#     @time p = simple_mc(Parameter(n,c,s,l), (nz, ne), io)
#     print(p)
# end

open("bernoulli_mc", "w") do io
    @time p = bernoulli_mc(Parameter(n,c,s,l), (nz, ne), io)
    print(p)
end

end
