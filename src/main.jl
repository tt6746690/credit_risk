include("CreditRisk.jl")

module tst
import Serialization: serialize, deserialize
import Plots: plot, plot!
import Random: seed!
import Profile: @profile
import BenchmarkTools: @btime, @benchmark

import Main.CreditRisk: Parameter
import Main.CreditRisk: simple_mc, bernoulli_mc
import Main.CreditRisk: get_estimates!
import Main.CreditRisk: make_replications


n = 2500
c = 4
s = 5
l = 0.2

nz = 200
ne = 200

seed!(0)

@benchmark p = bernoulli_mc(Parameter(n,c,s,l), (nz, ne))
# print(p)

end
