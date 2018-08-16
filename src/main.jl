include("CreditRisk.jl")

module tst
import Juno
import Serialization: serialize, deserialize
import Plots: plot, plot!
import Random: seed!
import Profile
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

nz = 2
ne = 2

seed!(0)


Profile.clear()
Profile.init(delay=0.01)
@profile bernoulli_mc(Parameter(n,c,s,l), (nz, ne))
Juno.profiletree()
Juno.profiler()
# print(p)

end
