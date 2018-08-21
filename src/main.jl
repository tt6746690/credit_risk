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
import Main.CreditRisk: simple_mc, bernoulli_mc, glassermanli_mc
import Main.CreditRisk: get_estimates!
import Main.CreditRisk: make_replications, plot_replications

seed!(0)

n = 2500
c = 4
s = 5
l = 0.2
nz = 200
ne = 200

Profile.clear()
Profile.init(n=10^8, delay=0.01)
Juno.@profiler bernoulli_mc(Parameter(n,c,s,l), (nz, ne))
Juno.profiletree()
Juno.profiler()

end
