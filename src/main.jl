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
import Main.CreditRisk: make_replications, plot_replications


n = 2500
c = 4
s = 5
l = 0.4

nz = 1000
ne = 1000

seed!(0)

for i in 1:10
    @time p = bernoulli_mc(Parameter(n,c,s,l), (nz, ne))
    print(p)
end

# make_replications((3, 3), "replications2.txt")
# plot_replications("replications2.txt")

# Profile.init(delay=0.1)
# @profile bernoulli_mc(Parameter(n,c,s,l), (nz, ne))
# Juno.profiletree()
# Juno.profiler()

end
