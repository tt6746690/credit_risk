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
nz = 1000
ne = 1000

# @time p = bernoulli_mc(Parameter(n,c,s,l), (nz,ne))
# display(p)  # 0.005

@time p = glassermanli_mc(Parameter(n,c,s,l), (nz,ne))
display(p)


# Profile.clear()
# Profile.init(n=10^8, delay=0.01)
# Juno.@profiler glassermanli_mc(Parameter(n,c,s,l), (nz, ne))
# Juno.profiletree()
# Juno.profiler()

end
