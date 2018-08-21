include("CreditRisk.jl")

module tst
import Juno
import Serialization: serialize, deserialize
import Plots: plot
import Random: seed!
import Profile
import Profile: @profile
import BenchmarkTools: @btime, @benchmark

import Main.CreditRisk: Parameter
import Main.CreditRisk: simple_mc, bernoulli_mc, glassermanli_mc
import Main.CreditRisk: get_estimates!
import Main.CreditRisk: make_replications, plot_replications, make_replications_b

seed!(0)

n = 2500
c = 4
s = 5
l = 0.5
nz = 2000
ne = 2000

θ = 2
μ = Vector(1:s)

make_replications_b((40, 30), "bernoulli2000.txt")

; @time p = bernoulli_mc(Parameter(n,c,s,l), (nz,ne))
; display(p)

# nrep = 10

# ls = range(0; stop=0.2, length=11)
# make_replications((ls, nrep), "gl1.txt")

# ls = range(0.22; stop=0.4, length=10)
# make_replications((ls, nrep), "gl2.txt")

# ls = range(0.42; stop=0.6, length=10)
# make_replications((ls, nrep), "gl3.txt")

# @time p = bernoulli_mc(Parameter(n,c,s,l), (nz,ne))
# display(p)  # 0.005
#
# @time p = glassermanli_mc(Parameter(n,c,s,l), (nz,ne), (μ, θ))
# display(p)


# Profile.clear()
# Profile.init(n=10^8, delay=0.01)
# Juno.@profiler glassermanli_mc(Parameter(n,c,s,l), (nz, ne))
# Juno.profiletree()
# Juno.profiler()

end
