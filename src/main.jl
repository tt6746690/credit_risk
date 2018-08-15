include("CreditRisk.jl")

module tst
import Plots: plot, plot!
import Random: seed!
import Profile: @profile
import BenchmarkTools: @btime, @benchmark

import Main.CreditRisk: Parameter
import Main.CreditRisk: simple_mc, bernoulli_mc
import Main.CreditRisk: get_estimates!


n = 2500
c = 4
s = 5
l = 0.2

nz = 1000
ne = 1000

seed!(0)
#
# open("/dev/null", "w") do io
#     @time p = simple_mc(Parameter(n,c,s,l), (nz, ne), io)
#     print(p)
# end
#
# nl = 50
# replication_per_l = 30
# 50 * 30 * 20s ≈ 8.33hr to finish

nl = 2
replication_per_l = 1
devnull = open("/dev/null", "w")

open("replications.txt", "w") do io
    for l = range(0; stop=0.8, length=nl)
        for _ in 1:replication_per_l
            @time p = bernoulli_mc(Parameter(n,c,s,l), (nz, ne), devnull)
            println(io, string(l) * "," * string(p))
            flush(io)
        end
    end
end

close(devnull)

# open("/dev/null", "w") do io
#     @time p = bernoulli_mc(Parameter(n,c,s,l), (nz, ne), io)
#     print(p)
# end

#
# mc = Vector{Float64}()
# get_estimates!(mc, "bernoulli_mc")
# plt = plot(mc[1:50000])
# display(plt)


end
