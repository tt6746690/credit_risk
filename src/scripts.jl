import Serialization: serialize

include("parameter.jl")

function make_replications((nl, replication_per_l), filename::String)
    n = 2500
    c = 4
    s = 5
    l = 0.2

    nz = 1000
    ne = 1000

    open(filename, "w") do io
        lp = []
        for l = range(0; stop=0.8, length=nl)
            for _ in 1:replication_per_l
                @time p = bernoulli_mc(Parameter(n,c,s,l), (nz, ne))
                append!(lp, (l, p))
            end
        end
        serialize(io, lp)
    end
end
