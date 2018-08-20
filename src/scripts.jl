import Plots: scatter
import Serialization: serialize, deserialize

include("parameter.jl")

function make_replications((nl, replication_per_l), filename::String)
    n = 2500
    c = 4
    s = 5

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


function plot_replications(filename::String)
    open(filename, "r") do io
        xys = deserialize(io)
        print(xys)
        xs, ys = zip(xys...)
        plt = scatter(xs, ys;
            ylims=(-8, 0.2),
            ylabel="log MC estimates",
            xlabel="tail l",
            plot_title="log (bernoulli) MC estimates vs tail l")
        display(plt)
    end
end
