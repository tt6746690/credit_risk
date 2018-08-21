import Plots: scatter
import Serialization: serialize, deserialize

include("parameter.jl")

function make_replications((nl, nrep), filename::String)
    n = 2500
    c = 4
    s = 5
    nz = 1000
    ne = 1000

    open(filename, "w") do io
        println(io, "l,bernoulli,glassermanli,mu,theta")
        for l = range(0; stop=1, length=nl)
            for rep in 1:nrep
                @time p1 = bernoulli_mc(Parameter(n,c,s,l), (nz, ne))
                for μ in range(0; stop=1.5, length=4)
                    for θ in range(0; stop=1.5, length=4)
                        @time p2 = glassermanli_mc(Parameter(n,c,s,l), (nz, ne), (fill(μ, s), θ))
                        println(io, string(l) * "," * string(p1) * "," * string(p2) * "," * string(μ) * "," * string(θ))
                        flush(io)
                    end
                end
            end
        end
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
