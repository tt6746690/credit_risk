import Plots: scatter
import Serialization: serialize, deserialize

include("parameter.jl")


function make_replications_b((nl, nrep), filename::String)
    n = 2500
    c = 4
    s = 5
    nz = 2000
    ne = 2000

    open(filename, "w") do io
        println(io, "l,bernoulli")
        for l = range(0; stop=0.8, length=nl)
            for rep in 1:nrep
                @time p = bernoulli_mc(Parameter(n,c,s,l), (nz,ne))
                println(io, string(l) * "," * string(p))
                flush(io)
            end
        end
    end
end

function make_replications((ls, nrep), filename::String)
    n = 2500
    c = 4
    s = 5
    nz = 1000
    ne = 1000

    open(filename, "w") do io
        println(io, "l,glassermanli,mu,theta")
        for l = ls
            for rep in 1:nrep
                for μ in range(0; stop=2, length=2)
                    θ = 0
                    @time p2 = glassermanli_mc(Parameter(n,c,s,l), (nz, ne), (fill(μ, s), θ))
                    println(io, string(l) * "," * string(p2) * "," * string(μ) * "," * string(θ))
                    flush(io)
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
