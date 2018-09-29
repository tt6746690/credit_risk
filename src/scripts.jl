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
        println(io, "l,glassermanli,mu")
        for l = ls
            for rep in 1:nrep

                parameter = Parameter(n,c,s,l)
                outerlevel = OuterLevelTwisting(n, c, s)
                twist!(outerlevel, parameter)
                μ = get_result(outerlevel)

                @time p = glassermanli_mc(parameter, (nz, ne), (μ, nothing))
                println(io, string(l) * "," * string(p) * "," * string(μ))
                flush(io)
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
