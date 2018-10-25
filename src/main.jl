import Pkg
Pkg.activate(".")

include("CreditRisk.jl")

module tst

using Revise
using ArgParse
import CreditRisk


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--n"
            arg_type = Int64
            default = 2500
        "--c"
            arg_type = Int64
            default = 4
        "--s"
            arg_type = Int64
            default = 5
        "--l"
            arg_type = Float64
            default = 0.2
        "--nz"
            arg_type = Int64
            default = 100
        "--ne"
            arg_type = Int64
            default = 100
        "--n_init"
            arg_type = Int64
            default = 5
        "--filename"
            default = "./default.txt"
        "--a"
            default = "bernoulli"
    end

    return parse_args(s)
end

args = parse_commandline()
n = args["n"]
c = args["c"]
s = args["s"]
l = args["l"]
nz = args["nz"]
ne = args["ne"]
n_init = args["n_init"]
algo = args["a"]
filename = args["filename"]
param = CreditRisk.Parameter(n,c,s,l)

println("Parsed args:")
for (arg,val) in args
    println("  $arg  =>  $val")
end

open(filename, "w") do io
    if algo == "bernoulli"
        estimate = CreditRisk.bernoulli_mc(param, (nz, ne), io)
        println(estimate)
    elseif algo == "glassermanli"
        estimate = CreditRisk.glassermanli_mc(param, (nz, ne), (nothing, nothing, n_init), io)
        println(estimate)
    else
        estimate = CreditRisk.bernoulli_mc(param, (nz, ne), io)
        println(estimate)
    end
end

end
