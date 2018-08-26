module CreditRisk

include("utils.jl")
include("parameter.jl")
include("algorithm.jl")
include("scripts.jl")

export get_estimates!, msexpr, checksize
export Parameter
export simple_mc, bernoulli_mc, glassermanli_mc
export make_replications, plot_replications, make_replications_b


end
