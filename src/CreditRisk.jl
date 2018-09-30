module CreditRisk

import Random: rand!
import Base.MathConstants: e
import LinearAlgebra: mul!, ⋅, dot
import Distributions: cdf, Normal, MvNormal
import Statistics: mean, quantile

import Serialization: serialize, deserialize
import BenchmarkTools: @benchmark

import Optim
import Optim: optimize
import Optim: BFGS, LBFGS, ConjugateGradient, GradientDescent
import Optim: MomentumGradientDescent, AcceleratedGradientDescent

import PyPlot: surf, savefig, plot
import Plots: contour, pdf, scatter

include("utils.jl")
include("parameter.jl")
include("algorithm.jl")
include("sampler.jl")
include("scripts.jl")

export msexpr, checksize, diff!, normcdf, invnormcdf
export Parameter, unpack
export InnerLevelTwisting, OuterLevelTwisting, twist!, init_Ψ, get_result, set_result!
export simple_mc, bernoulli_mc, glassermanli_mc, onelvl_mc
export slicesample
export make_replications, plot_replications, make_replications_b

end
