module CreditRisk

import Random: rand!
import Base.MathConstants: e
import LinearAlgebra: mul!, ⋅, dot
import Distributions: pdf, cdf, Normal, MvNormal, TDist
import Statistics: mean, quantile, var

import Serialization: serialize, deserialize
import BenchmarkTools: @benchmark

import Optim
import Optim: optimize
import Optim: BFGS, LBFGS, ConjugateGradient, GradientDescent
import Optim: MomentumGradientDescent, AcceleratedGradientDescent

import PyPlot: surf, savefig, plot
import Plots: contour, scatter

# General utilities
include("utils.jl")
include("sampler.jl")

# problem setup & helpers
include("parameter.jl")
include("confidence_interval.jl")

# different MC algorithms
include("mc_simple.jl")
include("mc_bernoulli.jl")
include("mc_glassermanli.jl")
include("mc_zero_variance.jl")

include("scripts.jl")

export msexpr, checksize, diff!, normcdf, invnormcdf, invtcdf
export Parameter, unpack
export slicesample
export InnerLevelTwisting, OuterLevelTwisting, twist!, init_Ψ, get_result, set_result!
export simple_mc, bernoulli_mc, glassermanli_mc, onelvl_mc, onelvlISCLT_mc
export make_replications, plot_replications, make_replications_b

end
