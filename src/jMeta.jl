# Basic setup
__precompile__()
module jMeta

# Call necessary packages
using LinearAlgebra, DataFrames, Distributions, Plots, Printf, BSplineKit, StatsModels

# Inputs & Outputs
import Base: show
export meta, predict, forest, funnel, s, modelcols, predFrame, @formula, DataFrame, rve

# Define model type
abstract type model end

# Load Extra Packages
include("Models.jl")
include("Utils.jl")
include("Display.jl")
include("Tests.jl")

end