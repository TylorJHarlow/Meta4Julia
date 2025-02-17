# Basic setup
__precompile__()
module Meta4Julia

# Call necessary packages
using LinearAlgebra, DataFrames, Distributions, Plots, Printf, BSplineKit, StatsModels

# Inputs & Outputs
import Base: show
export meta, predict, forest, funnel, s, modelcols, predFrame

# Define model type
abstract type model end

# Load Extra Packages
include("Models.jl")
include("Utils.jl")
include("Display.jl")

end