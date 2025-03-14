
# Define Basic Model
struct mdl <: model
    df::DataFrame
    Type::String
    Mu::Float64
    Weights::Vector{Float64}
    V::Float64
    CI::Vector{Float64}
    PI::Vector{Float64}
    SE::Vector{Float64}
    Q::Float64
    Qp::Float64
    Tau2::Float64
    I2::Float64
end

# Define Meta regression model
mutable struct MRmdl <: model
    df::DataFrame
    Type::String
    Formula::FormulaTerm
    Mu::Float64
    Weights::Vector{Float64}
    V::Float64
    CI::Vector{Float64}
    PI::Vector{Float64}
    SE::Vector{Float64}
    Beta::Vector{Float64}
    CovBeta::Matrix{Float64}
    SEBeta::Vector{Float64}
    Q::Float64
    Qp::Float64
    Tau2::Float64
    I2::Float64
    Clusters::Symbol
end

# This main function of this package offers to main methods of choice. The first is a basic 
# fixed effects meta analysis. Functionality for this method is limitied, in part to deter its 
# use. It does not support either subgroup/moderator analyses or metaregression.
#
# The second method available accepts a given formula and solves using restricted maximum 
# liklihood estimation. If desired, the maximum number of iterations  and tolerance for the 
# REML can be specified after any moderators of interest. REML method is taken from Veronik 
# et al 2015:  https://doi.org/10.1002/jrsm.1164
#
# Formulas are interpreted through the nomenclature of the StatsModels package:
# https://juliastats.org/StatsModels.jl/stable/
#
# Additional functions m() and s() are provided such that categorical moderators and spline
# regression can be incorporated in the meta-regression approach.
#
# Written Tylor J Harlow 2/6/2025


# Meta analysis function
# Basic Fixed effects version
function meta(df::DataFrame ; α::Float64=0.05, d=:d, v=:v)

    # Unpack
    d = df[!,d]
    v = df[!,v]

    # Basic Fixed effects meta
    w = 1 ./ v
    μ = sum(w .* d) ./ sum(w)

    # Q-statistic
    Q = sum(w .* (d .- μ) .^ 2)

    # Evaluate Chisq
    Qp = cdf(Chisq(k - p), Q)

    # Get early estimate τ2
    k = length(d)
    c = sum(w) - (sum(w .^ 2) / sum(w))
    τ2 = max(0.0, (Q - (k - 1)) / c)

    # Get confidence intervals
    z = quantile(Normal(0, 1), 1 - α / 2)
    V = 1 ./ sum(w)
    CI = [μ - (z * sqrt(V)), μ + (z * sqrt(V))]
    PI = [μ - (z * sqrt(V + τ2)), μ + (z * sqrt(V + τ2))]

    # Get I2
    I2 = 100 * (Q - (k - 1)) ./ Q

    # Send out
    return mdl(df, raw"Fixed Effects Meta-Analysis", μ, w, V, CI, PI, sqrt.(v), Q, Qp, τ2, I2)
end

# Larger, meta-regression approach
function meta(df::DataFrame,formula::FormulaTerm; v = :v, clusters = :clusters, α::Float64=0.05, iter::Int=100, tol::Float64=1e-5)

    # Check standard errors vs. variance
    v = df[!,v]

    # Call Meta Regression
    d, v, w, μ, Q, Qp, τ2, I2, β, covβ, seβ = reml(df,v,formula,iter,tol)

    # Get confidence intervals
    z = quantile(Normal(0, 1), 1 - α / 2)
    V = 1 ./ sum(w)
    CI = [μ - (z * sqrt(V)), μ + (z * sqrt(V))]
    PI = [μ - (z * sqrt(V + τ2)), μ + (z * sqrt(V + τ2))]

    # Send out
    return MRmdl(df, raw"Meta-Regression", formula, μ, w, V, CI, PI, sqrt.(v), β, covβ, seβ, Q, Qp, τ2, I2, clusters)

end
