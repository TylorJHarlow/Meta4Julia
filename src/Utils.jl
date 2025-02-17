

# Run Restricted maximum liklihood estiamtion on meta regression model
function reml(df::DataFrame, v::Vector{Float64}, type::String, formula::FormulaTerm, iter::Int, tol::Float64)

    # This function computes metaregression between the moderators specified in 
    # mods and the outcome variables specified in the DataFrame df. REML is used
    # in the estiamtion of the β coefficients, with the DerSimonian & Laird being
    # used as the initial estimate.

    # Get outcome variables & moderator matrix
    sch = schema(formula,df)
    ts = apply_schema(formula,sch)
    d,M = modelcols(ts,df)

    # Get DerSimonian & Laird estimate
    w = 1 ./ v
    μ = sum(w .* d) / sum(w)

    # Q-statistic
    Q = sum(w .* (d .- μ) .^ 2)

    # Get early estimate τ2
    k = length(d)
    c = sum(w) - (sum(w .^ 2) / sum(w))
    τ2 = max(0.0, (Q - (k - 1)) / c)

    cval = 0
    oldτ2 = -1.0
    while abs(τ2 - oldτ2) > tol && (cval < iter)

        # Update
        w = 1 ./ (v .+ τ2)
        W = diagm(w)

        # Get beta coefficients
        m = (M' * W * M)
        n = (M' * W * d)
        β = m \ n

        # Get predicted values
        d̄ = M * β

        # Update τ2
        f = sum(w .^ 2 .* ((d .- d̄) .^ 2 .- v)) / sum(w .^ 2) + 1 / sum(w)

        # Update values
        cval += 1
        oldτ2 = τ2
        τ2 = max(0.0, f)

    end

    # Update
    w = 1 ./ (v .+ τ2)
    μ = sum(w .* d) / sum(w)
    W = diagm(w)

    # Q-statistic
    Q = sum(w .* (d .- μ) .^ 2)

    # Get beta coefficients
    m = (M' * W * M)
    n = (M' * W * d)
    β = m \ n

    # Other outputs
    covβ = pinv(m)
    seβ = sqrt.(diag(covβ))

    # Send out the goods
    return d, v, w, μ, Q, τ2, β, covβ, seβ
end


## Organizing spline objects in meta formulas
# Define bSpline object
struct bSpline
    knots::Int
    order::Int
end

# Basic function for calling bSpline struct
function s(kn::Int, ord::Int)
    return bSpline(kn, ord)
end

# Function assigns & returns splien basis evaluated for a given moderator
function (s::bSpline)(x::AbstractVector{<:Real})

    # Get spline basis
    xmin,xmax = extrema(x)
    knotrange = LinRange(xmin,xmax,s.knots)
    basis = BSplineBasis(BSplineOrder(s.order), knotrange)

    # 3) Evaluate the basis for each x[i]
    B = zeros(Float64, length(x), length(basis))
    for (row, x_val) in enumerate(x)
        idx, vals = evaluate_all(basis, x_val)
        for (offset, b_val) in enumerate(vals)
            B[row, idx - offset + 1] = b_val
        end
    end
    return B
end

# Module-specific StatsModels.modelcols
function StatsModels.modelcols(ft::FunctionTerm{Main.meta4julia.bSpline, T}, d::NamedTuple) where T
    
    # Extract literal arguments: number of knots and order
    kn = ft.f.knots
    ord = ft.f.order
    dterm = ft.args[1].sym
    
    # Call functions
    B = s(kn,ord)(d[dterm])
    
    # Convert each column of B into a NamedTuple
    copy.(B)
end


# Generate dataframe for predicted outcomes
function predFrame(data::StepRangeLen{T},vars::String) where T

    # Initialize DataFrame
    df = DataFrame()

    # Loop through and compose dataframe
    df[!,vars] = collect(data)
    return df
end

# Generate dataframe for predicted outcomes
function predFrame(data::Vector{T},vars::Vector{String}) where T

    # Initialize DataFrame
    df = DataFrame()

    # Loop through and compose dataframe
    for (i,d) in enumerate(data)
        r = setdiff(1:length(data),i)
        for j in r
            if i == 1
                d = repeat(permutedims(d),size(data[j],1))
                d = vec(d)
            else 
                d = repeat(d,size(data[j],1))
            end
        end
        df[!,vars[i]] = d
    end
    return df
end


# Evaluate predictions for a meta-regression model
function predict(df::DataFrame, formula::FormulaTerm{T}, mdl::MRmdl; α::Float64 = 0.05) where T

    # Unpack variables
    Beta = mdl.Beta
    CovBeta = mdl.CovBeta
    tau2 = mdl.Tau2

    # Construct moderator matrix from new data
    l = size(df,1)
    df[!,formula.lhs.sym] .= ones(l,1);
    sch = schema(formula,df)
    ts = apply_schema(formula,sch)
    M = modelcols(ts.rhs,df)

    # Get predictions & standard errors
    ȳ = M * Beta
    V = diag(M * CovBeta * M')
    SE = sqrt.(V)

    # Grab z-threshold & create desired confidence intervals
    z = quantile(Normal(0, 1), 1 - α / 2)
    CI = hcat(ȳ .- z .* SE, ȳ .+ z .* SE)
    PI = hcat(ȳ .- z .* sqrt.(V .+ tau2), ȳ .+ z .* sqrt.(V .+ tau2))

    # # Back-transform y, CI, PI
    # if backTrans ~= nothing
    #     ȳ = backTrans(ȳ)
    #     CI = backTrans(CI)
    #     PI = backTrans(PI)
    # end
    # check(df,ȳ,CI,PI)

    # Return
    return ȳ, CI, PI

end

# Check size matches & variable types
function check()
end