## MAIN FUNCTION
# Run Restricted maximum liklihood estiamtion on meta regression model
function reml(df::DataFrame, v::Vector{Float64}, formula::FormulaTerm, iter::Int, tol::Float64)

    # This function computes metaregression between the moderators specified in 
    # mods and the outcome variables specified in the DataFrame df. REML is used
    # in the estimation of the β coefficients, with the DerSimonian & Laird being
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
    # w = 100 .* (w ./ sum(w))
    μ = sum(w .* d) / sum(w)
    W = diagm(w)

    # Get fixed effects weights & beta for Q-statistic
    fw = 1 ./ v
    FW = diagm(fw)
    fM = (M' * FW * M)
    stXfWX = inv(fM)
    fB = stXfWX * (M' * FW) * d
    res = d .- (M * fB)
    Q = sum(fw .* (res .^ 2))

    # Estiamte I2
    k,p = size(M)
    vtot = (k - p) ./ (sum(fw) - sum(fw.^2) ./ sum(fw))
    I2 = 100 * τ2 / (τ2 + vtot)

    # Evaluate Chisq
    Qp = cdf(Chisq(k - p), Q)

    # Get beta coefficients
    m = (M' * W * M)
    n = (M' * W * d)
    β = m \ n

    # Other outputs
    covβ = pinv(m)
    seβ = sqrt.(diag(covβ))

    # Send out the goods
    return d, v, w, μ, Q, Qp, τ2, I2, β, covβ, seβ
end



## ORGANIZING SPLINE OBJECTS
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
function StatsModels.modelcols(ft::FunctionTerm{Main.jMeta.bSpline, T}, d::NamedTuple) where T
    
    # Extract literal arguments: number of knots and order
    kn = ft.f.knots
    ord = ft.f.order
    dterm = ft.args[1].sym
    
    # Call functions
    B = s(kn,ord)(d[dterm])
    
    # Convert each column of B into a NamedTuple
    copy.(B)
end



## ROBUST VARIANCE ESTIMATION
# Correlated effects robust variance estimation
function rve(mdl::model; method::String = "CR0", smallsample::Bool = true)
    # Unpack necessary 
    df = mdl.df
    w = mdl.Weights
    clusters = df[!,mdl.Clusters]
    formula = mdl.Formula
    β = mdl.Beta
    
    # Get outcome variables & moderator matrix
    sch = schema(formula,df)
    ts = apply_schema(formula,sch)
    d,M = modelcols(ts,df)

    # Compute residuals
    R = d .- (M * β)

    # Get bread for sandwich
    bread = pinv(M' * diagm(w) * M)

    # Initialize variance covariance matrix
    k,p = size(M)
    meat = zeros(Float64,p,p)

    # Switch case for determing method
    clusts = unique(clusters)
    C = length(clusts)
    if method == "CR0"
        Ac = I
    elseif method == "CR1"
        c = sqrt( C / (C - 1) )
        Ac = c * I
    elseif method == "CR2"
        sqrW = diagm(sqrt.(w))
        H = sqrW * M * bread * M' * sqrW
    end

    # Loop through clusters & construct meat
    for c in clusts

        # Subset matrices into unique elements
        idx = findall(x -> x == c, clusters)
        Mc = M[idx,:]
        Rc = R[idx]
        Wc = diagm(w[idx])

        # Switch statment for CR2 method
        if method == "CR2"

            # Estimate from Imbens & Kolesar 2016
            Phi = (Rc * Rc')
            i = CartesianIndices(Phi) 
            f = filter(x -> x.I[1] ≠ x.I[2], i)
            n = size(Phi,1)
            b = reshape(Phi'[f], n - 1, n)'
            if length(idx) > 1; ρc = sum(b) ./ length(b);
            else; ρc = 1; end
            Phi = ρc .* ones(Float64,n,n)
            Phi[diagind(Phi)] .= 1

            # Cholesky factorization of error covariance & finalize Pustejovsky & Tipton 2018
            D = cholesky(Phi).U
            Hc = H[idx,idx]
            B = D * (I - Hc) * (Rc * Rc') * (I - Hc)' * D'
            Ac = D' * sqrt(Symmetric(pinv(B))) * D
            Ac = Ac * (I - Hc)
            Mc = (I - Hc) * Mc
        end

        # Add meat to sandhwich
        meat += Mc' * Wc * Ac * (Rc * Rc') * Ac * Wc * Mc
    end

    # Small sample correction factor
    if smallsample
        f = (C / (C - 1)) * ((k - 1) / (k - p))
        meat .*= f
    end

    # Get Variance covariance matrix
    vc = bread * meat * bread
    se = sqrt.(diag(vc))

    # Return variance covariance matrix & updated standard errors
    mdl.SEBeta = se
    return mdl
end



## GENERATING PREDICTIONS
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