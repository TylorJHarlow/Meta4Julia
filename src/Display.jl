# Show functions for each e < model
function show(io::IO, mdl::model)
    # Basic idea: print out a table with the effect size & CI,
    # plus Q, τ², I², etc.

    println(io, "============================================================")
    println(io, "", mdl.Type)
    println(io, "============================================================")

    # Pooled effect:
    mu   = mdl.Mu
    se   = sqrt(mdl.V)               # since mdl.V is variance of mu
    ci   = mdl.CI                    # e.g. [lower, upper]
    pi   = mdl.PI                    # e.g. [lower, upper]
    qval = mdl.Q
    tau2 = mdl.Tau2
    i2   = mdl.I2

    # Print summary row
    z = (mu / se)
    pz = 2 * (1 - cdf(Normal(0,1), abs(z)))  # approximate p-value
    
    println(io, @sprintf("Estimate:                                %8.3f", mu))
    println(io, @sprintf("Std Error:                               %8.3f", se))
    println(io, @sprintf("95%% CI:                                 [%5.3f, %5.3f]", ci[1], ci[2]))
    println(io, @sprintf("95%% PI:                                 [%5.3f, %5.3f]", pi[1], pi[2]))
    println(io, @sprintf("z-value:                                 %8.3f (p = %.4f)", z, pz))

    # Q, τ², I²
    println(io, @sprintf("Q-statistic:                             %8.3f", qval))
    println(io, @sprintf("tau² (tau^2):                            %8.3f", tau2))
    println(io, @sprintf("I² (I-sq):                              %8.2f%%", i2))

    println(io, "------------------------------------------------------------")
    println(io, @sprintf("Number of studies: %d", size(mdl.df,1)))
    println(io, "============================================================\n")
end

function show(io::IO, mdl::MRmdl)
    println(io, "======================================================================")
    println(io, "", mdl.Type)
    println(io, "======================================================================")

    # Overall effect (intercept?), though in meta-regression, 'Mu' 
    # might be the intercept or a special estimate from the iteration.
    mu   = mdl.Mu   # might be Beta[1], or a separate pooled effect
    se   = sqrt(mdl.V)
    ci   = mdl.CI
    qval = mdl.Q
    tau2 = mdl.Tau2
    i2   = mdl.I2

    # Summarize the intercept if you like (optional):
    println(io, @sprintf("Overall:                                                 %8.3f", mu))
    println(io, @sprintf("Std Error:                                               %8.3f", se))
    println(io, @sprintf("95%% CI:                                                 [%5.3f, %5.3f]", ci[1], ci[2]))
    println(io, "")

    # Now handle regression coefficients
    Beta  = mdl.Beta     # vector of length p
    seBeta = mdl.SEBeta  # vector of length p
    covBeta = mdl.CovBeta

    println(io, "Regression Coefficients:")
    println(io, "---------------------------------------------------------------------")
    println(io, @sprintf("%-12s  %10s  %10s  %15s  %10s", 
        "Parameter", "Estimate", "Std.Err", "[95% CI]", "p-value"))
    println(io, "---------------------------------------------------------------------")

    # for each coefficient β_j
    z_crit = quantile(Normal(0,1), 0.975)
    for j in 1:length(Beta)
        b  = Beta[j]
        sb = seBeta[j]
        z  = b / sb
        p  = 2 * (1 - cdf(Normal(0,1), abs(z)))
        lower = b - z_crit * sb
        upper = b + z_crit * sb

        if j == 1
            param_name = @sprintf("Intercept", j - 1)
        else
            param_name = @sprintf("B[%d]", j - 1)  # or more descriptive if you know the mod names
        end
        println(io, @sprintf("%-12s  %10.3f  %10.3f    [%5.3f, %5.3f]  %10.4g", 
            param_name, b, sb, lower, upper, p))
    end 
    println(io, "---------------------------------------------------------------------")
    println(io, @sprintf("Q-statistic:    %8.3f", qval))
    println(io, @sprintf("Tau² (tau^2):   %8.3f", tau2))
    println(io, @sprintf("I² (I-sq):     %8.2f%%", i2))
    println(io, @sprintf("Number of studies: %d", size(mdl.df,1)))
    println(io, "======================================================================\n")
end

# Figures 
function forest(mdl::model,α::Float64,color::Int = 1,dir::String = "/Forest.png")

    # Unpack data frame
    df = mdl.df
    d = df.d
    v = df.v
    V = mdl.V
    CI = mdl.CI

    # Creat y vector
    y = collect(length(d):-1:1)

    # Get CI
    z = quantile(Normal(0,1), 1 - α/2)
    ci = z .* sqrt.(v) 
    ci = cat(dims = 2,d .- ci, d .+ ci)

    # Get necessary parameters
    small = minimum(CI[:])
    bigg = maximum(CI[:])
    lims = round(maximum(abs.([small,bigg])),RoundUp)
    lims = [-lims,lims]'

    # Set ideal marker size
    s = (1 ./ v)
    c = (7 / mean(s));
    s .= c .* s

    # Initiate figure
    f = scatter(d,y,color = color, markersize = s,markerstrokewidth = 0, dpi = 600,legend = false,size = (600, 600))
    for i in y
        plot!(ci[i,:],[y[i], y[i]],linecolor = color, linewidth = 2,grid = false)
        annotate!(lims[1] - (1 * (lims[2] - lims[1])),y[i],text(df.authors[i],color,:left,10))
        annotate!(lims[2] + (0.5 * (lims[2] - lims[1])),y[i],text(string(round.(d[i],digits = 2)) * " " * string(round.(ci[i,:],digits = 2)),color,:left,10))
    end
    plot!(zeros(Float64,2),[0.0, Float64(length(d))],linecolor = color, linewidth = 1,linestyle = :dash)
    scatter!([mdl.Mu], -ones(Float64,1),color = color,markersize = 15,markerstrokewidth = 0)
    plot!(CI,-ones(Float64,2),linecolor = color, linewidth = 2)
    ylims!(-3, maximum(y) + 1)
    xlim = [lims[1] - (1 * (lims[2] - lims[1])),lims[2] + (1.2 * (lims[2] - lims[1]))]
    xlims!(xlim[1],xlim[2])
    plot!(xlim,zeros(Float64,2),linecolor = color,linewidth = 2,yaxis = false)
    annotate!(lims[1] - (1 * (lims[2] - lims[1])),-1,text(mdl.Type * ":", color, :left, 10))
    annotate!(lims[2] + (0.5 * (lims[2] - lims[1])),-1,text(string(round.(mdl.Mu,digits = 2)) * " " * string(round.(CI,digits = 2)), color, :left, 10))
    annotate!(lims[1] - (1 * (lims[2] - lims[1])),-2,text("τ^2: " * string(mdl.Tau2) * ", I^2: " * string(mdl.I2), color, :left, 10))
    xticks!(lims[1]:0.5:lims[2])

    # Return plot
    savefig(f,dir);

end

function funnel(mdl::model,α::Float64 = 0.05,color::Int = 1,dir::String = "/Eggers.png")

    # Unpack data frame
    df = mdl.df
    d = df.d
    v = df.v

    # Get necessary parameters
    small = minimum(d)
    bigg = maximum(d)
    lims = round(maximum(abs.([small,bigg])),RoundUp)
    lims = [-lims,lims]'

    # Set ideal marker size
    s = (1 ./ v)
    c = (7 / mean(s));
    s .= c .* s

    # Get quantile
    z = quantile(Normal(0,1), 1 - α/2)

    # Get grid to fill
    ygrid = range(0, maximum(sqrt.(v)), length=100)
    l  = 0 .- z .* ygrid
    r = 0 .+ z .* ygrid
    x = vcat([0.0], r, reverse(l))
    y = vcat([0.0], ygrid, reverse(ygrid))

    # Run Egger's test
    reg = lm(@formula(v ~ d), df);

    # Initiate figure
    f = scatter(d,v,color = color, markersize = s, markerstrokewidth = 0, dpi = 600, legend = false)
    ylims!(0,round(maximum(v),digits = 1))
    plot!(x, y, linecolor = color)
    plot!(x,y,c = color,fillrange = zeros(Float64,length(x)), fillalpha = 0.2, linealpha = 0)
    yflip!(true)

    # Return regression results &  plot
    savefig(f,dir);
end

function eggers(mdl::model,color::Int = 1,dir::String = "/Eggers.png")
    # Unpack data frame
    df = mdl.df
    d = df.d
    v = df.v

    # Get necessary parameters
    small = minimum(d)
    bigg = maximum(d)
    lims = round(maximum(abs.([small,bigg])),RoundUp)
    lims = [-lims,lims]'

    # Set ideal marker size
    s = (1 ./ v)
    c = (7 / mean(s));
    s .= c .* s

    # New x
    x = -5:0.2:5

    # Run Egger's test
    reg = lm(@formula(d ~ v), df);

    # Initiate figure
    f = scatter(v,d,color = color, markersize = s, markerstrokewidth = 0, dpi = 600, legend = false)
    plot!(x,GLM.predict(reg,DataFrame(v = x)),linecolor = color, linewidth = 2)
    pr = GLM.predict(reg,DataFrame(v = x),interval = :prediction, level = 0.95);
    xlims!(0,round(maximum(v),digits = 1))
    ylims!(round(minimum(d),digits = 1),round(maximum(d),digits = 1))
    
    # Return regression results &  plot
    savefig(f,dir);
    return reg

end