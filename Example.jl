# Load in necessary modules
using .jMeta, RData, DataFrames, StatsModels

# Select chosen meta-analysis to analyze
study = "/Applications/Toolbox/metadat-master/data/dat.baskerville2012.rda"
obj = load(study)
df = DataFrame(obj["dat.baskerville2012"])

# Change Standard errors to variance
df[!,:v] .= df[!,:se].^2

# Define formula 
f = @formula(smd ~ 1)

# Run simple equal-effects meta-analysis
mdl = meta(df, f; v = :v)
show(mdl)