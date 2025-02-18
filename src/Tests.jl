# Basic functions for testing model outputs

# Loads in Baskerville 2012 & tests intercept model
function testFixedEffects()

    # Select chosen meta-analysis to analyze
    study = "metadat-master/data/dat.baskerville2012.rda"
    obj = load(study)
    df = DataFrame(obj["dat.baskerville2012"])

    # Run simple equal-effects meta-analysis
    mdl = meta(df; di = :smd, vi = :se)

    # Return output for checking
    return mdl.Mu

end

# Loads in Baskerville 2012 & tests 
function testClusterRobust()

    # Select chosen meta-analysis to analyze
    study = "metadat-master/data/dat.baskerville2012.rda"
    obj = load(study)
    df = DataFrame(obj["dat.baskerville2012"])

    # Define formula 
    f = @formula(smd ~ 1)

    # Run simple equal-effects meta-analysis
    mdl = meta(df, f; cluster = :country, vi = :se)

    # Return output for checking
    return mdl.SEBeta

end