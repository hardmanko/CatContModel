
setwd("~/../Programming/R/CatContModel/examples/withinItem") #or wherever you are working

library(CatContModel)

data = read.delim("withinItem_data.txt")

# Set up a basic configuration
config = list(iterations=500, modelVariant="withinItem", maxCategories=15)


# MH tuning steps

# 1. Override MH tuning defaults
mhTuning = list()
mhTuning$catSelectivity = 3
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.5
mhTuning$pContWithin = 0.2
mhTuning$pContWithin_cond = 0.05
mhTuning$pMem_cond = 0.10
mhTuning$catMu = 4

# 2. Run with those MH tuning values
results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)

# 3. Examine MH acceptance rates. Go back to step 1 until the acceptance rates look good.
examineMHAcceptance(results)


# Once the MH tuning is good, sample more iterations
tempRes = continueSampling(results, 2500)
results = tempRes$combinedResults

# And save the results
#saveRDS(results, file="withinItem_results.RDS")

# So that you can read them back in later
results = readRDS("withinItem_results.RDS")



# Examine convergence for some stuff
post = convertPosteriorsToMatrices(results, "catActive")
activeCats = apply(post$catActive, 3, mean) * results$config$maxCategories
plot(activeCats, type='l')

plot(results$posteriors$`contSD_cond[3]`, type='l')


# Convergence looks pretty fast, so do only 500 burn-in iterations
results = removeBurnIn(results, 500)


# Examine the results

plotParameterSummary(results)

testConditionEffects(results)

posteriorMeansAndCredibleIntervals(results)

posteriorPredictivePlot(results, results$pnums, alpha=0.1)


# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("withinItem_parameters.txt")

compareTrueAndRecovered(results, trueParam)


