
# Example of using the Zhang & Luck (2008) model. This model only
# has continuous responding.
# This model is really simple compared with the categorical
# models, so it runs really really fast.

setwd("~/../Programming/R/CatContModel/examples/ZL") #or wherever you are working

library(CatContModel)

data = read.delim("ZL_data.txt")

# Set up a basic configuration
config = list(iterations=500, modelVariant="ZL")


# MH tuning steps

# 1. Override MH tuning defaults
mhTuning = list()
mhTuning$contSD = 1.5
mhTuning$contSD_cond = 0.7
mhTuning$pMem = 0.2
mhTuning$pMem_cond = 0.1

# 2. Run with those MH tuning values
results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)

# 3. Examine MH acceptance rates. Go back to step 1 until the acceptance rates look good.
examineMHAcceptance(results)


# Once the MH tuning is good, sample more iterations
results = continueSampling(results, 4500)$combinedResults

# And save the results
#saveRDS(results, file="ZL_results.RDS")

# So that you can read them back in later
results = readRDS("ZL_results.RDS")


# Examine convergence for some stuff

plot(results$posteriors[[ "pMem.mu" ]], type='l')
plot(results$posteriors[[ "contSD.mu" ]], type='l')

plot(results$posteriors[[ "pMem_cond[3]" ]], type='l')
plot(results$posteriors[[ "contSD_cond[3]" ]], type='l')


# Convergence looks pretty fast, so do only 500 burn-in iterations
results = removeBurnIn(results, 500)


# Examine the results

plotParameterSummary(results)

testMainEffectsAndInteractions(results)

testConditionEffects(results)

posteriorMeansAndCredibleIntervals(results)

posteriorPredictivePlot(results, results$pnums, alpha=0.1)


# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("ZL_parameters.txt")
trueParam = trueParam[ trueParam$param != "nCat", ]

compareTrueAndRecovered(results, trueParam)


