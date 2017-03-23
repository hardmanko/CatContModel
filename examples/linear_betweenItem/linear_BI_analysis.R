
setwd("~/../Programming/R/CatContModel/examples/linear_betweenItem") #or wherever you are working

library(CatContModel)

data = read.delim("linear_BI_data.txt")


# Set up a basic configuration
config = list(iterations = 500, modelVariant = "betweenItem", dataType = "linear", maxCategories = 7)


# MH tuning steps

# 1. Override MH tuning defaults
mhTuning = list()
mhTuning$catSelectivity = 4
mhTuning$catSD = 1.6
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.2
mhTuning$pMem_cond = 0.15


# 2. Run with those MH tuning values
results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)

# 3. Examine MH acceptance rates. Go back to step 1 until the acceptance rates look good.
examineMHAcceptance(results)


# Once the MH tuning is good, sample more iterations
continueResults = continueSampling(results, 2500)
results = continueResults$combinedResults


# And save the results
#saveRDS(results, file="linear_BI_results.RDS")

# So that you can read them back in later
results = readRDS("linear_BI_results.RDS")


# Examine convergence for some stuff
post = convertPosteriorsToMatrices(results, "catActive")
activeCats = apply(post$catActive, 3, mean) * results$config$maxCategories
plot(activeCats, type='l')


plot(results$posteriors[[ "pMem_cond[2]" ]], type='l')
plot(results$posteriors[[ "pContBetween_cond[2]" ]], type='l')


# Convergence looks pretty fast, so do only 500 burn-in iterations
results = removeBurnIn(results, 500)


# Then double-check convergence with the Geweke diagnostic.
library(coda)

pmat = convertPosteriorsToMatrix(results)
pmat = mcmc(pmat) #convert to coda format

gr = geweke.diag(pmat, 0.3, 0.3)
qqnorm(gr$z) #the z-scores should follow a standard normal distribution.
abline(0,1)


# Examine the results

plotParameterSummary(results)

testConditionEffects(results)

posteriorMeansAndCredibleIntervals(results)

posteriorPredictivePlot(results, results$pnums, alpha=0.1)

testMainEffectsAndInteractions(results)


# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("linear_BI_parameters.txt")

comp = compareTrueAndRecovered(results, trueParam)
whichRound = c("true", "rec", "cor", "slope", "dif", "percentDif")
comp[ , whichRound ] = round(comp[ , whichRound ], 2)
comp



# Fit the ZL model to this data that was generated from the betweenItem model.
# We want to see that the ZL models fits the data poorly.

zlConfig = list(iterations=3000, modelVariant="ZL", dataType = "linear", iterationsPerStatusUpdate = 200)

zlMh = list()
zlMh$contSD = 1.5
zlMh$contSD_cond = 0.7
zlMh$pMem = 0.2
zlMh$pMem_cond = 0.1

zlResults = runParameterEstimation(zlConfig, data, mhTuningOverrides = zlMh)

examineMHAcceptance(zlResults)

zlResults = removeBurnIn(zlResults, 500)

plotParameterSummary(zlResults)

posteriorMeansAndCredibleIntervals(zlResults)

# You can compare the fit of the ZL and betweenItem models with WAIC
# The data were generated under the betweenItem model, so you would expect it to win.
# Lower WAIC wins.
waic = calculateWAIC(results, subsamples = 10, subsampleProportion = NULL)
zlWaic = calculateWAIC(zlResults, subsamples = 10, subsampleProportion = NULL)

waic
zlWaic
