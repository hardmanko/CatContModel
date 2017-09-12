
# This example is for a case when you have linear (non-circular) data.
# Please see the Linear Data section in the manual. 
#
# In this example, the range of possible study angles is [110, 250] and the 
# range of possible response angles is [90, 270]. It is usually a good idea
# to allow participants to respond outside of the range of studied values.
# Hard edges at the ends of the range of study values can result in a large
# mass of responses right at the edge of the range.
#
# Although the values 90, 270, etc might look like degrees in a half circle, 
# the model works by treating them as being linear (no wrap-around at 360).

setwd("~/../Programming/R/CatContModel/examples/linear_betweenItem") #or wherever you are working

library(CatContModel)

data = read.delim("linear_BI_data.txt")

# Set up a basic configuration. 
# Note dataType = "linear" (defaults to "circular").
config = list(iterations = 500, modelVariant = "betweenItem", 
							dataType = "linear", maxCategories = 7)

# Optional: Specify the response range. If not specified, it is taken to be range(data$response).
config$responseRange = c(90, 270)


# MH tuning steps

# 1. Override MH tuning defaults
mhTuning = list()
mhTuning$catMu = 3
mhTuning$catSelectivity = 4
mhTuning$catSD = 2
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.5
mhTuning$pContBetween_cond = 0.2
mhTuning$pMem_cond = 0.15


# 2. Run with those MH tuning values
results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)

# 3. Examine MH acceptance rates. Go back to step 1 until the acceptance rates look good.
examineMHAcceptance(results)


# Once the MH tuning is good, sample more iterations
continueResults = continueSampling(results, 3000)
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

mei = testMainEffectsAndInteractions(results)
mei[ mei$bfType == "10", ]

testConditionEffects(results)

posteriorMeansAndCredibleIntervals(results)

posteriorPredictivePlot(results, "8", alpha=0.4)
posteriorPredictivePlot(results, results$pnums, alpha=0.1)



# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("linear_BI_parameters.txt")

comp = compareTrueAndRecovered(results, trueParam)
whichRound = c("true", "rec", "cor", "slope", "dif", "percentDif")
comp[ , whichRound ] = round(comp[ , whichRound ], 2)
comp



# Fit the ZL model to this data that was generated from the betweenItem model.
# We want to see that the ZL models fits the data poorly.

zlConfig = list(iterations = 3000, modelVariant = "ZL", dataType = "linear", iterationsPerStatusUpdate = 200)

zlMh = list()
zlMh$contSD = 2
zlMh$contSD_cond = 1
zlMh$pMem = 0.3
zlMh$pMem_cond = 0.12

zlResults = runParameterEstimation(zlConfig, data, mhTuningOverrides = zlMh)

examineMHAcceptance(zlResults)

zlResults = removeBurnIn(zlResults, 500)

plotParameterSummary(zlResults)

# You can compare the fit of the ZL and betweenItem models with WAIC
# The data were generated under the betweenItem model, so you would expect it to win.
# Lower WAIC wins.
waic = calculateWAIC(results, subsamples = 10, subsampleProportion = NULL)
zlWaic = calculateWAIC(zlResults, subsamples = 10, subsampleProportion = NULL)

waic[ waic$stat == "WAIC_2", ]
zlWaic[ zlWaic$stat == "WAIC_2", ]
