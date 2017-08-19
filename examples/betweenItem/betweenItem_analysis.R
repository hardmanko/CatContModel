
setwd("~/../Programming/R/CatContModel/examples/betweenItem") #or wherever you are working

library(CatContModel)

data = read.delim("betweenItem_data.txt")


# Set up a basic configuration
config = list(iterations=500, modelVariant="betweenItem", maxCategories=15)


# MH tuning steps

# 1. Override MH tuning defaults
mhTuning = list()
mhTuning$catSelectivity = 3
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.2
mhTuning$pMem_cond = 0.12


# 2. Run with those MH tuning values
results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)

# 3. Examine MH acceptance rates. Go back to step 1 until the acceptance rates look good.
examineMHAcceptance(results)


# Once the MH tuning is good, sample more iterations
continueResults = continueSampling(results, 3000)
results = continueResults$combinedResults

# And save the results
#saveRDS(results, file="betweenItem_results.RDS")

# So that you can read them back in later
results = readRDS("betweenItem_results.RDS")


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

gr = geweke.diag(pmat, 0.5, 0.5)
qqnorm(gr$z) #the z-scores should follow a standard normal distribution.
abline(0,1)


# Add a color generating function

results$colorGeneratingFunction = function(angle) {
	hsv((angle / 360) %% 1, 1, 1)
}


# Examine the results

plotParameterSummary(results)

testConditionEffects(results)

mei = testMainEffectsAndInteractions(results, subsamples = 20)
mei[ mei$bfType == "10", ]

posteriorMeansAndCredibleIntervals(results)

posteriorPredictivePlot(results, "1", alpha=0.3)
posteriorPredictivePlot(results, "5", alpha=0.3)


# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("betweenItem_parameters.txt")

comp = compareTrueAndRecovered(results, trueParam)
whichRound = c("true", "rec", "cor", "slope", "dif", "percentDif")
comp[ , whichRound ] = round(comp[ , whichRound ], 2)
comp


# Fit the ZL model to this data that was generated from the betweenItem model.
# We want to see that the ZL models fits the data poorly, which is done with WAIC.

zlConfig = list(iterations=3500, modelVariant="ZL", iterationsPerStatusUpdate = 200)

zlMh = list()
zlMh$contSD = 1.5
zlMh$contSD_cond = 0.7
zlMh$pMem = 0.2
zlMh$pMem_cond = 0.1

zlResults = runParameterEstimation(zlConfig, data, mhTuningOverrides = zlMh)

examineMHAcceptance(zlResults)

zlResults = removeBurnIn(zlResults, 500)

plotParameterSummary(zlResults)

# You can compare the fit of the ZL and betweenItem models with WAIC
# The data were generated under the betweenItem model, so you would expect it to win, which it does.
waic = calculateWAIC(results)
zlWaic = calculateWAIC(zlResults)

# The WAIC values of the betweenItem model should be lower than those of the ZL model
waic
zlWaic

