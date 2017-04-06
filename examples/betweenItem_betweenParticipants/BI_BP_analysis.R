
setwd("~/../Programming/R/CatContModel/examples/betweenItem_betweenParticipants") #or wherever you are working

library(CatContModel)

data = read.delim("BI_BP_data.txt")


# Set up a basic configuration
config = list(iterations=500, modelVariant="betweenItem", maxCategories=15)


# MH tuning steps

# 1. Override MH tuning defaults
mhTuning = list()
mhTuning$catSelectivity = 3
mhTuning$catSD = 1.5
mhTuning$contSD = 2.5
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.2
mhTuning$pMem_cond = 0.12


# 2. Run with those MH tuning values
results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)

# 3. Examine MH acceptance rates. Go back to step 1 until the acceptance rates look good.
examineMHAcceptance(results)


# Once the MH tuning is good, sample more iterations
continueResults = continueSampling(results, 2500)
results = continueResults$combinedResults

# And save the results
#saveRDS(results, file="BI_BP_results.RDS")

# So that you can read them back in later
results = readRDS("BI_BP_results.RDS")


# Run the model for each group individually
resList = list()
for (cond in c('A', 'B', 'C')) {
	
	subdata = data[ data$cond == cond, ]
	
	res = runParameterEstimation(config, subdata, mhTuningOverrides = mhTuning)
	
	resList[[ cond ]] = res
}

examineMHAcceptance(resList$C)

for (n in names(resList)) {
	continueResults = continueSampling(resList[[n]], 2500)
	resList[[n]] = continueResults$combinedResults
}

#saveRDS(resList, file="BI_BP_individualResults.RDS")

resList = readRDS("BI_BP_individualResults.RDS")

param = "pMem"

nRes = length(resList)
iter = resList[[1]]$config$iterations #check that all have same number of iter

priorMu = postMu = matrix(nrow=iter, ncol=nRes)
for (i in 1:length(names(resList))) {
	res = resList[[ names(resList)[i] ]]
	
	postMu[,i] = res$posteriors[[ paste0(param, ".mu") ]]
	
	pmu = res$priors[[ paste0(param, ".mu.mu") ]]
	psd = sqrt(res$priors[[ paste0(param, ".mu.var") ]])
	
	priorMu[,i] = rnorm(iter, pmu, psd)
}

library(CMBBHT)

factors = data.frame(group = c('A', 'B', 'C'))

testHypothesis(priorMu, postMu, factors, "group")






# Examine convergence for some stuff
post = convertPosteriorsToMatrices(results, "catActive")
activeCats = apply(post$catActive, 3, mean) * results$config$maxCategories
plot(activeCats, type='l')


plot(results$posteriors[[ "pMem_cond[C]" ]], type='l')
plot(results$posteriors[[ "pContBetween_cond[B]" ]], type='l')

ccf(results$posteriors$pContBetween.mu, results$posteriors[[ "pContBetween_cond[B]" ]])
ccf(results$posteriors[[ "pContBetween_cond[B]" ]], results$posteriors[[ "pContBetween_cond[C]" ]])


# Convergence looks pretty fast, so do only 500 burn-in iterations
results = removeBurnIn(results, 500)


# Then double-check convergence with the Geweke diagnostic.
library(coda)

pmat = convertPosteriorsToMatrix(results)
pmat = mcmc(pmat) #convert to coda format

gr = geweke.diag(pmat, 0.5, 0.5)
qqnorm(gr$z) #the z-scores should follow a standard normal distribution.
abline(0,1)


# Examine the results

plotParameterSummary(results)

testConditionEffects(results)

mei = testMainEffectsAndInteractions(results)
mei[ mei$bfType == "10", ]

posteriorMeansAndCredibleIntervals(results)

posteriorPredictivePlot(results, results$pnums, alpha=0.1)


# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("BI_BP_parameters.txt", stringsAsFactors = FALSE)

comp = compareTrueAndRecovered(results, trueParam)
whichRound = c("true", "rec", "cor", "slope", "dif", "percentDif")
comp[ , whichRound ] = round(comp[ , whichRound ], 2)
comp


# Fit the ZL model to this data that was generated from the betweenItem model.
# We want to see that the ZL models fits the data poorly, which is done with WAIC.

zlConfig = list(iterations=3000, modelVariant="ZL", iterationsPerStatusUpdate = 200)

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
# The data were generated under the betweenItem model, so you would expect it to win, which it does.
waic = calculateWAIC(results)
zlWaic = calculateWAIC(zlResults)



