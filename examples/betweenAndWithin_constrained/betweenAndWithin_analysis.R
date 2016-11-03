
library(CatContModel)

setwd("~/../Programming/R/CatContModel/examples/betweenAndWithin/")

data = read.delim("betweenAndWithin_data.txt")


# For this example, I am using a constrained version of the betweenAndWithin model variant.
# I am setting pContWithin to the constant value of 0.5 for all participants.
# I am also removing any condition effects for it: Constant parameters cannot change with condition.
# See conditionEffects in the config list.

# This list will be given as the constantValueOverrides argument of runParameterEstimation().
pContWithinConst = setConstantParameterValue(data, param = "pContWithin", value = 0.5)

config = list(iterations=500, maxCategories=15, modelVariant="betweenAndWithin")

# Configure the condition effects. No pContWithin condition effects.
config$conditionEffects = list(pMem = "all", 
															 pBetween = "all", 
															 pContBetween = "all", 
															 pContWithin = "none", 
															 contSD = "all")


mhTuning = list()
mhTuning$catMu = 12
mhTuning$catSelectivity = 3
mhTuning$pCatGuess = 1
mhTuning$catSD = 1.5
mhTuning$contSD = 2.5
mhTuning$contSD_cond = 1.5
mhTuning$pContBetween = 0.6
mhTuning$pContBetween_cond = 0.3

priors = list(catMuPriorSD = 12)

results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning, 
																 priorOverrides = priors, constantValueOverrides = pContWithinConst)

examineMHAcceptance(results)


tempRes = continueSampling(results, 2500)
results = tempRes$combinedResults

#saveRDS(results, file="betweenAndWithin_results.RDS")
results = readRDS("betweenAndWithin_results.RDS")

results$post = convertPosteriorsToMatrices(results)

#convergence
activeCats = apply(results$post$catActive, 3, mean) * results$config$maxCategories
plot(activeCats, type='l')

plot(results$posteriors[["contSD_cond[2]"]], type='l')
plot(results$posteriors[["pBetween_cond[2]"]], type='l')



noBurn = removeBurnIn(results, 500)

plotParameterSummary(noBurn)

posteriorMeansAndCredibleIntervals(noBurn)

testConditionEffects(noBurn, param=c("pMem", "pContBetween", "pBetween", "contSD"))



plot(noBurn$posteriors[["pBetween_cond[3]"]], noBurn$posteriors[["contSD_cond[3]"]])


source("../DataSimulatingFunctions.R")

trueParam = read.delim("betweenAndWithin_parameters.txt")

compareTrueAndRecovered(noBurn, trueParam)
