
# For this example, I am using a constrained version of the betweenAndWithin model variant.
# I am setting pContWithin to the constant value of 0.5 for all participants.
# This is very important for identifying the parameters because pBetween, pContBetween, 
# and pContWithin tend to trade off horribly without this constraint.
# I am also removing any condition effects for it: Constant parameters cannot change with condition.
# See conditionEffects in the config list.


library(CatContModel)

setwd("~/../Programming/R/CatContModel/examples/betweenAndWithin_constrained/")

data = read.delim("betweenAndWithin_data.txt")




# This list will be given as the constantValueOverrides argument of runParameterEstimation().
pContWithinConst = setConstantParameterValue(data, param = "pContWithin", value = 0.5)

config = list(iterations=500, maxCategories=15, modelVariant="betweenAndWithin")

# Configure the condition effects. No pContWithin condition effects.
config$conditionEffects = list(pMem = "all", 
															 pBetween = "all", 
															 pContBetween = "all", 
															 pContWithin = "none",  # <<< no effects
															 contSD = "all")


mhTuning = list()
mhTuning$catMu = 6
mhTuning$catSelectivity = 3
mhTuning$pCatGuess = 1
mhTuning$catSD = 1.5
mhTuning$contSD = 2.5
mhTuning$contSD_cond = 1.5
mhTuning$pContBetween = 0.6
mhTuning$pContBetween_cond = 0.3
mhTuning$pMem_cond = 0.12
mhTuning$pContBetween_cond = 0.25
mhTuning$pBetween = 0.7

priors = list(catMuPriorSD = 12)

results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning, 
																 priorOverrides = priors, constantValueOverrides = pContWithinConst)

examineMHAcceptance(results)


tempRes = continueSampling(results, 3000)
results = tempRes$combinedResults

#saveRDS(results, file="betweenAndWithin_results.RDS")

results = readRDS("betweenAndWithin_results.RDS")

results$post = convertPosteriorsToMatrices(results)

# Convergence
activeCats = apply(results$post$catActive, 3, mean) * results$config$maxCategories
plot(activeCats, type='l')

plot(results$posteriors[["contSD_cond[2]"]], type='l')
plot(results$posteriors[["pBetween_cond[2]"]], type='l')



results = removeBurnIn(results, 500)

plotParameterSummary(results)

posteriorMeansAndCredibleIntervals(results)

testConditionEffects(results)

mei = testMainEffectsAndInteractions(results)
mei[ mei$bfType == "10", ]


# Examine some parameter cross-correlations to determine how much parameters trade off
plot(results$posteriors[["pBetween_cond[3]"]], results$posteriors[["contSD_cond[3]"]])
plot(results$posteriors[["pBetween_cond[3]"]], results$posteriors[["pContBetween_cond[3]"]])
plot(results$posteriors[["pBetween.mu"]], results$posteriors[["pContBetween.mu"]])


source("../DataSimulatingFunctions.R")

trueParam = read.delim("betweenAndWithin_parameters.txt")

compareTrueAndRecovered(results, trueParam)
