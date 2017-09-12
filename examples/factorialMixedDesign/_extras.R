

guessFactorNames_BP(factors)

ces = getConditionEffects_BP(resList, param = "pMem", transform = TRUE)

#match(colnames(ces$prior), factors$key)

all(factors$key == ces$colKeys$key)
all(colnames(ces$prior) == ces$colKeys$key)

testSingleEffect_BP(resList, "pMem", "BP_first", factors=factors, transform = TRUE)
testSingleEffect_BP(resList, "pMem", "BP_second", factors=factors, transform = FALSE)

testSingleEffect_BP(resList, "contSD", "BP_first", factors=factors, transform = TRUE)
testSingleEffect_BP(resList, "contSD", c("BP_first", "BP_second"), factors=factors)



colnames(ces$prior) #should line up with the below
factors = data.frame(group = rep(c('A', 'B', 'C'), each=2),
										 num = rep(1:2, 3))

testHypothesis(ces$prior, ces$post, factors, "group")
testHypothesis(ces$prior, ces$post, factors, "num")
testHypothesis(ces$prior, ces$post, factors, c("group", "num"))

apply(ces$post, 2, mean)

effP = getEffectParameters(ces$post, factors, "group")
effP = getEffectParameters(ces$prior, factors, "group")

plotEffectParameterSummary(summarizeEffectParameters(effP))

groupEffectParameters(effP)

plotEffectParameterSummary(summarizeEffectParameters(ces$prior))
plotEffectParameterSummary(summarizeEffectParameters(ces$post))

param = "pContBetween"

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



mei = testMainEffectsAndInteractions_BP(resList, factors = factors, subsamples = 10)
mei = testMainEffectsAndInteractions_BP(resList, factors = factors, subsamples = 2, doPairwise = TRUE)
mei = testMainEffectsAndInteractions_BP(resList, subsamples = 2)

mei[ mei$bfType == "10", ]

mei$hasEffect = mei$bf > 3
mei[ mei$bfType == "10" & mei$hasEffect, ]
mei[ mei$bfType == "01" & mei$hasEffect, ]



cem = getConditionEffects_BP(resList, "contSD", transform = TRUE)

cmbbhtf = factors[ , guessFactorNames_BP(factors)$all ]

effp = CMBBHT::getEffectParameters(cem$post, cmbbhtf, c("BP_second"))
plotEffectParameterSummary(summarizeEffectParameters(effp))


m = CMBBHT:::makeDesignMatrix(cmbbhtf, dmFactors = c("BP_first", "BP_second"), contrastType = "contr.sum")
X = m$mat
S = solve(t(X) %*% X) %*% t(X)

zapsmall(S %*% conditionEffects$pContBetween)


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



# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("FMD_parameters.txt", stringsAsFactors = FALSE)

comp = compareTrueAndRecovered(results, trueParam)
whichRound = c("true", "rec", "cor", "slope", "dif", "percentDif")
comp[ , whichRound ] = round(comp[ , whichRound ], 2)
comp




library(CMBBHT)


###################################
# Participant level

source("~/../Programming/R/CatContModel/private/BetweenParticipants_participantLevelCombination.R")

ces = getConditionEffects_BP_participant(resList, "pMem", transform = TRUE, excludeMissing = TRUE)

factors = ces$colKeys
factors$BP_first = 1
factors$BP_first[factors$bpGroup %in% c('C', 'D')] = 2
factors$BP_second = 1
factors$BP_second[factors$bpGroup %in% c('B', 'D')] = 2
factors$Factor = factors$cond



#ff = makeFactors_BP_participant(resList)

guessFactorNames_BP_participant(factors)

tf = function(prior, post) {
	testFunction_SDDR(prior, post, min_pKept = 0)
}
tf = create_TF_EPInterval(-0.2, 0.2)

res = tf(ces$prior, ces$post)

testSingleEffect_BP_participant(resList, "pMem", "BP_first:BP_second", factors = factors)
testSingleEffect_BP_participant(resList, "pMem", "WP_cond", factors = factors)

mei = testMainEffectsAndInteractions_BP_participant(resList, factors = factors, subsamples = 1, transform=FALSE)

subFact = factors[ , c("BP_first", "BP_second", "WP_cond") ]

isDesignFullyCrossed(subFact)

fp = getEffectParameters(ces$prior, subFact, "BP_first:BP_second")
summary = summarizeEffectParameters(fp)
summary
plotEffectParameterSummary(summary)

