
setwd("~/../Programming/R/CatContModel/examples/factorial_betweenItem")

source("../DataSimulatingFunctions.R")


conditionEffects = list()


#factors: a1 a2 a3 b1 b2 b3

conditionEffects$pMem = c(0, -0.5, -1.5, 0, -0.5, -1.5) #strong main effect, no interaction
conditionEffects$pContBetween = c(0, 0, 0, 0.1, 0.2, 0.3) #weak main effect, weak interaction
conditionEffects$contSD = c(0, 5, 10, 20, 15, 10) #strong interaction, plus main effect of letters

# Don't use condition effects for these parameters.
nCond = 6
conditionEffects$pCatGuess = rep(0, nCond) 
conditionEffects$catSelectivity = rep(0, nCond)
conditionEffects$catSD = rep(0, nCond)


set.seed(127)

nPart = 20
maxCat = 15

partParam = list()

#Participant probabilities are in the manifest space
partParam$pMem = rtnorm(nPart, 0.8, 0.1, 0, 1) 
partParam$pContBetween = rtnorm(nPart, 0.6, 0.2, 0, 1)
partParam$pCatGuess = rtnorm(nPart, 0.5, 0.2, 0, 1)

partParam$contSD = rtnorm(nPart, 11, 5, 1.5)
partParam$catSelectivity = rtnorm(nPart, 10, 3, 1.5)
partParam$catSD = rtnorm(nPart, 7, 3, 1.5)

partParam$nCat = floor(rtnorm(nPart, 8, 3, 0, maxCat))
partParam$catMu = matrix(NA, nrow=nPart, ncol=maxCat)
for (i in 1:nPart) {
	if (partParam$nCat[i] > 0) {
		partParam$catMu[i,1:(partParam$nCat[i])] = makeCategories(partParam$nCat[i], maxCat, noiseScale=0.2)
	}
}


betweenSim = sampleSimulatedData(conditionEffects, partParam, trialsPerCondition=rep(100, nCond), 
																 modelVariant = "betweenItem")

betweenData = betweenSim$simulatedData
betweenParameters = convertParameterListToDataFrame(betweenSim$combinedParam, maxCat = maxCat)

#Rename the conditions with factor-based names.
betweenData$cond = c('a1', 'a2', 'a3', 'b1', 'b2', 'b3')[ betweenData$cond ]
betweenParameters$cond = c('a1', 'a2', 'a3', 'b1', 'b2', 'b3')[ betweenParameters$cond ]


write.table(betweenData, file="factorial_BI_data.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(betweenParameters, file="factorial_BI_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")
