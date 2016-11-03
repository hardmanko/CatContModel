

setwd("~/../Programming/R/CatContModel/examples/betweenAndWithin")

source("../DataSimulatingFunctions.R")



conditionEffects = list()

conditionEffects$pMem = c(0, -0.5, -1.5) #in the latent space
conditionEffects$pContWithin = c(0, 0.5, -0.5) #in the latent space
conditionEffects$contSD = c(0, 5, 15)

# Don't use condition effects for these parameters.
nCond = 3
conditionEffects$pCatGuess = rep(0, nCond) 
conditionEffects$catSelectivity = rep(0, nCond)
conditionEffects$catSD = rep(0, nCond)


set.seed(127)

nPart = 20
maxCat = 15

partParam = list()

#Participant probabilities are in the manifest space
partParam$pMem = rtnorm(nPart, 0.8, 0.1, 0, 1) 
partParam$pContWithin = rtnorm(nPart, 0.4, 0.2, 0, 1)
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





withinSim = sampleSimulatedData(conditionEffects, partParam, trialsPerCondition=c(150, 150, 150), 
																modelVariant = "withinItem", dataType = "circular")

withinData = withinSim$simulatedData
withinParameters = convertParameterListToDataFrame(withinSim$combinedParam, maxCat = maxCat)

write.table(withinData, file="withinItem/withinItem_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

write.table(withinParameters, file="withinItem/withinItem_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")

