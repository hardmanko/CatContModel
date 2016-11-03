
setwd("~/../Programming/R/CatContModel/examples/betweenAndWithin")

source("../DataSimulatingFunctions.R")

###
# Between and within, fixed pContWithin = 0.5
###

conditionEffects = list()

conditionEffects$pMem = c(0, -0.5, -1) #in the latent space
conditionEffects$pBetween = c(0, 1, 0)
conditionEffects$pContBetween = c(0, 1, -0.5) #in the latent space
conditionEffects$pContWithin = c(0, 0, 0)
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
partParam$pBetween = rtnorm(nPart, 0.5, 0.1, 0, 1) 
partParam$pContWithin = rep(0.5, nPart)
partParam$pContBetween = rtnorm(nPart, 0.5, 0.2, 0, 1)
partParam$pCatGuess = rtnorm(nPart, 0.5, 0.2, 0, 1)

partParam$contSD = rtnorm(nPart, 11, 5, 1.5)
partParam$catSelectivity = rtnorm(nPart, 10, 3, 1.5)
partParam$catSD = rtnorm(nPart, 8, 3, 1.5)

partParam$nCat = floor(rtnorm(nPart, 8, 3, 0, maxCat))
partParam$catMu = matrix(NA, nrow=nPart, ncol=maxCat)
for (i in 1:nPart) {
	if (partParam$nCat[i] > 0) {
		partParam$catMu[i,1:(partParam$nCat[i])] = makeCategories(partParam$nCat[i], maxCat, noiseScale=0.2)
	}
}


betweenAndWithin = sampleSimulatedData(conditionEffects, partParam, trialsPerCondition=c(150, 150, 150), 
																			 modelVariant = "betweenAndWithin", dataType = "circular")

bwData = betweenAndWithin$simulatedData
bwParameters = convertParameterListToDataFrame(betweenAndWithin$combinedParam, maxCat = maxCat)


write.table(bwData, file="betweenAndWithin_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

write.table(bwParameters, file="betweenAndWithin_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")
