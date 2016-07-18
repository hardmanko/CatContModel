
setwd("~/../Programming/R/CatContModel/examples/linear_betweenItem/")

source("../DataSimulatingFunctions.R")

###
# BetweenItem model
###


conditionEffects = list()

conditionEffects$pMem = c(0, -0.5, -1.5) #in the latent space
conditionEffects$pContBetween = c(0, 0.7, -1.2) #in the latent space
conditionEffects$contSD = c(0, 5, 15)

#The model does not have condition effects for these parameters, so they must be 0
nCond = 3
conditionEffects$pCatGuess = rep(0, nCond) 
conditionEffects$catSelectivity = rep(0, nCond)
conditionEffects$catSD = rep(0, nCond)


set.seed(127)

nPart = 20
maxCat = 7

partParam = list()

#Participant probabilities are in the manifest space
partParam$pMem = rtnorm(nPart, 0.8, 0.1, 0, 1) 
partParam$pContBetween = rtnorm(nPart, 0.6, 0.2, 0, 1)
partParam$pCatGuess = rtnorm(nPart, 0.5, 0.2, 0, 1)

partParam$contSD = rtnorm(nPart, 11, 5, 1.5)
partParam$catSelectivity = rtnorm(nPart, 13, 3, 1.5)
partParam$catSD = rtnorm(nPart, 12, 4, 1.5)

partParam$nCat = floor(rtnorm(nPart, 3.7, 1.5, 2, 5.5))
partParam$catMu = matrix(NA, nrow=nPart, ncol=maxCat)
for (i in 1:nPart) {
	if (partParam$nCat[i] > 0) {
		partParam$catMu[i,1:(partParam$nCat[i])] = makeCategories(partParam$nCat[i], maxCat, noiseScale=0.2, catRange=c(120, 240))
	}
}


betweenSim = sampleSimulatedData(conditionEffects, partParam, trialsPerCondition=c(150, 150, 150), 
																 modelVariant = "betweenItem", dataType = "linear", 
																 responseRange = c(90, 270), studyRange=c(110,250))

betweenData = betweenSim$simulatedData
betweenParameters = convertParameterListToDataFrame(betweenSim$combinedParam, maxCat = maxCat)


write.table(betweenData, file="linear_BI_data.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(betweenParameters, file="linear_BI_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")
