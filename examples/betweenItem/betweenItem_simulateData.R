
setwd("D:/Programming/R/CatContModel/examples/betweenItem")

source("../DataSimulatingFunctions.R")

# Set condition effects to constant values
conditionEffects = list()

# Probability parameter condition effects are in the latent space
conditionEffects$pMem =         c(0, -0.5, -1.5) # latent 
conditionEffects$pContBetween = c(0, 0.7, -1.2)  # latent
conditionEffects$contSD =       c(0, 5, 15)      # manifest

# Don't use condition effects for these parameters.
nCond = length(conditionEffects$pMem)
conditionEffects$pCatGuess =      rep(0, nCond) 
conditionEffects$catSelectivity = rep(0, nCond)
conditionEffects$catSD =          rep(0, nCond)


# Sample participant parameters
set.seed(127)

nPart = 20
maxCat = 15

partParam = list()

# Participant parameters, including probabilities, are in the manifest space
partParam$pMem =         rtnorm(nPart, 0.8, 0.1, 0, 1) 
partParam$pContBetween = rtnorm(nPart, 0.6, 0.2, 0, 1)
partParam$pCatGuess =    rtnorm(nPart, 0.5, 0.2, 0, 1)

partParam$contSD =         rtnorm(nPart, 11, 5, 1.5)
partParam$catSelectivity = rtnorm(nPart, 10, 3, 1.5)
partParam$catSD =          rtnorm(nPart, 7,  3, 1.5)

partParam$nCat = floor(rtnorm(nPart, 8, 3, 0, maxCat))
partParam$catMu = matrix(NA, nrow=nPart, ncol=maxCat)
for (i in 1:nPart) {
	if (partParam$nCat[i] > 0) {
		partParam$catMu[i,1:(partParam$nCat[i])] = makeCategories(partParam$nCat[i], maxCat, noiseScale=0.2)
	}
}

betweenSim = sampleSimulatedData(conditionEffects, partParam, 
                                 trialsPerCondition=c(150, 150, 150), 
																 modelVariant = "betweenItem")

betweenData = betweenSim$simulatedData
betweenParameters = convertParameterListToDataFrame(betweenSim$combinedParam, maxCat = maxCat)

# Save data and parameters for analysis
write.table(betweenData, file="betweenItem_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

write.table(betweenParameters, file="betweenItem_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")

