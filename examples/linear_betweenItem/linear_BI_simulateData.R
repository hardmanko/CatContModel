###
# BetweenItem model with linear data


setwd("D:/Programming/R/CatContModel/examples/linear_betweenItem/")

source("../DataSimulatingFunctions.R")


# The range of responses is slightly wider than the range of study values.
studyRange = c(110, 250)
responseRange = c(90, 270)

# The categories don't go all the way to the edge of the study values
trueCategoryRange = c(120, 240)



# Set values for condition effect parameters
conditionEffects = list()
nCond = 3

# Probability parameter condition effects are in the latent space
conditionEffects$pMem =         c(0, -1, -0.5)  # latent
conditionEffects$pContBetween = c(0, 0.7, -0.6) # latent
conditionEffects$contSD =       c(0, 2, -2)     # manifest (weak main effect here)

# Don't use condition effects for these parameters.
conditionEffects$pCatGuess =      rep(0, nCond) 
conditionEffects$catSelectivity = rep(0, nCond)
conditionEffects$catSD =          rep(0, nCond)



# Sample values for participant parameters
set.seed(382)

nPart = 20
maxCat = 6

partParam = list()

#Participant probabilities are in the manifest space
partParam$pMem =         rtnorm(nPart, 0.8, 0.1, 0, 1) 
partParam$pContBetween = rtnorm(nPart, 0.6, 0.2, 0, 1)
partParam$pCatGuess =    rtnorm(nPart, 0.5, 0.2, 0, 1)

partParam$contSD =         rtnorm(nPart, 11, 5, 1.5)
partParam$catSelectivity = rtnorm(nPart, 10, 3, 1.5)
partParam$catSD =          rtnorm(nPart, 8, 3, 1.5)

partParam$nCat = floor(rtnorm(nPart, 3.4, 1.5, 2, 5.5))
partParam$catMu = matrix(NA, nrow=nPart, ncol=maxCat)
for (i in 1:nPart) {
	if (partParam$nCat[i] > 0) {
		partParam$catMu[i,1:partParam$nCat[i]] = makeCategories(partParam$nCat[i], maxCat, 
		                                                          noiseScale=0.2, catRange=trueCategoryRange, linear=TRUE)
	}
}

betweenSim = sampleSimulatedData(conditionEffects, partParam, trialsPerCondition=c(150, 150, 150), 
																 modelVariant = "betweenItem", dataType = "linear", 
																 responseRange = responseRange, studyRange=studyRange)

betweenData = betweenSim$simulatedData
betweenParameters = convertParameterListToDataFrame(betweenSim$combinedParam, maxCat = maxCat)


write.table(betweenData, file="linear_BI_data.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(betweenParameters, file="linear_BI_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")
