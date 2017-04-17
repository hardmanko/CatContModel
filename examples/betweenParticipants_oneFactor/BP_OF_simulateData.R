
setwd("~/../Programming/R/CatContModel/examples/betweenParticipants_oneFactor")

source("../DataSimulatingFunctions.R")


conditionEffects = list()

conditionEffects$pMem = c(0, -0.5, -1.5) #in the latent space
conditionEffects$pContBetween = c(0, 0.7, -1.2) #in the latent space
conditionEffects$contSD = c(0, 5, 15)

# Don't use condition effects for these parameters.
nCond = 3
conditionEffects$pCatGuess = rep(0, nCond) 
conditionEffects$catSelectivity = rep(0, nCond)
conditionEffects$catSD = rep(0, nCond)


set.seed(127)

nPart = 45
maxCat = 12

partParam = list()

#Participant probabilities are in the manifest space
partParam$pMem = rtnorm(nPart, 0.8, 0.1, 0, 1) 
partParam$pContBetween = rtnorm(nPart, 0.6, 0.2, 0, 1)
partParam$pCatGuess = rtnorm(nPart, 0.5, 0.2, 0, 1)

partParam$contSD = rtnorm(nPart, 11, 5, 1.5)
partParam$catSelectivity = rtnorm(nPart, 10, 3, 1.5)
partParam$catSD = rtnorm(nPart, 7, 3, 1.5)

partParam$nCat = floor(rtnorm(nPart, 7, 3, 2, maxCat))
partParam$catMu = matrix(NA, nrow=nPart, ncol=maxCat)
for (i in 1:nPart) {
	if (partParam$nCat[i] > 0) {
		partParam$catMu[i,1:(partParam$nCat[i])] = makeCategories(partParam$nCat[i], maxCat, noiseScale=0.2)
	}
}


betweenSim = sampleSimulatedData(conditionEffects, partParam, trialsPerCondition=c(300, 300, 300), 
																 modelVariant = "betweenItem")

betweenData = betweenSim$simulatedData
betweenParameters = convertParameterListToDataFrame(betweenSim$combinedParam, maxCat = maxCat)


#Filter out participants to make between-participant design.
betweenData$cond = LETTERS[betweenData$cond]
betweenParameters$cond = LETTERS[betweenParameters$cond]
groups = list(A = 1:15, B = 16:30, C = 31:45)

betweenParameters$cond[ betweenParameters$param %in% c("nCat", "catMu") ] = "ALL_CONDS"

keepData = rep(FALSE, nrow(betweenData))
keepParam = rep(FALSE, nrow(betweenParameters))

for (gr in LETTERS[1:3]) {
	
	gp = groups[[gr]]
	gr = c(gr, "ALL_CONDS")
	
	keepData = keepData | (betweenData$pnum %in% gp & betweenData$cond %in% gr)
	keepParam = keepParam | (betweenParameters$pnum %in% gp & betweenParameters$cond %in% gr)
	
}

betweenData = betweenData[ keepData, ]
betweenParameters = betweenParameters[ keepParam, ]

#Check results
unique(betweenData[, c('pnum', 'cond')])
unique(betweenParameters[, c('pnum', 'cond')])

write.table(betweenData, file="BP_OF_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

write.table(betweenParameters, file="BP_OF_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")

