
setwd("~/../Programming/R/CatContModel/examples/factorialMixedDesign")

source("../DataSimulatingFunctions.R")


conditionEffects = list()

# factors are                     A1,A2,  B1, B2,     C1, C2,  D1, D2
conditionEffects$pMem         = c(0, 0,   0.5, -0.5,  -1, 1,   -0.5, 0.5) #in the latent space
conditionEffects$contSD       = c(0, 8,   5, 13,       0, 8,    7, 15)
conditionEffects$pContBetween = c(0, 0,   0.5, 0.6,   -1, -1,  -0.5, -0.4) #in the latent space


# Don't use condition effects for these parameters.
nCond = 8
conditionEffects$pCatGuess = rep(0, nCond) 
conditionEffects$catSelectivity = rep(0, nCond)
conditionEffects$catSD = rep(0, nCond)


set.seed(127)

nPart = 60
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


betweenSim = sampleSimulatedData(conditionEffects, partParam, trialsPerCondition=rep(200, nCond), 
																 modelVariant = "betweenItem")

betweenData = betweenSim$simulatedData
betweenParameters = convertParameterListToDataFrame(betweenSim$combinedParam, maxCat = maxCat)


#Filter out participants to make between-participant design.
condNames = paste(rep(LETTERS[1:4], each=2), 1:2, sep="_")
betweenData$cond = condNames[ betweenData$cond ]
betweenParameters$cond = condNames[ betweenParameters$cond ]

betweenParameters$cond[ betweenParameters$param %in% c("nCat", "catMu") ] = "ALL_CONDS"

groups = list(A = 1:15, B = 16:30, C = 31:45, D = 46:60)

keepData = rep(FALSE, nrow(betweenData))
keepParam = rep(FALSE, nrow(betweenParameters))

for (gr in names(groups)) {
	
	gp = groups[[gr]]
	gc = c(paste(gr, 1:2, sep="_"), "ALL_CONDS")
	
	keepData = keepData | (betweenData$pnum %in% gp & betweenData$cond %in% gc)
	keepParam = keepParam | (betweenParameters$pnum %in% gp & betweenParameters$cond %in% gc)
	
}

betweenData = betweenData[ keepData, ]
betweenParameters = betweenParameters[ keepParam, ]

#Check results
unique(betweenData[, c('pnum', 'cond')])
unique(betweenParameters[, c('pnum', 'cond')])

write.table(betweenData, file="FMD_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

write.table(betweenParameters, file="FMD_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")

