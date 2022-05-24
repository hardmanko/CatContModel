
setwd("D:/Programming/R/CatContModel/examples/halfCircle_betweenItem")

source("../DataSimulatingFunctions.R")


# Set condition effects
conditionEffects = list()

# Probability parameter condition effects are in the latent space
conditionEffects$pMem =         c(0, -0.8, -1.2) # latent
conditionEffects$pContBetween = c(0, 1, -0.6)    # latent
conditionEffects$contSD =       c(0, 12, 6)      # manifest

# Don't use condition effects for these parameters.
nCond = length(conditionEffects$pMem)
conditionEffects$pCatGuess =      rep(0, nCond) 
conditionEffects$catSelectivity = rep(0, nCond)
conditionEffects$catSD =          rep(0, nCond)


# Sample participant parameters
set.seed(957)

nPart = 20 # number of participants
maxCat = 4 # Use exactly 4 categories for each participant for this example.

partParam = list()

#Participant probabilities are in the manifest space
partParam$pMem =         rtnorm(nPart, 0.8, 0.1, 0, 1) 
partParam$pContBetween = rtnorm(nPart, 0.5, 0.2, 0, 1)
partParam$pCatGuess =    rtnorm(nPart, 0.4, 0.2, 0, 1)


# When data are sampled, they will use the full range of a circle [0,360].
# But then the data will be divided by 2 to fit in [0,180].
# Thus sample standard deviation parameters that are relatively large so that
# when the data are divided, they are more typical of observed parameter values.
partParam$contSD =         rtnorm(nPart, 19, 6, lower=1.5)
partParam$catSelectivity = rtnorm(nPart, 22, 5, lower=1.5)
partParam$catSD =          rtnorm(nPart, 16, 4, lower=1.5)



# For the halfCircle example, assume that the categories follow the standard
# pattern for orientation categories (8 total, 4 cardinal, 4 off angle).
# Then, because only half the space is being used, there are only 4 categories that show in the data.
# All participants have basically the same categories for this example, with slight perturbation.
partParam$nCat = rep(maxCat, nPart)
partParam$catMu = matrix(NA, nrow=nPart, ncol=maxCat)

# Choose angles along half the circle
partParam$catMu[,1] = 0
partParam$catMu[,2] = 45
partParam$catMu[,3] = 90
partParam$catMu[,4] = 135

# Perturb catMu a little (as much to add noise to the data as anything)
partParam$catMu = partParam$catMu + runif(length(partParam$catMu), min=-3, max=3)

# Force true category values to be in [0,180] (any that went negative)
partParam$catMu[ partParam$catMu < 0 ] = partParam$catMu[ partParam$catMu < 0 ] + 180

# Check range
range(partParam$catMu)


# The sampling code assumes that you are sampling in a full circle space.
# Multiply catMu by 2 before sampling data.
partParam$catMu = partParam$catMu * 2


# Sample data
simSample = sampleSimulatedData(conditionEffects, partParam, 
                                 trialsPerCondition=c(150, 150, 150), 
																 modelVariant = "betweenItem")


# Divide data by 2 to put it into [0,180]
simData = simSample$simulatedData
simData$study = simData$study / 2
simData$response = simData$response / 2

write.table(simData, file="halfCircle_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

# Don't modify any parameters that were used to simulate data.
# The saved parameter values for standard deviation parameters and catMu are 
# double what they should be with respect to the data.
paramDf = convertParameterListToDataFrame(simSample$combinedParam, maxCat = maxCat)

write.table(paramDf, file="halfCircle_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")

