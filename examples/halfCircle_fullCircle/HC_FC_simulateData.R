# This example is for a within-participants design that has two tasks:
# 1. A half-circle task with data in [0,180] degrees.
# 2. A full-circle task with data in [0,360] degrees.
# The data will be sampled separately for the two tasks.
# The ZL model will be used to simplify things. (No categories.)
# Further, there will be two conditions within each task.



setwd("D:/Programming/R/CatContModel/examples/halfCircle_fullCircle")

source("../DataSimulatingFunctions.R")


modelVariant = "ZL"
trialsPerCondition = c(100, 100) # Per participant

nPart = 30 # number of participants
maxCat = 0 # ZL means no categories
nCond = 2 # 2 conditions per task

set.seed(957)


# 1. Sample half-circle data.

# Set condition effects
hcCondEff = list()

# Probability parameter condition effects are in the latent space
hcCondEff$pMem =         c(0, -0.8) # latent
hcCondEff$contSD =       c(0, 15) # manifest

# Don't use condition effects for these parameters.
nCond = length(hcCondEff$pMem)
hcCondEff$pCatGuess =      rep(0, nCond) 
hcCondEff$catSelectivity = rep(0, nCond)
hcCondEff$catSD =          rep(0, nCond)
hcCondEff$pContBetween =   rep(0, nCond)


# Sample participant parameters
hcPartParam = list()

hcPartParam$pMem =   rtnorm(nPart, 0.8, 0.1, 0, 1)
hcPartParam$contSD = rtnorm(nPart, 19, 6, lower=1.5)

hcPartParam$pContBetween =   rep(0, nPart)
hcPartParam$pCatGuess =      rep(0, nPart)
hcPartParam$catSelectivity = rep(0, nPart)
hcPartParam$catSD =          rep(0, nPart)

hcPartParam$nCat = rep(maxCat, nPart)
hcPartParam$catMu = matrix(NA, nrow=nPart, ncol=maxCat)


# The sampling code assumes that you are sampling in a full circle space.

hcSimulated = sampleSimulatedData(hcCondEff, hcPartParam, 
                                 trialsPerCondition = trialsPerCondition, 
																 modelVariant = modelVariant)


# Divide data by 2 to put it into [0,180]
hcSimulated$simulatedData$study = hcSimulated$simulatedData$study / 2
hcSimulated$simulatedData$response = hcSimulated$simulatedData$response / 2

write.table(hcSimulated$simulatedData, file="halfCircle_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

# Don't modify any parameters that were used to simulate data.
# The saved parameter values for standard deviation parameters are 
# double what they should be with respect to the data.
hcParamDf = convertParameterListToDataFrame(hcSimulated$combinedParam, maxCat = maxCat)

write.table(hcParamDf, file="halfCircle_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")



# 2. Sample full circle data.

fcCondEff = hcCondEff
fcCondEff$pMem = c(0, -0.2) # latent
fcCondEff$contSD = c(0, 4)  # manifest


fcPartParam = hcPartParam

# Make participants correlated between tasks by adding a little noise to participant parameter values 
fcPartParam$pMem = hcPartParam$pMem + runif(nPart, -0.1, 0.07)
fcPartParam$pMem = pmax(0, pmin(fcPartParam$pMem, 1)) # clamp between 0 and 1

cor(hcPartParam$pMem, fcPartParam$pMem)


fcPartParam$contSD = hcPartParam$contSD + runif(nPart, -3, 5)
fcPartParam$contSD = pmax(1.5, fcPartParam$contSD) # clamp greater than 1.5

cor(hcPartParam$contSD, fcPartParam$contSD)


fcSimulated = sampleSimulatedData(fcCondEff, fcPartParam, 
                                  trialsPerCondition = trialsPerCondition, 
                                  modelVariant = modelVariant)


write.table(fcSimulated$simulatedData, file="fullCircle_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

# Don't modify any parameters that were used to simulate data.
# The saved parameter values for standard deviation parameters are 
# double what they should be with respect to the data.
fcParamDf = convertParameterListToDataFrame(fcSimulated$combinedParam, maxCat = maxCat)

write.table(fcParamDf, file="fullCircle_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")
