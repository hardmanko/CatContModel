
setwd("~/../Programming/R/CatContModel/examples/ZL")

source("../DataSimulatingFunctions.R")


nPart = 20

conditionEffects = list()
conditionEffects$pMem = c(0, -1.5, -1)
conditionEffects$contSD = c(0, 2, 10)

set.seed(124)

partParam = list()
partParam$pMem = rtnorm(nPart, 0.8, 0.1, 0.1, 0.9)

partParam$contSD = rtnorm(nPart, 13, 4, 1.5, 100)


zlSim = sampleSimulatedData(conditionEffects, partParam, trialsPerCondition=c(150, 150, 150), modelVariant = "ZL")


zlData = zlSim$simulatedData
zlParameters = convertParameterListToDataFrame(zlSim$combinedParam, maxCat = maxCat)

write.table(zlData, file="ZL_data.txt", row.names=FALSE, quote=FALSE, sep="\t")

write.table(zlParameters, file="ZL_parameters.txt", row.names=FALSE, quote=FALSE, sep="\t")


