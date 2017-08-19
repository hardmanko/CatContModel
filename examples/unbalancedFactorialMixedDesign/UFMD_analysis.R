
# This example is of an unbalanced factorial mixed design.
# It shows how to work with
# 1. A mixed design, with both within-participant and between-participant factors
# 2. An unbalanced design, where neither the within-participant nor the between-participant factors are fully crossed.


setwd("~/../Programming/R/CatContModel/examples/unbalancedFactorialMixedDesign") #or wherever you are working

library(CatContModel)

data = read.delim("UFMD_data.txt")

data$group = gsub("_\\d", "", data$cond)
data$num = gsub("\\S_", "", data$cond)

# Modify the basic data to make the design unbalanced

data = data[ data$cond != "B_2", ]

data$cond = as.character(data$cond)

data$num[ data$cond == "D_1" ] = 3
data$num[ data$cond == "D_2" ] = 4

data$cond[ data$cond == "D_1" ] = "C_3"
data$cond[ data$cond == "D_2" ] = "C_4"

data$group[ data$group == "D" ] = "C"

# Check out the resulting design
unique(data[ , c("cond", "group", "num") ])




# Set up the same basic configuration for each group
configs = list()
for (group in unique(data$group)) {
	configs[[ group ]] = list(iterations=300, modelVariant="betweenItem", maxCategories=12)
}

# Modify the basic config by specifying the factors of the design.
# This is important to do, even for group B, where the factors are trivial.

configs$A$factors = data.frame(cond = c("A_1", "A_2"),
															 letters = c("a", "b"),
															 numbers = c(1, 1),
															 stringsAsFactors = FALSE)

configs$B$factors = data.frame(cond = "B_1",
															 letters = "a",
															 numbers = 1,
															 stringsAsFactors = FALSE)

configs$C$factors = data.frame(cond = paste0("C_", 1:4),
															 letters = c("a", "a", "b", "b"),
															 numbers = c(1, 2, 1, 2),
															 stringsAsFactors = FALSE)

# For group C, specify which factors may vary
configs$C$conditionEffects = list(pMem = "letters",
																	contSD = "numbers",
																	pContBetween = "all")


# MH tuning
mhTuning = list()
mhTuning$catMu = 4
mhTuning$catSelectivity = 2
mhTuning$catSD = 1.5
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.15
mhTuning$pMem_cond = 0.12


# Run the model for each group individually
# Depending on your data, you may want or need to use different MH tuning for each group
groups = list()
for (group in unique(data$group)) {
	
	subdata = data[ data$group == group, ]
	
	# Remove extra columns
	subdata$group = NULL 
	subdata$num = NULL
	
	res = runParameterEstimation(configs[[ group ]], subdata, mhTuningOverrides = mhTuning)
	
	groups[[ group ]] = res
}

# Check MH acceptance rates. You should do this for each group individually.
examineMHAcceptance(groups$A)
examineMHAcceptance(groups$B)
examineMHAcceptance(groups$C)

# Once the MH tuning seems good, sample more iterations
for (n in names(groups)) {
	continueResults = continueSampling(groups[[n]], 2700)
	groups[[n]] = continueResults$combinedResults
}

# Commented out so you don't accidentally run it and overwrite results.
#saveRDS(groups, file="UFMD_groups.RDS")



############################
# TESTING: REMOVE
setwd("~/../Programming/R/CatContModel/examples/unbalancedFactorialMixedDesign")
library(CatContModel)
############################

# Read the results back in
groups = readRDS("UFMD_groups.RDS")

#groups$B = NULL

bpRes = mergeGroupResults.BP(groups)

bpRes = removeBurnIn(bpRes, 10)

bpRes$colorGeneratingFunction = function(angle) {
	hsv((angle / 360) %% 1, 1, 1)
}

#####################################


testSingleEffect(bpRes, "pContBetween", "numbers:letters")

CatContModel:::testMEI_singleParameter.Generic(bpRes, "pMem")
CatContModel:::testMEI_singleParameter.Generic(bpRes$groups$C, "pMem")

mei = CatContModel:::testMainEffectsAndInteractions.Generic(bpRes, subsamples = 5)
mei = CatContModel:::testMainEffectsAndInteractions.Generic(bpRes$groups$C, subsamples = 5)

mei[ mei$bfType == "10", ]


ceht = testConditionEffects.Generic(bpRes)
ceht = testConditionEffects.Generic(bpRes$groups$C)

plotFactorialLineChart(bpRes$groups$C, "contSD")


cef.bp = getConditionEffects(bpRes, "pMem")
cef.wp = getConditionEffects(bpRes$groups$C, "pMem")

param = "contSD"
plotDf = plotFactorialLineChart(bpRes, param)
plotDf = plotFactorialLineChart(bpRes, param, factorOrder = c("letters", "numbers"))
plotDf = plotFactorialLineChart(bpRes, param, factorOrder = c("letters", "BP_Factor"))


plotFactorialLineChart(bpRes$groups$C, "pMem")


plotParameterSummary(bpRes, factorOrder = c("letters", "BP_Factor", "numbers"), asPdf=TRUE, pdfScale = 2)
plotParameterSummary(bpRes, asPdf=TRUE, pdfScale = 2)

plotParameterSummary(bpRes$groups$C, asPdf=TRUE, pdfScale = 2)


plotCatMu(bpRes)

bpRes$groups$C$colorGeneratingFunction = function(angle) {
	hsv((angle / 360) %% 1, 1, 1)
}

plotCatMu(bpRes$groups$C)

plotFactorialLineChart(bpRes, "pMem", factorOrder = c("letters", "BP_Factor"))


plotFactorialLineChart(bpRes, "catSD")

plotHistogram(bpRes, "catSD")

pp = getAllParameterPosteriors(bpRes, "catSD", manifest = TRUE, format = "data.frame")


getFactorNameToType(fact)


ceg = CatContModel:::getConditionEffects.BP(bpRes, "contSD", addMu = TRUE, manifest = TRUE)

hist(ceg$prior)


plotParameterSummary.BP(bpRes$groups$C, factorOrder = c("letters", "numbers"), asPdf=TRUE, pdfScale = 2)






f = bpRes$groups$C$config$factors
f
normalizeFactors(f)

normalizeFactors(bpRes$config$factors)

#config = bpRes$config
config = bpRes$groups$C$config
param = "pMem"

f = config$factors
#f = f[ c(1, 2, 3, 5, 7, 8), ]
#f = f[ c(1, 3, 5, 7), ]
f = f[ 1:3, ]
config$factors = f

config$conditionEffects$pMem = "none"

getFactorsForConditionEffect.WP(config, param)
getFactorsForConditionEffect.BP(bpRes, param)




updateFactorsForConditionEffects(bpRes, "pMem")
updateFactorsForConditionEffects(bpRes, "contSD")
updateFactorsForConditionEffects(bpRes, "pContBetween")

updateFactorsForConditionEffects(bpRes$groups$C, "pMem")


getvaryingFactorNames(bpRes, "pMem")
getFactorsForConditionEffect(bpRes, "pMem")


bpRes$config$factors
getFactorsForConditionEffect(bpRes, "pMem")
getFactorsForConditionEffect.BP(bpRes, "contSD")
getFactorsForConditionEffect.BP(bpRes, "pContBetween")


getFactorTypeToName(bpRes$config$factors)

bpCopy = bpRes
fact = bpCopy$config$factors
fact$letters = "a"
bpCopy$config$factors = fact
bpCopy = CatContModel:::backPropogateBPFactorsToWPFactors(bpCopy)

bpCopy$groups$A$config$factors
bpCopy$groups$B$config$factors
bpCopy$groups$C$config$factors



bpCopy = bpRes
bpCopy
makeBPConditionEffects(bpCopy)





subs = CatContModel:::calculateWAIC_subsamples(bpRes, 1, 0.1)




