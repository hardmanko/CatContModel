
# This example is of an unbalanced factorial mixed design.
# It shows how to work with
# 1. A mixed design, with both within-participant and between-participant factors
# 2. An unbalanced design, where neither the within-participant nor the between-participant factors are fully crossed.


setwd("D:/Programming/R/CatContModel/examples/unbalancedFactorialMixedDesign") #or wherever you are working

library(CatContModel)

data = read.delim("UFMD_data.txt")

# What are the conditions?
data$cond = as.character(data$cond)
unique(data$cond)

# Read out the group letter and number from cond
data$group = gsub("_\\d", "", data$cond)
data$num = gsub("\\S_", "", data$cond)


# Modify the basic data to make the design unbalanced in interesting ways

# Drop B_2
data = data[ data$cond != "B_2", ] 

# Rename group D to group C and change D_1 and D_2 to C_3 and C_4
data$group[ data$group == "D" ] = "C"
data$num[ data$cond == "D_1" ] = 3
data$num[ data$cond == "D_2" ] = 4

# Redo cond from group and num
data$cond = paste(data$group, data$num, sep="_") 

# Check out the resulting design
unique(data[ , c("cond", "group", "num") ])



# Set up the same basic configuration for each group
configs = list()
for (group in unique(data$group)) {
	configs[[ group ]] = list(iterations=300, modelVariant="betweenItem", maxCategories=12)
}

# Modify the basic config by specifying the factors of the design.
# This is important to do, even for group B, where the factors are trivial.
#
# There are two within-participants factors: letters and numbers.
# There is one between-participant factor: BP_Group. (Not "group", which is a reserved factor name.)

configs$A$factors = data.frame(cond = c("A_1", "A_2"),
															 letters = c("a", "b"),
															 numbers = c(1, 1),
															 BP_Group = "A",
															 stringsAsFactors = FALSE)

configs$B$factors = data.frame(cond = "B_1",
															 letters = "a",
															 numbers = 1,
															 BP_Group = "B",
															 stringsAsFactors = FALSE)

configs$C$factors = data.frame(cond = paste0("C_", 1:4),
															 letters = c("a", "a", "b", "b"),
															 numbers = c(1, 2, 1, 2),
															 BP_Group = "C",
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
	continueResults = continueSampling(groups[[n]], 3000)
	groups[[n]] = continueResults$combinedResults
}

# Commented out so you don't accidentally run it and overwrite results.
#saveRDS(groups, file="UFMD_groups.RDS")


##########################
# Read the results back in
groups = readRDS("UFMD_groups.RDS")

bpRes = combineGroupResults.BP(groups)

bpRes = removeBurnIn(bpRes, 500)

bpRes$colorGeneratingFunction = function(angle) {
	hsv((angle / 360) %% 1, 1, 1)
}

##############
# Plotting and analysis goes more or less as usual with unbalanced designs

plotParameterSummary(bpRes, asPdf=TRUE, pdfScale = 2)

plotParameter(bpRes, "pMem")
plotParameter(bpRes, "catSD")
plotParameter(bpRes, "catActive")
plotParameter(bpRes, "catMu")

plotParameterLineChart(bpRes, "catSD")
plotParameterHistogram(bpRes, "catSD")


##########
# Main effects and interactions
mei = testMainEffectsAndInteractions(bpRes, subsamples = 5)

mei[ mei$bfType == "10", ]


#########
# Due to how ugly this design is, you might want to compare cells of the design directly.
ce = testConditionEffects(bpRes)
ce

# Just bayes factors in favor of sameness (no difference)
ce[ ce$param == "pMem" & ce$bfType == "01", ]

# To help you understand what cells are being compared, use the key column
# of factors and compare to the key column in "ce".
bpRes$config$factors

#########
# Other stuff
posteriorMeansAndCredibleIntervals(bpRes)

waic = calculateWAIC(bpRes)



############################################################
# Imagine that you want to do an analysis with only group C
# You can use the results object stored in bpRes$groups$C
grpC = bpRes$groups$C

testMainEffectsAndInteractions(grpC)


