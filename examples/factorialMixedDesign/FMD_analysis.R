
# This example is for a factorial mixed design. In the simulated data, 
# there are 2 between-participants factors and 1 within-participants factor.
# Each factor has two levels.

setwd("~/../Programming/R/CatContModel/examples/factorialMixedDesign") #or wherever you are working

library(CatContModel)

data = read.delim("FMD_data.txt")

# Get group information from cond names
data$group = gsub("_\\d", "", data$cond)

# Examine the design, which is 4 groups with a within-participants factor with 2 levels
unique(data[ , c("group", "cond") ])

# Set up the same basic configuration for each group
configs = list()
for (group in unique(data$group)) {
	configs[[ group ]] = list(iterations=500, modelVariant="betweenItem", maxCategories=12)
}

# Modify the basic config by specifying the factors of the design.
# You must do this before running parameter estimation.
#
# The numbers factor is within-participants (varying within groups).
# The letters and LETTERS factors are between-participants (varying between groups).

configs$A$factors = data.frame(cond    = c("A_1", "A_2"),
															 numbers = c(1, 2),
															 letters = c('a', 'a'),
															 LETTERS = c('A', 'A'),
															 stringsAsFactors = FALSE)

configs$B$factors = data.frame(cond    = c("B_1", "B_2"),
															 numbers = c(1, 2),
															 letters = c('a', 'a'),
															 LETTERS = c('B', 'B'),
															 stringsAsFactors = FALSE)

configs$C$factors = data.frame(cond    = c("C_1", "C_2"),
															 numbers = c(1, 2),
															 letters = c('b', 'b'),
															 LETTERS = c('A', 'A'),
															 stringsAsFactors = FALSE)

configs$D$factors = data.frame(cond    = c("D_1", "D_2"),
															 numbers = c(1, 2),
															 letters = c('b', 'b'),
															 LETTERS = c('B', 'B'),
															 stringsAsFactors = FALSE)



# MH tuning configuration
# Although it is the same for each group in this example, with real data
# you may need to use different MH tuning values for different groups.
mhTuning = list()
mhTuning$catMu = 4
mhTuning$catSelectivity = 2
mhTuning$catSD = 1.5
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.15
mhTuning$pMem_cond = 0.12


# Run the model for each group individually and store the results for each group in "groups"
groups = list()
for (grp in unique(data$group)) {
	
	# Select only the data for the current group
	subdata = data[ data$group == grp, ]
	
	# Remove extra group column
	subdata$group = NULL 
	
	res = runParameterEstimation(configs[[ grp ]], subdata, mhTuningOverrides = mhTuning)
	
	groups[[ grp ]] = res
}

# Check MH acceptance rates. You should do this for each group individually.
examineMHAcceptance(groups$A)
examineMHAcceptance(groups$B)
examineMHAcceptance(groups$C)
examineMHAcceptance(groups$D)

# Once the MH tuning seems good, sample more iterations
for (n in names(groups)) {
	continueResults = continueSampling(groups[[n]], 3000)
	groups[[n]] = continueResults$combinedResults
}

# Commented out so you don't accidentally run it and overwrite results.
#saveRDS(groups, file="FMD_groups.RDS")



# Read the results back in
groups = readRDS("FMD_groups.RDS")

# Once all of the groups' results objects are in a *named* list,
# you can put them into the correct format with combineGroupResults.BP
bpRes = combineGroupResults.BP(groups)

# Remove burn-in iterations
bpRes = removeBurnIn(bpRes, 500)
# bpRes has two elements at the top level:
# groups: A copy of the groups arguments of mergeGroupResults.BP
# config: Like results$config, but applicable to all groups (subject to some limitations)



##############
# Like with single-group designs, you can just add 
# a colorGeneratingFunction to the results object.
bpRes$colorGeneratingFunction = function(x) {
	hsv((x %% 360) / 360, 1, 1)
}

########################
# Visualizing parameters

# Overall plot
plotParameterSummary(bpRes, asPdf = TRUE)

# Plot the individual parameters of the overall parameter summary.
# This gives you much more control than using plotParameterSummary.

plotParameterLineChart(bpRes, "pContBetween")

plotParameterHistogram(bpRes, "pMem", breaks = seq(0, 1, 0.1))

plotCatMu(bpRes)

# Plot any parameter with default settings
plotParameter(bpRes, "catActive")

#############
# Test main effects and interactions for the factors of the design.
# Due to the very large number of tests, a small number of subsamples is used
# for the purposes of speeding up this example.
mei = testMainEffectsAndInteractions(bpRes, subsamples = 5)

# Examine test results as usual.
mei

# Just contSD Bayes factors in favor of the alternative
mei[ mei$param == "contSD" & mei$bfType == "10", c("factor", "bf") ]

# And plot just that parameter to check the reasonableness of the results
plotParameterLineChart(bpRes, "contSD", factorOrder = c("LETTERS", "letters", "numbers"))


#############
# A subsetted analysis: You only want to use one level of the numbers factor

# Copy original results
resCopy = bpRes

# Copy over the factors
factors = resCopy$config$factors

# Select only desired rows
factors = factors[ factors$numbers == 1, ]

# Remove now-useless factor with only 1 level
factors$numbers = NULL

# Copy back over
resCopy$config$factors = factors

# Test effects
mei = testMainEffectsAndInteractions(resCopy, subsamples = 5)

mei

# Also make plots for subsetted design
plotParameterSummary(resCopy, asPdf = TRUE)



