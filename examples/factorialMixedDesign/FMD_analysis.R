
# This example is for a factorial mixed design. In the simulated data, 
# there are 2 between-participants factors and 1 within-participants factor.
# Each factor has two levels.

setwd("~/../Programming/R/CatContModel/examples/factorialMixedDesign") #or wherever you are working

library(CatContModel)

data = read.delim("FMD_data.txt")

data$group = gsub("_\\d", "", data$cond)
data$num = gsub("\\S_", "", data$cond)

# Set up a basic configuration
config = list(iterations=500, modelVariant="betweenItem", maxCategories=12)

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
	subdata$cond = subdata$num
	
	# Remove extra columns
	subdata$group = NULL 
	subdata$num = NULL
	
	res = runParameterEstimation(config, subdata, mhTuningOverrides = mhTuning)
	
	groups[[ group ]] = res
}

# Check MH acceptance rates. You should do this for each group individually.
examineMHAcceptance(groups$A)
examineMHAcceptance(groups$B)
examineMHAcceptance(groups$C)
examineMHAcceptance(groups$D)

# Once the MH tuning seems good, sample more iterations
for (n in names(groups)) {
	continueResults = continueSampling(groups[[n]], 3500)
	groups[[n]] = continueResults$combinedResults
}

# Commented out so you don't accidentally run it and overwrite results.
#saveRDS(groups, file="FMD_resultsList_newPriors.RDS")



# TESTING REMOVE
setwd("~/../Programming/R/CatContModel/examples/factorialMixedDesign") #or wherever you are working
library(CatContModel)
# END TESTING

# Read the results back in
groups = readRDS("FMD_resultsList_newPriors.RDS")

# Once all of the groups' results objects are in a *named* list,
# you can put them into the correct format with mergeGroupResults.BP
bpRes = mergeGroupResults.BP(groups)

# Remove burn-in iterations
bpRes = removeBurnIn(bpRes, 500)
# bpRes has two elements at the top level:
# groups: A copy of the groups arguments of mergeGroupResults_BP
# config: Like results$config, but applicable to all groups (subject to some limitations)

removeBurnIn2 = function (res, burnIn) 
{
	rval = NULL
	if (resultIsType(res, "WP")) {
		rval = removeBurnIn.WP(res, burnIn)
	}
	else if (resultIsType(res, "BP")) {
		rval = removeBurnIn.BP(res, burnIn)
	}
	rval
}

wtf = function(res) {
	class(res)
}

removeBurnIn2(bpRes, 500)
wtf(bpRes)

resultIsType(res, "BP")

####################
# Tweak the factors

# mergeGroupResults_BP creates a default config$factors.
# The default assumes that there is one group factor with each
# group being its own level of the factor. Our design has two
# between-participants factors, so factors must be modified.
# Copy the factors to modify it.
factors = bpRes$config$factors

# Note that between-participants designs have new columns in factors.
# The important new ones are
# group: The name of the data set in bpRes$groups
# key: A unique identifier. key is just "group:cond"
factors


# There are four between-participants groups.
# Imagine that the between-participants part of the design is 2x2. 
# The four groups are divided into these two factors as follows.
factors$BP_first = rep(c(1,2), each=4)
factors$BP_second = rep(c(1,1,2,2), 2)

# The column BP_Factor is made by default and it treats each group as its own
# level of a factor. You don't need it because BP_first and BP_second
# are better descriptors of the design than BP_Factor.
factors$BP_Factor = NULL

# Rename the generically-named Factor to WP_cond to make it clear that it is the
# within-participants factor.
factors$WP_cond = factors$Factor
factors$Factor = NULL

# Copy the modified factors back into bpRes
bpRes$config$factors = factors

# N.B.: Don't modify the key, group, or cond columns.

##############
# Like with single-group designs, you can just add 
# a colorGeneratingFunction to the results object.
bpRes$colorGeneratingFunction = function(x) {
	hsv((x %% 360) / 360, 1, 1)
}

plotParameterSummary.BP(bpRes, factorOrder = c("WP_cond", "BP_first", "BP_second"), asPdf = TRUE)

# Plot the individual parts of the overall parameter summary
# This gives you much more control than using plotParameterSummary
plotFactorialLineChart.BP(bpRes, param = "contSD", factorOrder = c("WP_cond", "BP_first", "BP_second"))


plotCatMu.BP(bpRes)

plotHistogram(bpRes, "pMem", breaks = seq(0, 1, 0.1))



# Test main effects and interactions for the factors of the design.

# Once you have created the mixed-design factors object, you can use it when
# you call the mixed-design version of testMainEffectsAndInteractions, shown here.
# (Subsamples is set to a small value for speed of the example, but you should use more.)
mei = testMainEffectsAndInteractions.BP(bpRes, subsamples = 5)

# Examine test results as usual.
mei

# Just contSD Bayes factors in favor of the alternative
mei[ mei$param == "contSD" & mei$bfType == "10", c("factor", "bf") ]




#############
# A subsetted analysis: You only want to use one level of the WP_cond factor

# Copy original results
rlCopy = bpRes

# Select only desired rows
rlCopy$config$factors = rlCopy$config$factors[ rlCopy$config$factors$WP_cond == 1, ]

# Remove now-useless factor with only 1 level
rlCopy$config$factors$WP_cond = NULL

# Test effects
mei = testMainEffectsAndInteractions_BP(rlCopy, subsamples = 5)

mei

pdf("~/../Desktop/subset_plot.pdf", width=10, height=10)
plotParameterSummary_BP(rlCopy, factorOrder = c("BP_first", "BP_second"))
dev.off()



results = bpRes$groups$A

singleRL = mergeGroupResults_BP(list(A = bpRes$groups$A))
singleRL$config$factors$BP_Factor = NULL
testMainEffectsAndInteractions_BP(singleRL, subsamples = 5)

plotParameterSummary_BP(singleRL)

t1 = Sys.time()
plotParameterSummary(bpRes$groups$A)
Sys.time() - t1

t1 = Sys.time()
plotParameterSummary_WP(bpRes$groups$A)
Sys.time() - t1



waic = calculateWAIC_BP(bpRes)


ce1 = testConditionEffects(results)
ce2 = testConditionEffects_BP(results)

ce1$bf / ce2$bf

