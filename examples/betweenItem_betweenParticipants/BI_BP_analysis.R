# This example is for a one-factor between-participants (BP) design.
# It shows how to 


setwd("~/../Programming/R/CatContModel/examples/betweenItem_betweenParticipants") #or wherever you are working

library(CatContModel)

# There is data from three groups in this data set,
# with group membership given by the cond variable.
data = read.delim("BI_BP_data.txt")

# You can see that pnum X cond does not overlap.
unique( data[, c("pnum", "cond") ])


# Set up a basic configuration
config = list(iterations=500, modelVariant="betweenItem", maxCategories=15)

# If you were doing this for real, you might want to use different MH tunings
# for the different groups.
mhTuning = list()
mhTuning$catSelectivity = 3
mhTuning$catSD = 1.5
mhTuning$contSD = 2.5
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.2
mhTuning$pMem_cond = 0.12


# Run the model for each group individually
resList = list()
for (cond in c('A', 'B', 'C')) {
	
	# Select out only the data for this one group.
	subdata = data[ data$cond == cond, ]
	
	res = runParameterEstimation(config, subdata, mhTuningOverrides = mhTuning)
	
	resList[[ cond ]] = res
}

# Check MH acceptance for the different groups
examineMHAcceptance(resList$A)

# Sample more
for (n in names(resList)) {
	continueResults = continueSampling(resList[[n]], 2500)
	resList[[n]] = continueResults$combinedResults
}

#saveRDS(resList, file="BI_BP_individualResults.RDS")

resList = readRDS("BI_BP_individualResults.RDS")

for (n in names(resList)) {
	resList[[n]] = removeBurnIn(resList[[n]], 500)
}

############
# Using the CMBBHT package, perform tests of the between-participants main effect.

library(CMBBHT)

# Pick a parameter. It can be any parameter in the model other than catMu and catActive.

param = "pContBetween" 

# You should be sure that all results have same number of iterations
iter = resList[[1]]$config$iterations 

# Get matrices of prior and posterior samples of the grand mean parameters (e.g. pMem.mu)
priorMu = postMu = matrix(nrow=iter, ncol=length(resList))
for (i in 1:length(names(resList))) {
	res = resList[[ names(resList)[i] ]]
	
	postMu[,i] = res$posteriors[[ paste0(param, ".mu") ]]
	
	# Get prior parameter values
	pmu = res$priors[[ paste0(param, ".mu.mu") ]]
	psd = sqrt(res$priors[[ paste0(param, ".mu.var") ]])
	
	# And sample from the prior
	priorMu[,i] = rnorm(iter, pmu, psd)
}


# Create the factors data frame for CMBBHT
factors = data.frame(group = c('A', 'B', 'C'))

# And test the main effect of group
testHypothesis(priorMu, postMu, factors, "group")



# This package does not have very good support for between-participants designs,
# but I'm working on adding some features. For now, most analyses are awkward.



