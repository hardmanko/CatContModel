
setwd("~/../Programming/R/CatContModel/examples/factorial_betweenItem") #or wherever you are working

library(CatContModel)

data = read.delim("factorial_BI_data.txt")


# Set up a basic configuration
config = list(iterations=500, modelVariant="betweenItem", maxCategories=15)


# Create the factors data.frame. There are two factors: "letters" and "numbers".
config$factors = data.frame(letters = c('a',  'a',  'a',  'b',  'b',  'b'),
														numbers = c('1',  '2',  '3',  '1',  '2',  '3'),
														cond =    c('a1', 'a2', 'a3', 'b1', 'b2', 'b3'),
														stringsAsFactors = FALSE)


# The names of the cond (conditions) in factors must match those in the data.
all(sort(unique(data$cond)) == sort(config$factors$cond))


# Configure the parameters with condition effects.
# Parameters that are not mentioned have no condition effects.
config$conditionEffects = list(
	pMem = "numbers",         # only use effect of numbers, not letters
	pContBetween = "letters", # only use effect of letters, not numbers
	contSD = "all"            # use all effects (letters and numbers)
)



# MH tuning steps

# 1. Override MH tuning defaults
mhTuning = list()
mhTuning$catSelectivity = 3
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.2
mhTuning$pMem_cond = 0.15


# 2. Run with those MH tuning values
results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)


# 3. Examine MH acceptance rates. Go back to step 1 until the acceptance rates look good.
examineMHAcceptance(results)


# Once the MH tuning is good, sample more iterations
continueResults = continueSampling(results, 2500)
results = continueResults$combinedResults

# And save the results
#saveRDS(results, file="factorial_BI_results.RDS")

# So that you can read them back in later
results = readRDS("factorial_BI_results.RDS")


# Examine convergence for some stuff
post = convertPosteriorsToMatrices(results, "catActive")
activeCats = apply(post$catActive, 3, mean) * results$config$maxCategories
plot(activeCats, type='l')


plot(results$posteriors[[ "pMem_cond[a2]" ]], type='l')
plot(results$posteriors[[ "pContBetween_cond[b2]" ]], type='l')


# Convergence looks pretty fast, so do only 500 burn-in iterations
results = removeBurnIn(results, 500)


# Examine the results

# Because of figure legends, you need a large plotting surface, so plot to a pdf.
pdf("parameterSummary.pdf", 10, 10)
plotParameterSummary(results)
dev.off()

# There are many, many pairs of comparisons, so it's better to look at main effects
# and interactions to determine what effects are present.
testConditionEffects(results)

# Test main effects and interactions
devfun = var

mei = testMainEffectsAndInteractions(results, summarize = FALSE, devianceFunction = devfun, subsamples = 100, subsampleProportion = 1)

meiSum = summarizeBFResults(mei)

# Examine everything at once, but this is kind of overwhelming.
meiSum

# Examine omnibus tests and bayes factors in favor of the effect.
meiSum[ meiSum$levels == "Omnibus" & meiSum$bfType == "10", ]

# Select omnibus tests with reasonably strong BFs (either for or against an effect).
meiSum[ meiSum$levels == "Omnibus" & meiSum$bf > 3, ]

# pMem has a main effect of numbers. Compare levels of the factor.
meiSum[ meiSum$param == "pMem" & meiSum$factor == "numbers", ]



#There are many condition, so only plot two of them
posteriorPredictivePlot(results, results$pnums, conditions=c("a1", "b3"), alpha=0.1)


posteriorMeansAndCredibleIntervals(results)


# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("factorial_BI_parameters.txt")

trueParam$cond = results$config$factors$cond[ trueParam$cond ]

compareTrueAndRecovered(results, trueParam)



