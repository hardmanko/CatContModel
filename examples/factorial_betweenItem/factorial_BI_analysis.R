
setwd("D:/Programming/R/CatContModel/examples/factorial_betweenItem") #or wherever you are working

library(CatContModel)

data = read.delim("factorial_BI_data.txt")


# Set up a basic configuration
config = list(iterations=500, modelVariant="betweenItem", maxCategories=15)


# Create the factors data.frame. There are two factors: "letters" and "numbers".
config$factors = data.frame(letters = c('a',  'a',  'a',  'b',  'b',  'b'),
														numbers = c('1',  '2',  '3',  '1',  '2',  '3'),
														cond =    c('a1', 'a2', 'a3', 'b1', 'b2', 'b3'),
														stringsAsFactors = FALSE)


# The names of the cond (conditions) in factors should match those in the data.
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
mhTuning$pMem = 0.2
mhTuning$pMem_cond = 0.08
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.5
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.12


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


# Convergence looks pretty fast, so remove only 500 burn-in iterations
results = removeBurnIn(results, 500)


#####################
# Examine the results

# Because of figure legends in factorial designs, you need a large 
# plotting surface, so plot to a pdf.
plotParameterSummary(results, asPdf = TRUE)


# Test main effects and interactions
mei = testMainEffectsAndInteractions(results, subsamples = 20)

# Examine bayes factors in favor of the effect.
mei[ mei$bfType == "10", ]

# Select tests with reasonably strong BFs (either for or against an effect).
mei[ mei$bf > 3, ]


# You can also compare individual conditions, but there are many, many pairs 
# of conditions to compare, so it's typically better to look at main effects
# and interactions to determine what effects are present.
testConditionEffects(results)



#####################
# We (should) find an interaction between letters and numbers for contSD. 
# Because of the interaction, we can't interpret the main effects of letters or numbers.
# One analysis strategyso is to perform tests of the effect of letters for each level 
# of numbers (which is called simple effects in ANOVA).

subRes = results # Copy the results so you can freely modify it.

# Pick a level of numbers to use
newFactors = subRes$config$factors
newFactors = newFactors[ newFactors$numbers == 2, ] # <<< pick a level
subRes$config$factors = newFactors

# Run a test of letters for contSD.
testSingleEffect(subRes, "contSD", "letters")

# Are the results reasonable when compared to a plot of the parameters?
plotParameterLineChart(subRes, "contSD")

#######################
# Let's say that we want to know which levels of the numbers factor are 
# different for the pMem parameter. We want to do pairwise comparisons of
# factor levels. This is different than what is done by testConditionEffects,
# which does pairwise comparisons of conditions (cells), not of factor levels.
subRes = results

# Drop the letters factor (because we are collaping across it)
subRes$config$factors$letters = NULL

# Pick the levels of the numbers factor to use
usedFactorLevels = data.frame(numbers = c(1,2)) 

testSingleEffect(subRes, "pMem", "numbers", usedFactorLevels = usedFactorLevels)




# Plot the posterior predictive distribution to visually examine model fit.
# There are many conditions, so only plot two of them
posteriorPredictivePlot(results, results$pnums, conditions=c("a1", "b3"), alpha=0.1)

# Get posterior means and credible intervals for conditions.
posteriorMeansAndCredibleIntervals(results)


####################
# Examine the success of the parameter estimation

source("../DataSimulatingFunctions.R")

trueParam = read.delim("factorial_BI_parameters.txt")

comp = compareTrueAndRecovered(results, trueParam)
whichRound = c("true", "rec", "cor", "slope", "dif", "percentDif")
comp[ , whichRound ] = round(comp[ , whichRound ], 2)
comp




#############
# Working directly with CMBBHT
#
# You can do customized hypothesis tests by working directly with the CMBBHT package.
# See the documentation for that package for more information about using it.
library(CMBBHT)

fact = results$config$factors

# The factors data.frame used by CMBBHT cannot have extra columns
# beyond the real factors of the design, so remove the cond column.
fact$cond = NULL 

# Get prior and posterior condition effects
condEff = getConditionEffects(results, "contSD")
names(condEff)

# Check that condition names align properly between fact and colKeys
all(fact$cond == condEff$colKeys$cond)


# A function used to perform hypothesis tests
testFun = function(prior, post) {
	
	I_M0 = function(eff) {
		# The null hypothesis/model is that all absolute 
		# effect parameters are less than 2.
		all(abs(eff) < 2)
	}
	
	CMBBHT::testFunction_EPA(prior, post, I_M0)
}

CMBBHT::testHypothesis(condEff$prior, condEff$post, fact, "numbers", testFunction = testFun)


# Continue working with CMBBHT to directly analyze 
# the posterior effect parameters

# Option 1: Use CMBBHT to get effect parameters
cleanFact = fact
cleanFact$cond = NULL # No extra columns
postEP = CMBBHT::getEffectParameters(condEff$post, cleanFact, "numbers")

# Option 2: Use getMEIParameters helper function
ep = getMEIParameters(results, "contSD", "numbers")
postEP = ep$post

# Analyze the parameters
head(postEP)
summary = CMBBHT::summarizeEffectParameters(postEP)
summary
CMBBHT::plotEffectParameterSummary(summary)


