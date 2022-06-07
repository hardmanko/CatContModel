
setwd("D:/Programming/R/CatContModel/examples/betweenItem") #or wherever you are working

library(CatContModel)

# Read in your data.
data = read.delim("betweenItem_data.txt")
head(data)

# Set up a basic configuration. More iterations can be sampled later.
config = list(modelVariant="betweenItem", iterations=500, maxCategories=15)

# Check your config and see default values for settings you didn't provide.
verifyConfigurationList(config, data)


########################
# Parameter Estimation and MH Tuning:
# Nearly all parameters are sampled using the Metropolis-Hastings (MH) method.
# It is desirable, but not necessary, for the MH acceptance rates to be moderate (near 0.5).
# MH tuning is done with 3 steps.

# 1. Override MH tuning defaults
mhTuning = list()
mhTuning$catSelectivity = 3
mhTuning$contSD = 2
mhTuning$contSD_cond = 1.2
mhTuning$pContBetween = 0.4
mhTuning$pContBetween_cond = 0.2
mhTuning$pMem_cond = 0.12


# 2. Run with those MH tuning values
results = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)

# 3. Examine MH acceptance rates. Go back to step 1 until the acceptance rates look good.
examineMHAcceptance(results)


# Once the MH tuning is good, sample more iterations
continueResults = continueSampling(results, 3000)
results = continueResults$combinedResults

# And save the results. I always save results before removing burn-in iterations.
# saveRDS(results, file="betweenItem_results.RDS")

# So that you can read them back in later
results = readRDS("betweenItem_results.RDS")


########################
# Examine convergence by plotting parameter chains to determine how many burn-in iterations to remove.
# It is good to check the slowest parameters to converge.


# Condition effect parameters are often slowest to converge.
plot(results$posteriors[[ "pMem_cond[2]" ]], type='l')
plot(results$posteriors[[ "pContBetween_cond[2]" ]], type='l')


# It takes a couple of steps to plot catActive in a useful way.
# catActive is slow to converge because it depends on the catMu finding their way to the category locations.
post = convertPosteriorsToMatrices(results, "catActive")
activeCats = apply(post$catActive, 3, mean) * results$config$maxCategories
plot(activeCats, type='l')


# Make a huge pdf with convergence plots for all parameter chains
convergencePlots(results, "Convergence_BetweenItem.pdf")

########################
# Remove burn-in iterations
# Convergence is usually reasonably fast, so I typically do 500 to 1000 burn-in iterations
results = removeBurnIn(results, 500)


########################
# Double-check convergence with the Geweke diagnostic.
# Geweke is done after removing burn-in iterations.
library(coda)

pmat = convertPosteriorsToMatrix(results)
pmat = mcmc(pmat) #convert to coda format

gr = geweke.diag(pmat, 0.2, 0.5)
qqnorm(gr$z) #the z-scores should follow a standard normal distribution.
abline(0,1)

geweke.plot(pmat) # Another way of plotting the Geweke diagnostic


########################
# After checking for convergence and removing burn-in, examine the results.
# These are standard analyses that most researchers would do.

# Add a color generating function for plots
results$colorGeneratingFunction = function(angle) {
	hsv((angle / 360) %% 1, 1, 1)
}

# Plot summaries of the parameters (see asPdf argument if plot is hard to read)
plotParameterSummary(results)


# Plot individual parts of the overall summary
plotParameterLineChart(results, "pMem")
plotParameterHistogram(results, "catSD")
plotParameterHistogram(results, "catActive")
plotCatMu(results)


# Numerical summary of parameters
posteriorMeansAndCredibleIntervals(results)


# Hypothesis tests of differences between task conditions (pairwise tests between levels of data$cond).
testConditionEffects(results)


# Hypothesis tests of main effects and interactions
mei = testMainEffectsAndInteractions(results, subsamples = 20)
mei[ mei$bfType == "10", ]


# Test model validity by comparing data sampled from the fitted model to real data.
posteriorPredictivePlot(results, "1", alpha=0.3)
posteriorPredictivePlot(results, "5", alpha=0.3)



##############################
# The data for this example was generated from the betweenItem model.
# Test model variant specificity by testing if the ZL model fits the data well.
# If ZL fits as well as betweenItem, then betweenItem is not capturing enough 
# of a meaningful pattern in the data to be a useful model.
# With real data, researchers might want to check that the simpler ZL model does not adequately describe their data.

zlConfig = list(iterations=3500, modelVariant="ZL", iterationsPerStatusUpdate = 200)

zlMh = list()
zlMh$contSD = 1.5
zlMh$contSD_cond = 0.7
zlMh$pMem = 0.2
zlMh$pMem_cond = 0.1

zlResults = runParameterEstimation(zlConfig, data, mhTuningOverrides = zlMh)

examineMHAcceptance(zlResults)

# I don't show convergence diagnostics here, but they are technically relevant. 
# In my experience, ZL converges very fast, so I didn't bother.
zlResults = removeBurnIn(zlResults, 500)

plotParameterSummary(zlResults)

# You can compare the fit of the ZL and betweenItem models with WAIC
# The data were generated under the betweenItem model, so you would expect it to win.
biWaic = calculateWAIC(results)
zlWaic = calculateWAIC(zlResults)

# The WAIC values of the betweenItem model should be lower than those of the ZL model.
biWaic
zlWaic


### 
# Model comparison with WAIC is the best way to select a model.
# However, the following function does some tests to see if
# categories are being used by participants.

testCategoricalResponding(results)


##############################
# Examine the success of the parameter estimation vs known true values that were used to simulate the data.
source("../DataSimulatingFunctions.R")

trueParam = read.delim("betweenItem_parameters.txt")

comp = compareTrueAndRecovered(results, trueParam, roundDigits = 2)
comp

median(abs(comp$percentDif))
