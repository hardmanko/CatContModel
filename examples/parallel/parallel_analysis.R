
setwd("D:/Programming/R/CatContModel/examples/parallel") #or wherever you are working

library(CatContModel)

# This example uses the same data and parameters as the betweenItem example.
data = read.delim("betweenItem_data.txt")


# TODO
#
# This example uses the ZL model variant to demonstrate parallel features because it is fast, 
# but parallel stuff is mostly applicable to the slow model variants.
# All model variants work the same with regard to running in parallel.
#
config = list(modelVariant="betweenItem", iterations=50)
config$modelVariant = "ZL"

# Start by running the model with runParameterEstimation in the main R process to make sure everything is working.
# It is harder to get output and do debugging when running in parallel.
# Use a small number of iterations for testing.
startRes = runParameterEstimation(config, data)


# MH tuning can be done in the main process before switching to parallel, 
# You'll need at least 500 iterations for useful MH tuning.
# See other examples for more details about MH tuning.
examineMHAcceptance(startRes)


#####################
# Parallel Parameter Estimation
# Once the configuration is tested and working, the main estimation run can be done in parallel using runParameterEstimation_parallel.
# See that function's documentation for more configuration options. Only the most important options are shown here.

# Choose how many parallel nodes will run at once.
# By default, these nodes are run on the local machine in separate R processes.
# See the "cluster" argument for a way to use a cluster that you set up yourself.
numNodes = 4


# Choose a file (or files) to store output from the parallel nodes.
# The files are overwritten each time a new run starts.
# File names can have a directory in front (which must already exist)
outputFile = NULL    # Empty ("") or NULL outputFile means do not store output.
outputFile = "progress/runProgress.txt"  # File name with .txt extension means to save output from one node only. (All nodes should finish at similar times.)
outputFile = "progress/runProgress_"     # An incomplete file name (no .txt extension) means to make a file for each node. (Named like runProgress_1.txt.)


# Set the number of iterations per node (not total iterations).
config$iterations = 3000

# Run parallel parameter estimation. 
resList = runParameterEstimation_parallel(numNodes=numNodes, 
                                         config = config, 
                                         data = data,
                                         outputFile = outputFile)

# resList is a list with one element per node. Each element is the same value as returned by runParameterEstimation.
length(resList)
names(resList[[1]])

# The seed that was used, whether provided by the user or sampled, is in $info for each node's element.
resList[[1]]$info



#####################
# continueSampling

# If you would like to sample additional iterations from the parallel results, use continueSampling_parallel.
# Provide the results list as the first argument.
continueRes = continueSampling_parallel(resList, iterations=1000, outputFile = "progress/continueProgress_")

# The return value from continueSampling_parallel has the same format as from runParameterEstimation_parallel:
length(continueRes)
names(continueRes[[1]])



#####################
# Parallel results lists can be saved and loaded as RDS files:
# saveRDS(continueRes, "ParallelResults_BI.RDS")
# continueRes = readRDS("ParallelResults_BI.RDS")
resList = continueRes

#####################
# Pre burn-in convergence


convergencePlots_parallel(resList, "convergence_BI.pdf")

convergencePlots_parallel(resList, "convergence_BI_200_500.pdf", whichIterations = 200:500)


plot(resList[[1]]$posteriors$`pMem_cond[2]`, type='l')



#####################
# To convert the parallel results to a single, usable results object, you must:

# 1: Remove burn in iterations from each results object.
postBurnIn = list()
for (i in 1:length(resList)) {
  postBurnIn[[i]] = removeBurnIn(resList[[i]], 1000)
}

# 2: Merge the post-burn-in results objects.
usableRes = mergeResults(resList=postBurnIn)


# Analyze the merged results as you normally would with a non-parallel results object.
posteriorMeansAndCredibleIntervals(usableRes)
testConditionEffects(usableRes)
# etc...



