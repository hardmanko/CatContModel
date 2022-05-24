# Half Circle example
#
# This example shows two ways to deal with tasks in which participants 
# remember the orientation of a line (or bar) that rotates around a central point. 
# The line sticks out in both directions from the center, so it points in two directions at once, in a sense.
# For example, a bar that is 30 degrees to the right of center at the top is
# 30 degrees left of center at the bottom. 
# Thus, participants measure the orientation of the line in a range from 0 to 180 degrees, a half circle.
#
# Two solutions are explored:
# Solution 1. Treat the data are linear, constrained between 0 and 180 degrees.
# Solution 2. Multiply the data by 2 so that they fill a full circle, 0 to 360 degrees.
#
# For both solutions, the method of analysis is shown with comments.
# At the end, the accuracy of the recovered parameter estimates is  
# compared with the true parameters that were used to sample the data.

setwd("D:/Programming/R/CatContModel/examples/halfCircle_betweenItem") #or wherever you are working

library(CatContModel)


# See halfCircle_simulateData.R for notes on how the data were sampled.
data = read.delim("halfCircle_data.txt")

head(data)

# The data are in a half-circle [0,180]
range(data$study)
range(data$response)

# Data scatterplot
plot(data$study, data$response, pch=20, col=rgb(0,0,0,0.3))

# Response histogram
hist(data$response, breaks=100)


########################
# Solution 1: Treat data as linear.
# For CatContModel, linear data have endpoints and the data do not wrap around.
# The big problem with Solution 1 is that data in a half circle do not have endpoints.
# Solution 1 will produce results that are distorted in ways that may be hard to see,
# but the analysis is very straightforward.

dataLin = data # copy the data before modification

# Try to move the typically observed orientation categories away from the edge (0 and 180).
# There are usually 8 categories in a whole circle, so rotate 1/16th of a circle to
# be between those 8 categories. 360/16 = 22.5 degrees.
# In real data, the assumption that there are 8 categories should be double checked by
# plotting a histogram of responses (with many bins).
rotateAmount = 360/16

# Do %% 180 to move values that were rotated past 180 back to a base of 0.
dataLin$study = (dataLin$study + rotateAmount) %% 180
dataLin$response = (dataLin$response + rotateAmount) %% 180

# Check the ranges and plot the data
range(dataLin$study)
range(dataLin$response)
plot(dataLin$study, dataLin$response)
hist(dataLin$response, breaks=100)

# If treating data as linear, set responseRange and catMuRange to the half circle space.
configLin = list(modelVariant="betweenItem", iterations=500, maxCategories=8,
                dataType="linear", 
                responseRange=c(0, 180), 
                catMuRange=c(0, 180))

# Run the model
resLin = runParameterEstimation(configLin, dataLin)

# MH tuning is not usually needed, but double check the acceptance rates just in case there is a big problem
examineMHAcceptance(resLin) 


continueResults = continueSampling(resLin, 3000)
resLin = continueResults$combinedResults

# saveRDS(resLin, file="halfCircle_results_linear.RDS")
resLin = readRDS("halfCircle_results_linear.RDS")


# Brief convergence analysis and remove burn in
plot(resLin$posteriors[[ "pMem_cond[3]" ]], type='l')
plot(resLin$posteriors[[ "pContBetween_cond[2]" ]], type='l')

resLin = removeBurnIn(resLin, 500)


###
# When treating half-circle data as linear, analyze the results normally.

plotParameterSummary(resLin)

posteriorMeansAndCredibleIntervals(resLin)

mei = testMainEffectsAndInteractions(resLin)
mei[ mei$bfType == "10", ]

ce = testConditionEffects(resLin)
ce[ ce$bfType == "10", ]

posteriorPredictivePlot(resLin)


########################
# Solution 2: Modify the data by multiplying by 2, changing the range from [0,180] to [0,360]
#
# This solution is good in that it treats the data as circular, just on a half-size circle.
#
# This solution has downsides when it comes to describing the results.
# In particular, because the range of the data has been increased by a factor of 2,
# the standard deviation parameters (contSD, catSD, and catSelectivity) and catMu
# will be double their true values.
#
# For plots and other descriptive statistics, these parameters must be scaled by dividing by 2.
# For hypothesis tests, however, the original, undivided parameter values must be used.
# The code below gives more detail.

data2x = data # copy the data before modification

# Modify the data by multiplying by 2
data2x$study = data2x$study * 2
data2x$response = data2x$response * 2

# Check the ranges and plot the data
range(data2x$study)
range(data2x$response)
plot(data2x$study, data2x$response)

# Set up a basic configuration and run the model
config2x = list(modelVariant="betweenItem", iterations=500, maxCategories=8)

res2x = runParameterEstimation(config2x, data2x)

examineMHAcceptance(res2x)

continueResults = continueSampling(res2x, 3000)
res2x = continueResults$combinedResults

# saveRDS(res2x, file="halfCircle_results_2x.RDS")
res2x = readRDS("halfCircle_results_2x.RDS")


# Brief convergence analysis and remove burn in
plot(res2x$posteriors[[ "pMem_cond[3]" ]], type='l')
plot(res2x$posteriors[[ "pContBetween_cond[2]" ]], type='l')

res2x = removeBurnIn(res2x, 500)



#####
# Descriptive Analysis
# To describe the standard deviation and catMu parameter posteriors with respect to the 
# original half-circle data space, those parameter posteriors must be divided by 2.

# This function to convert posterior parameter chains to/from half circle space
# The standard deviation parameters (contSD, catSD, and catSelectivity) and catMu 
# posteriors will be multiplied by multiple.
halfCircle_convertParameters = function(res, multiple) {

  fullParamNames = names(res$posteriors)
  
  # Select standard deviation parameters from all parameters
  sdParamNames = getSdParams(res)
  fullSdParam = NULL
  for (sdp in sdParamNames) {
    matchSd = grepl(sdp, fullParamNames, fixed=TRUE)
    fullSdParam = c(fullSdParam, fullParamNames[ matchSd ])
  }
  
  # Modify standard deviation parameters
  for (sdp in fullSdParam) {
    
    post = res$posteriors[[ sdp ]]
    
    if (grepl(".var", sdp, fixed=TRUE)) {
      # Variance parameters are in squared units, so sqrt, divide, then square
      post = sqrt(post)
      post = post * multiple
      post = post^2
      # Adjusting variance parameters like this may not be totally correct, but approximate.
    } else {
      # Participant and mean (.mu) parameters divide by 2
      post = post * multiple
    }
    
    res$posteriors[[ sdp ]] = post
  }
  
  # Modify catMu
  fullCatMu = fullParamNames[ grepl("catMu", fullParamNames, fixed=TRUE) ]
  for (cm in fullCatMu) {
    
    post = res$posteriors[[ cm ]]
    
    # The catMu parameters are free to wander outside of [0,360].
    # Bring them back before multiplying
    post = CatContModel::clampAngle(post, pm180 = FALSE, degrees = TRUE)
    
    post = post * multiple
    
    res$posteriors[[ cm ]] = post
  }

  res
}

# res2x was estimated with data that were multiplied by 2 to get to the full circle space.
# Get back to the half circle space by dividing by 2 (i.e. multiply by 1/2)
res2x_conv = halfCircle_convertParameters(res2x, multiple=1/2)

# Means and credible intervals are descriptive and use converted parameters
posteriorMeansAndCredibleIntervals(res2x_conv)

# Plots are also descriptive and use converted parameters.
# The catMu plot has some problems, but is still readable
plotParameterSummary(res2x_conv)

# You can plot individual panels of the overall plot with plotParameter
plotParameter(res2x_conv, param="contSD")
plotParameter(res2x_conv, param="catSelectivity")
plotParameter(res2x_conv, param="catSD")


# This function is a possible way to plot catMu for half circle data
plotCatMu_halfCircle = function(res, precision = 1, pnums = res$pnums) {
  
  post = convertPosteriorsToMatrices(res, param = c("catMu", "catActive"))
  
  catMu = post$catMu[ pnums, , ]
  catActive = post$catActive[ pnums, , ]
  
  plotData = CatContModel:::catMu_getPlotData(catMu, catActive, precision = precision, 
                                              dataType = res$config$dataType, 
                                              responseRange = res$config$responseRange)
  
  # Ignore all points past 180
  plotData = plotData[ plotData$x <= 180, ]
  
  # The last point is the first point, so copy the value from the first point.
  plotData$y[ plotData$x == 180 ] = plotData$y[1]
  
  CatContModel:::catMu_plotSetup(plotData, res$config$dataType)
  CatContModel:::catMu_plotData(plotData$x, plotData$y, col = plotData$color, type = "polygon")
  
  invisible(plotData)
}

plotCatMu_halfCircle(res2x_conv)


###
# Although converted parameters are used for describing the results, other kinds of 
# analysis will not work correctly if you use the converted parameter values.
# Here are some major examples of where you should use the original parameter values.

# Hypothesis tests depend on priors, so use unconverted parameters
mei = testMainEffectsAndInteractions(res2x)
mei[ mei$bfType == "10", ]

ce = testConditionEffects(res2x)
ce[ ce$bfType == "10", ]

# Posterior predictive (sampling data based on parameters) uses unconverted parameters
posteriorPredictivePlot(res2x)

# WAIC also uses unconverted parameters
calculateWAIC(res2x)

# There are many functions in this package that summarize results.
# When using each, think carefully about whether the function wants
# the original real posteriors or converted posteriors.
# Anything that uses priors must use original posteriors.




##############################
# Examine the success of parameter estimation for both solutions.

source("../DataSimulatingFunctions.R")

# The true parameter values are the values that were used to generate the data
# in the full circle space. So the "true" parameters are wrong by 2x with
# respect to the data that have been divided by 2 to fit in the half circle space.
# This is confusing, but I didn't want to modify the parameters 
# that were actually used to generate the data.

trueParam = read.delim("halfCircle_parameters.txt")


###
# Solution 1: Linear
# The linear results are in the [0,180] space while the "true" parameters are in the [0,360] space.
# It's predictable that the recovered are about 50% of the "true".
compareTrueAndRecovered(resLin, trueParam, roundDigits=2)

# For more direct comparison of parameters, convert the linear parameters by 2x.
# This conversion is only done to compare to "true" parameter values in this example.
# For normal analysis of resLin, no conversion is needed. 
resLinConverted = halfCircle_convertParameters(resLin, multiple=2)

trLinConv = compareTrueAndRecovered(resLinConverted, trueParam, roundDigits=2)
trLinConv


###
# Solution 2: Compare the original parameters (not converted) to true parameters.
tr2x = compareTrueAndRecovered(res2x, trueParam, roundDigits=2)
tr2x
# It is unsurprising if the parameters are recovered fairly accurately.
# The data were divided by 2 in halfCircle_simulateData.R, then multiplied by 2 in this file.
# Thus, the data and true parameters line up exactly.


###
# Lazy aggregation of error amounts between solutions.
# Solution 2 tends to be more accurate, but not by a lot.
mean(abs(trLinConv$percentDif))
mean(abs(tr2x$percentDif))

median(abs(trLinConv$percentDif))
median(abs(tr2x$percentDif))

# The aggregate error has been moderate in the runs I've done (between 5 and 10 percent) for both solutions.
# This is a tolerable amount of error, but slightly higher than other examples (like the betweenItem example).

