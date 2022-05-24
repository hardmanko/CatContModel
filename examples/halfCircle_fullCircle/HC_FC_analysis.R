# Half-Circle Full-Circle example
# 
# This example uses the Zhang & Luck (2008) model variant ("ZL") 
# with data from two tasks:
#
# 1. A half-circle task with data in [0,180] degrees.
# 2. A full-circle task with data in [0,360] degrees.
#
# There are two conditions in each task for a 2x2 design.
# One of the conditions can be dropped from each task to simplify the design.
#
# To ensure maximum comparability between tasks, the following approach is used:
# 1. Convert study:response pairs to response error (response - study):0 pairs.
# This can be done since the ZL model has no categories and can be fitted with just response error.
# 2. Treat response error in both tasks as linear, with bounded endpoints. Since the endpoints are in
# the tails of the response error distribution, treating data as linear rather than circular should 
# not bias the results.


setwd("D:/Programming/R/CatContModel/examples/halfCircle_fullCircle") #or wherever you are working

library(CatContModel)

# See HC_FC_simulateData.R for notes on how the data were sampled.
hcData = read.delim("halfCircle_data.txt")
fcData = read.delim("fullCircle_data.txt")


plot(hcData$study, hcData$response)
plot(fcData$study, fcData$response)


# Multiply half circle by 2 just so the circular distance function works correctly.
# The data are divided by 2 before fitting.
hcData$study = hcData$study * 2
hcData$response = hcData$response * 2

# Calculate signed response error
hcData$error = CatContModel::circDist(hcData$study, hcData$response)
fcData$error = circDist(fcData$study, fcData$response)


# Divide the half circle data, including error, by 2 to return to original space
hcData$study = hcData$study / 2
hcData$response = hcData$response / 2
hcData$error = hcData$error / 2

# Put the centered data into new data frames for fitting.
# Set study to 0 to reflect that response error has had study subtracted from it.
hcDataCentered = data.frame(pnum=hcData$pnum, cond=hcData$cond, study=0, response=hcData$error)
fcDataCentered = data.frame(pnum=fcData$pnum, cond=fcData$cond, study=0, response=fcData$error)

hist(hcDataCentered$response)
hist(fcDataCentered$response)

range(hcData$response)
range(fcData$response)

# For a simpler design, use only one of the conditions from each task.
#hcData = hcData[ hcData$cond == 1, ]
#fcData = fcData[ fcData$cond == 1, ]


# Make half circle config treating the data as linear.
hcConfig = list(modelVariant = "ZL", 
                iterations = 5000,
                dataType = "linear", 
                responseRange = c(-90,90))


# Full circle config is the same except for response range
fcConfig = hcConfig
fcConfig$responseRange = c(-180, 180)

# Run the model
hcRes = runParameterEstimation(hcConfig, hcDataCentered)
fcRes = runParameterEstimation(fcConfig, fcDataCentered)

# Do MH tuning if desired
examineMHAcceptance(hcRes)
examineMHAcceptance(fcRes)

# Do convergence here analysis if desired

# For this example, burn in assumed to be fast for ZL model
hcRes = removeBurnIn(hcRes, 500)
fcRes = removeBurnIn(fcRes, 500)


###
# Correlation analysis
# Are the participant-level parameters correlated between tasks?

# Get posterior means of participant parameters
hcPP = participantPosteriorSummary(hcRes)
fcPP = participantPosteriorSummary(fcRes)

# For correlations, just use one of the conditions because different conditions 
# just have the conditions effect offset from one another.
hcPP = hcPP[ hcPP$cond == 1, ]
fcPP = fcPP[ fcPP$cond == 1, ]

# Correlations between posterior means.
# You probably want to use a Bayesian correlation test here rather than something like cor.test.
# This just calculates the correlation coefficient.
cor(hcPP$mean[ hcPP$param == "pMem" ], fcPP$mean[ fcPP$param == "pMem" ])
cor(hcPP$mean[ hcPP$param == "contSD" ], fcPP$mean[ fcPP$param == "contSD" ])


###
# For other analyses, make a between participants results object 
# with combineGroupResults.BP then analyze that object.

bpRes = combineGroupResults.BP(list(hc=hcRes, fc=fcRes))

# Using default factor names:
# Factor is condition within task and
# BP_Group is task, full-circle or half-circle.
# You can get nice
plotParameterSummary(bpRes)

testMainEffectsAndInteractions(bpRes)

testConditionEffects(bpRes)






##############################
# Examine the success of parameter estimation
#
# The true parameter values used to generate the half circle data were
# in the full circle space. So the "true" parameters are wrong by 2x with
# respect to the data that have been divided by 2 to fit in the half circle space.
# This is confusing, but I didn't want to modify the parameters 
# that were actually used to generate the data.

# Function to convert standard deviation parameters between half circle and full circle space.
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

source("../DataSimulatingFunctions.R")

hcTrueParam = read.delim("halfCircle_parameters.txt")
fcTrueParam = read.delim("fullCircle_parameters.txt")

# Full circle
compareTrueAndRecovered(fcRes, fcTrueParam, roundDigits=2)

# For half circle, contSD is off by a factor of 2
compareTrueAndRecovered(hcRes, hcTrueParam, roundDigits=2)

# Multiply parameters by 2 just to compare to "true" parameter values
hcConv = halfCircle_convertParameters(hcRes, 2)
compareTrueAndRecovered(hcConv, hcTrueParam, roundDigits=2)

