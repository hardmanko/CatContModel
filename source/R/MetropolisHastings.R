
#' Examine the Metropolis-Hastings Acceptance Rates
#' 
#' Most of the parameters in the model do not have conjugate priors, so Metropolis-Hastings (MH) parameter updating steps are used 
#' to sample from their conditional posterior distributions. This function helps to examine the acceptance rate of the 
#' MH steps for different parameter groups.
#' 
#' When a category is inactive, `catMu` is always accepted. 
#' Thus, when considering the acceptance rate of the `catMu` parameters, you should only consider the acceptances of active categories. 
#' This function gives two acceptance rates for `catMu`: 
#' Active only, which you should use, and the total, including inactive categories, which is basically irrelevant.
#' 
#' When `recalculate=FALSE`, the results of this function are not changed by removing burn-in iterations, except that the `catMu` 
#' acceptance rates are recalculated in a way that makes them meaningless once burn-in iterations have been removed. 
#' In short: Only use this function before burn-in iterations are removed or use `recalculate=TRUE`.
#' 
#' @param results The results from the [`runParameterEstimation`] function.
#' @param recalculate If `FALSE`, values from the Gibbs sampler are used. If `TRUE`, acceptance rates are recalculated from posteriors.
#'
#' @return A `data.frame` containing a summary of the MH acceptance rates for each group of parameters.
#' 
#' @family MetropolisHastings 
#' 
#' @export
examineMHAcceptance = function(results, recalculate=FALSE) {
  
  # TODO: Rename? to summarizeMHAcceptance mhAcceptanceSummary
  
  if (!resultIsType(results, "WP")) {
    stop("Can only examine MH acceptance of WP results objects.")
  }

  if (recalculate) {
    
    acceptRates = calculateMHAcceptance(results)
    
  } else {
    # Use values returned by Gibbs sampler
    acceptRates = results$MH$acceptance
    
    #separate out active catMu acceptances.
    for (pnum in results$pnums) {
      
      for (cat in 1:results$config$maxCategories) {
        
        cmn = paste("catMu[", pnum, ",", cat, "]", sep="")
        can = paste("catActive[", pnum, ",", cat, "]", sep="")
        
        if (cmn %in% names(results$constantValueOverrides) || !(cmn %in% names(results$posteriors))) {
          next #skip constant catMu and, if somehow the catMu is not in the posteriors, skip it as well
        }
        
        totalPossibleAcceptances = results$runConfig$iterations
        
        nActive = sum(results$posteriors[[ can ]])
        nInactive = totalPossibleAcceptances - nActive
        
        totalAcceptances = results$MH$acceptance[ results$MH$acceptance$name == cmn, ]$acceptanceRate * totalPossibleAcceptances
        activeAcceptances = totalAcceptances - nInactive
        activeAcceptanceRate = activeAcceptances / nActive
        if (is.nan(activeAcceptanceRate) || is.na(activeAcceptanceRate) || activeAcceptanceRate < 0) {
          activeAcceptanceRate = 0
        }
        
        temp = data.frame(name=paste(cmn, "(Active)"), group="catMu (Active only)", acceptanceRate = activeAcceptanceRate)
        
        acceptRates = rbind(acceptRates, temp)
        
      }
    }
  }
  
  
  summarizeAcceptance = function(x) {
    qs = as.vector(stats::quantile(x, c(0, 0.025, 0.5, 0.975, 1)))
    qs = c(mean(x), qs)
    qs
  }

  summ = stats::aggregate(acceptanceRate ~ group, acceptRates, summarizeAcceptance)
  summ = do.call(data.frame, summ)
  
  names(summ) = c("paramGroup", "mean", "min", "2.5%", "median", "97.5%", "max")
  summ$paramGroup = as.character(summ$paramGroup)
  
  summ$paramGroup[ summ$paramGroup == "catMu" ] = "catMu (Total, ignore)"
  
  summ = summ[ order(summ$paramGroup), ]
  
  summ
}

calculateMHAcceptance = function(results) {
  
  acceptRates = NULL
  
  # All parameters are MH parameters except for P.mu and P.var.
  for (pn in names(results$posteriors)) {
    
    splitName = splitParamName(pn)
    
    # Ignore "so-called" "hyper" parameters (because not MH)
    if ("hyper" %in% splitName$types) {
      next
    }
 
    # Do fancy for active catMu. find sections of active catMu and calculate acceptance rate for each section.
    if (splitName$baseName == "catMu") {
      
      caName = paste0("catActive[", splitName$index[1], ",", splitName$index[2], "]")
      
      cm = results$posteriors[[pn]]
      ca = results$posteriors[[caName]]
      
      if (all(ca == 0)) {
        next # skip totally inactive
      }
      if (all(cm == cm[1])) {
        next # skip constant
      }
      
      # Find starts and ends of active sections
      allFirst1 = 1 + which((ca[1:(length(ca)-1)] != ca[2:length(ca)]) & (ca[2:length(ca)] == 1))
      if (ca[1] == 1) {
        allFirst1 = c(1, allFirst1)
      }
      allLast1 = which((ca[1:(length(ca)-1)] != ca[2:length(ca)]) & (ca[2:length(ca)] == 0))
      if (ca[length(ca)] == 1) {
        allLast1 = c(allLast1, length(ca))
      }
      
      if (length(allFirst1) != length(allLast1)) {
        stop("Unequal catActive active section lengths")
      }
      
      whichEq = (allFirst1 == allLast1)
      allFirst1 = allFirst1[ !whichEq ]
      allLast1 = allLast1[ !whichEq ]
      
      acceptSum = 0
      possibleSum = 0
      if (length(allFirst1) >= 1) {
        for (i in 1:length(allFirst1)) {
          
          sacm = cm[ allFirst1[i]:allLast1[i] ]
          
          accepted = sacm[2:length(sacm)] != sacm[1:(length(sacm) - 1)]
  
          acceptSum = acceptSum + sum(accepted)
          possibleSum = possibleSum + length(accepted)
  
        }
      }
      
      rate = acceptSum / possibleSum
      if (possibleSum == 0) {
        rate = 0
      }
      
      temp = data.frame(name = pn, group = "catMu (Active only)", acceptanceRate = rate)
      
      acceptRates = rbind(acceptRates, temp)
    }
    
    # Calcaulate acceptance rate.
    # Even if this parameter is catMu, do this to get acceptance for both active and inactive cats.
    pv = results$posteriors[[pn]]
    pv = c(results$startingValues[[pn]], pv)
    
    if (all(pv == pv[1])) {
      # Skip constant parameters
      next
    } else {
      accepted = pv[2:length(pv)] != pv[1:(length(pv) - 1)]
      
      rate = mean(accepted)
    }

    
    group = splitName$baseName
    if ("cond" %in% splitName$types) {
      group = paste0(group, "_cond")
    }
    if (splitName$baseName == "catMu") {
      group = "catMu (Total, ignore)"
    }
    
    temp = data.frame(name = pn, group = group, acceptanceRate = rate)
    
    acceptRates = rbind(acceptRates, temp)

  }
  
  acceptRates
}


# Helper to extract MH tuning values from targets data frame
getMHTuningList = function(targets, mhTuningOverrides=list()) {
  
  rval = list()
  for (i in 1:nrow(targets)) {
    pn = targets$parName[i]
    if (!is.null(mhTuningOverrides[[pn]])) {
      rval[[pn]] = mhTuningOverrides[[pn]]
      
    } else if (!is.na(targets$constantTuning[i])) {
      rval[[pn]] = targets$constantTuning[i]
      
    } else {
      rval[[pn]] = targets$currentTuning[i]
    }
  }
  rval
}



#' Make Metropolis-Hastings Configuration List
#' 
#' Configures optimization of Metropolis-Hastings tuning values. 
#' Optimization is done by [`optimizeMHTuning`], which is used by [`runParameterEstimation`] if the `mhOptim` argument enables MH optimization.
#' 
#' @param optimSamples Number of samples of `optimSampleIter` iterations to take while optimizing MH tuning values.
#' @param optimSampleIter Number of iterations to sample between updates to the MH tuning values.
#' @param optimRate Modifies how rapidly MH tuning parameters are updated. Values greater than 1 result in faster optimization.
#' @param targetAcceptance The target MH acceptance rate for parameters. See details for more.
#' @param chosenTuning When optimization is complete, final tuning values are chosen based on `chosenTuning`. See details.
#' @param startingTuningValues A named list of starting values for MH tuning values. Defaults to [`getDefaultMHTuning`].
#' @param constantTuningValues A named list of constant values for MH tuning values. See also the `mhTuningOverrides` argument of [`makeModelConfig`].
#' @param modelVariant If provided, only parameters used by the `modelVariant` will be in the returned object.
#' 
#' @details 
#' 
#' `targetAcceptance` can be a list that maps from parameter name (like "pMem" or "contSD_cond") to target acceptance 
#' rate for that parameter. If listing target acceptance rates, include a element named "default" with the default 
#' acceptance rate for unlisted parameters. Alternately, set targets for parameters by modifying the returned `targets` data frame.
#' 
#' `chosenTuning` specifies which tuning values are automatically selected when doing MH optimization through [`runParameterEstimation`].
#' If using [`optimizeMHTuning`], all sets of tuning values can be examined and/or selected.
#' 
#' The options for `chosenTuning` are:
#' + `"last"`: The MH tuning values that were updated following the last sample.
#' + `"best"`: For each parameter type (listed in `rval$targets`), the MH tuning value used for the best sample for that parameter (lowest difference from target). 
#' + `"bestHalfMean"`: For each parameter type, the mean of the MH tuning values for the best half of samples for that parameter.
#' 
#' @return A configuration list that can be used as the `mhOptimConfig` argument of [`optimizeMHTuning`] or the `mhOptim` argument of [`runParameterEstimation`]. 
#' Values in this list can be changed to configure MH optimization in detailed ways.
#' 
#' @family MetropolisHastings 
#' @export
#' 
#' @examples 
#' targetAcceptance = list(default = 0.4, pMem = 0.6, pMem_cond = 0.4)
#' 
#' startingTuning = list(contSD = 0.8, catSD_cond = 0.3)
#' 
#' constantTuning = list(pCatGuess = 0.8, catSelectivity = 2.5)
#' 
#' mhConfig = makeMHConfig(optimSamples=7, optimSampleIter=70,
#' chosenTuning = "best",
#' targetAcceptance = targetAcceptance,
#' startingTuningValues = startingTuning,
#' constantTuningValues = constantTuning,
#' modelVariant = "betweenItem"
#' )
makeMHConfig = function(optimSamples=10, optimSampleIter=100, optimRate=1,
                        targetAcceptance=0.5, 
                        chosenTuning=c("bestHalfMean", "best", "last"),
                        startingTuningValues=getDefaultMHTuning(),
                        constantTuningValues=NULL,
                        modelVariant=NULL) 
{
  
  if (is.list(targetAcceptance) && is.null(targetAcceptance$default)) {
    targetAcceptance$default = 0.5
  }

  cfg = list(optimSamples = optimSamples,
             optimSampleIter = optimSampleIter,
             optimRate = optimRate, 
             chosenTuning = chosenTuning[1])
  
  defaultMH = getDefaultMHTuning()
  
  # Remove unused parameters here, but need to know modelVariant to do that.
  usedParamNames = names(defaultMH)
  if (!is.null(modelVariant)) {
    modelParamNames = getParamNames(modelVariant=modelVariant, cond=TRUE)
    usedParamNames = modelParamNames[ modelParamNames %in% names(defaultMH) ]
  }
  
  invalidStartingNames = names(startingTuningValues)[ !(names(startingTuningValues) %in% names(defaultMH)) ]
  invalidConstantNames = names(constantTuningValues)[ !(names(constantTuningValues) %in% names(defaultMH)) ]
  
  if (length(invalidStartingNames) > 0) {
  	logWarning("Invalid names in startingTuningValues: ", paste(invalidStartingNames, sep=", "))
  }
  if (length(invalidConstantNames) > 0) {
  	logWarning("Invalid names in constantTuningValues: ", paste(invalidConstantNames, sep=", "))
  }
  
  targetsDF = data.frame()
  
  for (n in usedParamNames) {
    temp = data.frame(parName=n)
    if (grepl("_cond", n, fixed=TRUE)) {
      temp$targetStat = "mean"
    } else {
      temp$targetStat = "median"
    }
    
    if (is.list(targetAcceptance)) {
      temp$targetAcceptance = valueIfNull(targetAcceptance[[n]], targetAcceptance$default)
    } else {
      temp$targetAcceptance = targetAcceptance
    }
    
    temp$startingTuning = valueIfNull(startingTuningValues[[n]], defaultMH[[n]])
    temp$currentTuning = temp$startingTuning
    
    temp$constantTuning = valueIfNull(constantTuningValues[[n]], NA)
    
    # Set current and starting tuning values to constant values.
    if (!is.na(temp$constantTuning)) {
      temp$currentTuning = temp$constantTuning
      temp$startingTuning = temp$constantTuning # constant overrides starting
    } else {
      temp$currentTuning = temp$startingTuning
    }
    
    targetsDF = rbind(targetsDF, temp)
  }
  
  cfg$targets = targetsDF

  cfg
}


# targetStat can be mean, median, or any other column of the MH acceptance data frame
MHOptim_updateTuning = function(res, mhOptimConfig, proportionComplete=0) {
  
  mha = examineMHAcceptance(res)
  
  progress = NULL
  nextTuning = res$MH$tuning # start with previous values
  bestInfo = NULL
  
  for (n in names(res$MH$tuning)) {
    rowName = n
    if (n == "catMu") {
      rowName = "catMu (Active only)"
    }
    
    if (!(rowName %in% mha$paramGroup)) {
      next
    }
    
    targetRow = mhOptimConfig$targets[ mhOptimConfig$targets$parName == n, ]
    
    # Don't update constant values
    if (!is.na(targetRow$constantTuning)) {
      nextTuning[[n]] = targetRow$constantTuning
      next
    }
    
    acceptRow = which(mha$paramGroup == rowName)
    
    acceptanceRate = mha[acceptRow, targetRow$targetStat]
    
    dif = acceptanceRate - targetRow$targetAcceptance
    
    optimRate = mhOptimConfig$optimRate
    if (is.function(optimRate)) {
      optimRate = optimRate(proportionComplete)
    }
    
    # If acceptance is too high (dif +), increase step size
    # If acceptance is too low (dif -), decrease step size
    nextScale = (1 + dif * optimRate)
    
    # Clamp to ensure tuning values don't end up negative or huge
    nextScale = min(max(0.1, nextScale), 1.8)
    
    nextTuning[[n]] = res$MH$tuning[[n]] * nextScale
    
    progTemp = data.frame(parName=n, acceptanceRate=acceptanceRate, dif=dif, scale=nextScale, 
                          oldTuning=res$MH$tuning[[n]], newTuning=nextTuning[[n]])
    progress = rbind(progress, progTemp)
    
    bestTemp = data.frame(parName=n, acceptanceRate=acceptanceRate, dif=dif, tuning=res$MH$tuning[[n]])
    bestInfo = rbind(bestInfo, bestTemp)
    
  }
  
  rval = list(updatedTuning = nextTuning, acceptance=mha, progress=progress, bestInfo=bestInfo)
  
  rval
}





#' Optimize Metropolis-Hastings Tuning Values
#' 
#' Metropolis-Hastings (MH) tuning values affect the quality of posterior parameter chains.
#' Moderate MH acceptance rates are desirable (see [`examineMHAcceptance`]). 
#' 
#' @param data The data to fit.
#' @param modCfg A model configuration. See [`makeModelConfig`].
#' @param optimSamples Optimization argument passed to [`makeMHConfig`].
#' @param optimSampleIter Optimization argument passed to [`makeMHConfig`].
#' @param optimRate Optimization argument passed to [`makeMHConfig`].
#' @param targetAcceptance Optimization argument passed to [`makeMHConfig`].
#' @param chosenTuning Optimization argument passed to [`makeMHConfig`].
#' @param startingTuningValues Optimization argument passed to [`makeMHConfig`].
#' @param mhOptimConfig A configuration list like that returned by [`makeMHConfig`]. If `mhOptimConfig` is provided, all other optimization arguments are ignored.
#' 
#' @details 
#' 
#' This function runs `optimSamples` parameter estimations of `optimSampleIter` each.
#' When each sample is complete, MH acceptance rates are calculated and tuning values from the sample are updated based 
#' on the difference between the acceptance rate for the sample and the target acceptance rate.
#' 
#' The next sample continues the chain from the previous sample with the updated MH tuning values. 
#' This means that the model starts with overdispersed starting values and gradually converges during MH optimization.
#' Ideally, MH tuning values produce moderate MH acceptance rates for the converged model
#' 
#' 
#' @return A list with various components.
#' \tabular{ll}{
#' 	`chosenTuning` \tab The MH tuning values matching the `chosenTuning` (see [`makeMHConfig`]). \cr
#'  `tunings` \tab A list of the tuning values from the `chosenTuning` is taken. \cr
#' 	`optimSteps` \tab A `data.frame` with information about the steps of the optimization procedure. \cr
#'	`sampleResults` \tab A numbered list with one element per `optimSamples`. \cr
#'	`combinedResults` \tab All `sampleResults` combined into a single results object that can be evaluated for convergence or passed to [`continueSampling`].
#' }
#' 
#' @family MetropolisHastings
#' 
#' @export
optimizeMHTuning = function(data, modCfg,
                            optimSamples=10, optimSampleIter=100, optimRate=1,
                            targetAcceptance=0.5, 
                            chosenTuning=c("bestHalfMean", "best", "last"),
                            startingTuningValues=getDefaultMHTuning(), 
                            mhOptimConfig=NULL) 
{
  
  if (is.null(mhOptimConfig)) {
    mhOptimConfig = makeMHConfig(optimSamples = optimSamples, optimSampleIter = optimSampleIter, optimRate = optimRate,
                                 targetAcceptance = targetAcceptance, 
                                 chosenTuning = chosenTuning,
                                 startingTuningValues = startingTuningValues,
                                 constantTuningValues = modCfg$mhTuningOverrides, 
                                 modelVariant = modCfg$modelVariant)
  }
  
  # Set MH overrides from modCfg. 
  # This is primarily for when mhOptimConfig is provided. 
  # If NULL, modCfg$mhTuningOverrides are set by makeMHConfig.
  for (n in names(modCfg$mhTuningOverrides)) {
    mhOptimConfig$targets[ mhOptimConfig$targets$parName == n, c("startingTuning", "currentTuning", "constantTuning") ] = modCfg$mhTuningOverrides[[n]]
  }

  if (mhOptimConfig$optimSamples <= 0 || mhOptimConfig$optimSampleIter <= 0) {
  	logWarning("No MH optimization possible: optimSamples <= 0 or optimSampleIter <= 0.")
    rval = list(chosenTuning=getMHTuningList(mhOptimConfig$targets))
    return(rval)
  }
  
  if (all(!is.na(mhOptimConfig$targets$constantTuning))) {
  	logWarning("All MH tuning values are constant. No optimization possible.")
    rval = list(chosenTuning=getMHTuningList(mhOptimConfig$targets))
    return(rval)
  }
  
  # Copy targets to be modified as MH tunings are updated
  nextTargets = mhOptimConfig$targets
  
  # Only use startingParamValues on the first sample. Continue the chain after that.
  nextStartingValues = modCfg$startingParamValues
  
  sampleResults = vector('list', mhOptimConfig$optimSamples)
  bestInfo = NULL
  
  for (sampleIndex in 1:mhOptimConfig$optimSamples) {
    
    #
    sampleOptimCfg = mhOptimConfig
    sampleOptimCfg$targets = nextTargets
    
    # Silent after first sample
    runConfig = makeRunConfig(mhOptimConfig$optimSampleIter, verbose = (sampleIndex == 1))
    
    # Copy model config for this sample
    sampleModCfg = modCfg
    
    # Starting values to continue chain
    sampleModCfg$startingParamValues = nextStartingValues
    
    sampleModCfg$mhTuningOverrides = getMHTuningList(nextTargets, modCfg$mhTuningOverrides)
    
    # Run sample
    cat(paste0("\nRunning MH optimization sample ", sampleIndex, "/", sampleOptimCfg$optimSamples, "\n"))
    
    res = RPE_internal(data=data, modCfg=sampleModCfg, runConfig=runConfig)
    
    sampleResults[[sampleIndex]] = res
    
    
    # Update tuning values
    update = MHOptim_updateTuning(res, sampleOptimCfg, proportionComplete = sampleIndex / sampleOptimCfg$optimSamples)
    
    # Copy updated tuning values into nextTargets
    for (i in 1:nrow(nextTargets)) {
      if (nextTargets$parName[i] %in% names(update$updatedTuning)) {
        nextTargets$currentTuning[i] = update$updatedTuning[[ nextTargets$parName[i] ]]
      }
    }
    
    # Track best tuning values
    update$bestInfo$sample = sampleIndex
    bestInfo = rbind(bestInfo, update$bestInfo)
    
    # Copy last values for next starting values
    lastIteration = removeBurnIn(res, res$runConfig$iterations - 1)
    nextStartingValues = lastIteration$posteriors
    
  }
  
  
  bestInfo = bestInfo[ order(bestInfo$parName), ] # Reorder for users
  
  # Fill in the tuning types
  tunings = list(last = list(), best = list(), bestHalfMean = list())
  
  for (n in unique(bestInfo$parName)) {
    
    # last: The last used tuning values are the current tuning values
    tunings$last[[n]] = nextTargets$currentTuning[ nextTargets$parName == n ]
    
    # The best and bestHalfMean tunings depend on the tuning steps
    bt = bestInfo[ bestInfo$parName == n, ]
    
    bt = bt[ order(abs(bt$dif)), ] # Order by lowest abs dif from target
    
    # best: For each parameter individually, take the lowest absolute difference.
    tunings$best[[n]] = bt$tuning[1]
    
    # bestHalfMean: Median split on abs(dif) of acceptance rates. Take the mean tuning value of the top half.
    bestDifs = 1:(floor(nrow(bt) / 2))
    
    tunings$bestHalfMean[[n]] = mean(bt$tuning[ bestDifs ])
    
  }
  
  chosenTuning = tunings[[ mhOptimConfig$chosenTuning ]]
  
  combinedRes = combineResults(doIntegrityChecks = FALSE, resList = sampleResults)
  combinedRes$MH$config = mhOptimConfig
  combinedRes$MH$tuning = sampleResults[[length(sampleResults)]]$MH$tuning
  
  # TODO: Where to put this?
  for (n in names(chosenTuning)) {
    # Set chosenTuning for the combined results.
    combinedRes$config$mhTuningOverrides[[n]] = chosenTuning[[n]]
    
    # Replace the tuning values that were used for the last sample with chosenTuning so that continueSampling works.
    # The actual tunings used for each sample are in sampleResults.
    combinedRes$MH$tuning[[n]] = chosenTuning[[n]]
  }
  

  rval = list(
    chosenTuning = chosenTuning,
    tunings = tunings,
    optimSteps = bestInfo,
    combinedResults = combinedRes, # ready to continueSampling
    sampleResults = sampleResults
  )
  
  cat("\nMH optimization complete.\n")
  
  rval
}


