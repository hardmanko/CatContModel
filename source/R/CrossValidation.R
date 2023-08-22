
#' Fit Model for N-Fold Cross Validation Analysis
#' 
#' N-Fold Cross Validation is a low-assumptions procedure that is used to estimate model fit.
#' In N-Fold Cross Validation, the data are divided into N "folds", each fold leaving out `1/N` observations from the data set.
#' Model parameters are estimated for each fold independently.
#' This function subsamples the data to create the folds, then [crossValidationStats()] is used to calculate fit statistics based on the estimated parameters.
#' 
#' @details
#' For each fold, part of the data are left out. When estimating parameters for fold 2, for example, all data labeled "2" are excluded for that fold. When calculating fit statistics, how well the model fits the data labeled 2 is calculated. 
#' To summarize, model parameters are estimated with data labeled anything other than 2 then model fit (given the estimated parameters) is calculated with data labeled 2.
#' The fold numbers associated with each observation are in the CVFold column of the "labeledData" member of the returned list.
#' 
#' Details about making folds: With N folds, every observation (a study/response pair) will be labeled with a number from 1 to N. For each condition within a participant, the observations are sorted by study value. Beginning with the lowest study value, every Nth observation is given the same number. The starting N varies to help each fold to contain data that are evenly distributed across the range of study values.
#' 
#' 
#' 
#' @param data See [makeModelConfig()].
#' @param modCfg See [makeModelConfig()].
#' @param nFolds The number of sub-data-sets to make, each leaving out `1/nFolds` proportion of the data.
#' @param iterations See [runParameterEstimation()].
#' @param parallelConfig See [makeParallelConfig()]. The `numChains` argument is ignored and instead set equal to `nFolds`.
#' @param mhOptim See [runParameterEstimation()].
#' @param runConfig See [runParameterEstimation()].
#' 
#' @return An object with class "CCM_CrossValidation" that contains the `labeledData` and the `fits`.
#' 
#' @seealso [crossValidationStats()]
#' 
crossValidationFit = function(data, modCfg, nFolds, iterations, parallelConfig=NULL, mhOptim=TRUE, runConfig=NULL) {
  
  labeledData = CV_labelFoldData(data, nFolds)
  
  if (is.null(runConfig)) {
    runConfig = makeRunConfig(iterations)
  }
  
  fits = list()
  if (!is.null(parallelConfig) && parallelConfig$numNodes > 1) {
    
  	# Is this the right place to do this? It happens again in crossValidationFit_parallel
    #parallelConfig$numChains = nFolds
    
    # Fit in parallel
    parRes = crossValidationFit_parallel(labeledData, modCfg=modCfg, parallelConfig=parallelConfig, mhOptim=mhOptim, runConfig=runConfig)
    
    fits = parRes$chains
    parallelConfig = parRes$parallelConfig # Save parallelConfig for calculating stats
    
  } else {
    for (fold in 1:nFolds) {
      fits[[fold]] = CVF_runSingleFit(labeledData, thisFold=fold, modCfg=modCfg, runConfig=runConfig, mhOptim=mhOptim)
    }
  }
  
  # Save parallelConfig for crossValidationStats
  rval = list(labeledData = labeledData, fits=fits, parallelConfig=parallelConfig)
  class(rval) = c(class(rval), "CCM_CrossValidation") # or something
  rval
}



CVF_runSingleFit = function(labeledData, thisFold, modCfg, runConfig, mhOptim) {
  
  included = labeledData[ labeledData$CVFold != thisFold, ]
  excluded = labeledData[ labeledData$CVFold == thisFold, ]
  
  rval = runParameterEstimation(included, config=modCfg, iterations=NULL, runConfig=runConfig, mhOptim=mhOptim)
  
  rval$crossValidation = list(includedData = included, excludedData = excluded)
  
  rval
}





crossValidationFit_parallel = function(labeledData, modCfg, parallelConfig, mhOptim, runConfig) {
  
  nFolds = length(unique(labeledData$CVFold))

  # Override numChains with nFolds
  parallelConfig$numChains = nFolds
  
  jobs = list()
  for (fold in 1:nFolds) {
    
    jobs[[fold]] = makeParallelJob(labeledData = labeledData, thisFold = fold,
                        modCfg = modCfg, mhOptim = mhOptim, runConfig = runConfig)
    
  }
  
  CVF_parallelWrapper = function(job) {

    rval = CVF_runSingleFit(labeledData=job$labeledData, thisFold=job$thisFold, modCfg=job$modCfg, 
    												runConfig=job$runConfig, mhOptim=job$mhOptim)
    
    rval$runConfig$seed = job$seed # save seed in runConfig
    
    rval
  }
  
  clustRes = runParallelJobs(parallelConfig, jobs, jobFun=CVF_parallelWrapper)
  

  rval = list(chains = clustRes, parallelConfig = parallelConfig)
  # TODO: Class as Parallel? (which would be technically wrong, because it is CV in parallel)
  rval
}


#' Calculate Fit Statistics for Cross Validation
#' 
#' You must first use [crossValidationFit()] to estimate parameters, then this function to calculate fit statistics.
#' 
#' @details
#' Model fit is calculated by using the estimated model parameters to make predictions about the data that were left out when the model parameters were estimated. In the case of CatContModel, the cross validation statistics are calculating using the study/response pairs that were left out of the fold from which the model parameters were estimated. 
#' 
#' Likelihood: For each study/response pair, the likelihood for that observation is calculated using the estimated model parameters.
#' 
#' Response Error: For each study value, a distribution of likely responses is sampled from and the prediction error is calculated. Differences in prediction error between different models can be used for model selection.
#' 
#' @param cvFits Return value from [crossValidationFit()].
#' @param nIter Number of iterations to use from the total posterior distribution.
#' @param burnIn Number of burn-in iterations to remove (see [removeBurnIn()]).
#' @param doAgg If `TRUE`, the cross validation statistics are aggregated. If `FALSE`, raw statistics are returned.
#' @param inParallel If `TRUE`, calculations are done in parallel processes.
crossValidationStats = function(cvFits, nIter=10, burnIn=0, doAgg=FALSE, inParallel=TRUE) {
  
	# For left out data, for each observation, calculate:
	# 1a) Likelihood for the observation.
	# Sample N new responses from the posterior predictive:
	# 2a) Calculate raw error (circular mean of PP samples - observation)
	# 2b) MSE or similar.
	#
	# Lower MSE wins.
	# Lower log likelihood wins.
	
	if (!resultIsType(cvFits, "CrossValidation")) {
		stop("Can use only CrossValidation results objects.")
	}
	
  nFolds = length(cvFits$fits)
  
  allStats = NULL
  
  if (inParallel) {

    #parallelConfig = makeParallelConfig(nFolds)
    parallelConfig = cvFits$parallelConfig
    parallelConfig$numChains = nFolds
    
    CV_statsSingleFold_parallelWrapper = function(job) {
      CV_statsSingleFold(job$fit, nIter=job$nIter, burnIn=job$burnIn)
    }
    
    jobs = list()
    for (fold in 1:nFolds) {
      jobs[[fold]] = makeParallelJob(fit = cvFits$fits[[fold]], nIter = nIter, burnIn = burnIn)
    }
    
    statRes = runParallelJobs(parallelConfig, jobs, CV_statsSingleFold_parallelWrapper)
    
    for (fold in 1:nFolds) {
      allStats = rbind(allStats, statRes[[fold]])
    }
    
  } else {
    for (fold in 1:nFolds) {
      
      foldStats = CV_statsSingleFold(cvFits$fits[[fold]], nIter=nIter, burnIn=burnIn)
      
      allStats = rbind(allStats, foldStats)
    }
  }
  
  allStats$logLike = log(allStats$likelihood)
  allStats$ppError2 = allStats$ppError^2
  
  # Do this agg as well?
  #pciAgg = aggregate(likelihood ~ pnum * cond * iteration, allStats, mean)
  
  if (!doAgg) {
    return(allStats)
  }
  
  meanLikelihood = stats::aggregate(likelihood ~ CVFold, allStats, mean)

  # Collapse across iterations first
  llIterMean = stats::aggregate(logLike ~ CVFold * pnum * cond, allStats, mean)
  llSum = stats::aggregate(logLike ~ CVFold, llIterMean, sum)
  
  rawError = stats::aggregate(ppError ~ CVFold, allStats, mean)
  MSE = stats::aggregate(ppError2 ~ CVFold, allStats, mean)
  
  # Collect results
  byFold = meanLikelihood[ , "CVFold", drop=FALSE ]
  
  byFold$meanLikelihood = meanLikelihood$likelihood
  byFold$logLikeSum = llSum$logLike
  byFold$meanResponseError = rawError$ppError
  byFold$MSE = MSE$ppError2
  
  rval = byFold
  
  #if (includeRaw) {
  #  rval = list(byFold=byFold, rawStats = allStats)
  #}
  
  rval
}

CV_statsSingleFold = function(fit, nIter, burnIn) {
  
  if (burnIn > 0) {
    fit = removeBurnIn(fit, burnIn)
  }
  
  if (nIter > fit$runConfig$iterations) {
    nIter = fit$runConfig$iterations
  }
  
  allCalc = NULL

  excludedData = fit$crossValidation$excludedData
  
  for (pnum in unique(excludedData$pnum)) {
    for (cond in unique(excludedData$cond[excludedData$pnum == pnum])) {
      
      pcData = excludedData[excludedData$pnum == pnum & excludedData$cond == cond, ]
      
      # TODO
      iterSample = sample(1:fit$runConfig$iterations, nIter)
      
      for (iter in iterSample) {
        itParam = getSingleIterationParameters(fit, pnum=pnum, cond=cond, iteration=iter)
        
        
        # 1a) Likelihood for left out data (using what parameters? all iterations?)
        likeData = likelihood(itParam, pcData, 
                              modelVariant=fit$config$modelVariant, 
                              dataType=fit$config$dataType, 
                              responseRange=fit$config$responseRange,
                              minSD=fit$config$minSD)
        
        likeData$iteration = iter

        # 2) Sample posterior predictive response (PPR) for each left out observation. Compare left out response to PPR.
        
        ppSamp = sampleDataFromModel(pcData$study, itParam, 
                                         modelVariant = fit$config$modelVariant, 
                                         dataType = fit$config$dataType, 
                                         responseRange = fit$config$responseRange)
        
        likeData$ppResp = ppSamp$response
        
        allCalc = rbind(allCalc, likeData)
      }
      
    }
  }

  # TODO: For linear, circDist is approximately correct, but not actually correct.
  allCalc$ppError = circDist(allCalc$response, allCalc$ppResp, absDist = FALSE, degrees = TRUE)
  
  allCalc
}


CV_labelFoldData = function(data, nFolds, minCellCount=5) {
  
  foldStart = 1 # nFolds - 1?
  
  foldData = NULL
  
  for (pnum in unique(data$pnum)) {
    for (cond in unique(data$cond)) {
      
      thisData = data[ data$pnum == pnum & data$cond == cond, ]
      
      if (nrow(thisData) < nFolds) {
        stop(paste0("Not enough data for each fold to have an observation in a pnum * cond cell. pnum = ", pnum, ", cond = ", cond))
      }
      
      # Order by study value
      thisData = thisData[ order(thisData$study), ]
      
      # Take every nFolds observation
      for (i in 1:nrow(thisData)) {
      	thisData$CVFold[i] = (((foldStart - 1) + (i - 1)) %% nFolds) + 1
      }
      
      foldData = rbind(foldData, thisData)
      
      # Rotate starting fold so that each pnum x cond will start on a rotated fold number.
      # This helps each fold to get a more even distribution of study values.
      foldStart = (foldData$CVFold[nrow(foldData)] %% nFolds) + 1
    }
  }
  
  if (minCellCount > 0) {
    # Also a warning for small amounts of data per fold, maybe calculate later? What is number?
    agg = stats::aggregate(study ~ pnum * cond * CVFold, foldData, length)
    agg$ObsCount = agg$study
    
    smallCounts = agg$ObsCount < minCellCount
    if (any(smallCounts)) {
      warning("Small observation count for some cells of the design (a cell is pnum * cond * fold).")
    	print(agg[ smallCounts, ])
    }
  }
  
  foldData
}

