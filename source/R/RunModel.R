
#' Configure Parameter Estimation Run
#' 
#' Creates a configuration list that controls some aspects of parameter estimation.
#' 
#' @param iterations Iterations of the Gibbs sampler to run.
#' @param iterationsPerStatusUpdate Sets the number of iterations between updates to the progress bar.
#' @param verbose If `TRUE`, extra information is printed during the start of the run.
#' 
#' @return A configuration list.
#' 
#' @export
makeRunConfig = function(iterations, iterationsPerStatusUpdate=10, verbose=TRUE) {
  
  list(iterations=iterations,
       iterationsPerStatusUpdate = iterationsPerStatusUpdate,
       verbose = verbose,
       packageVersion = utils::packageVersion("CatContModel")
       )
  
}

#' Estimate Parameters of the Models
#' 
#' This function runs the Gibbs sampler for the given data and model variant.
#' In addition to this function documentation, you should read the manual and look at some of the examples 
#' (see the Examples section of the manual).
#'
#' @param data The data to fit the model to. See details.
#' @param config A list of configuration options for the model. See [`makeModelConfig`]. 
#' @param iterations Number of iterations of the Gibbs sampler to run (after MH optimization is complete). Each iteration is a sample from the posterior distribution of the parameters.
#' @param parallel Configuration for parallel parameter estimation. Either a list like returned by [`makeParallelConfig`] or an integer giving the number of chains to run in parallel.
#' @param mhOptim Any of `TRUE`, `FALSE`, `NULL`, or the result of [`makeMHConfig`]. See details.
#' @param runConfig See [`makeRunConfig`]. If `runConfig` is provided, `iterations` is ignored.
#' 
#' @details 
#' The data must be in a `data.frame` with 4 columns: `pnum`, `cond`, `study`, and `response`.
#' 
#' + `pnum`: Participant number or identifier. May be a string.
#' + `cond`: Experimental condition. May be a string.
#' + `study`: Studied stimulus value.
#' + `response`: Participant response.
#' 
#' If `config$dataType == "circular"`, then `study` and `response` should be degrees in the interval `[0, 360)`. 
#' If `config$dataType == "linear"`, then `study` and `response` can be in any units. 
#' For default priors and MH tuning values, it is best for linear data to have the same range as circular data 
#' (i.e. approximately 360 units long, but 100 to 1000 units should be fine).
#' 
#' 
#' If `mhOptim` is `TRUE` or the return value of [`makeMHConfig`], Metropolis-Hastings tuning values will be optimized at the start of parameter estimation. 
#' If `TRUE`, the defaults of [`makeMHConfig`] are used.
#' MH optimization iterations are kept, but can be discarded as burn-in iterations with [`removeBurnIn`]. 
#' If the model converges during MH optimization, there is no need to discard all optimization iterations as burn-in.
#' If using parallel estimation, optimization is done in each chain separately.
#' 
#' The `iterations` argument specifies the number of iterations to sample after MH optimization is complete. If using MH optimization, `iterations` can be 0.
#' 
#' 
#' @return A list containing the results of the Gibbs samplers:
#' \tabular{ll}{
#' 	`data` \tab The data.\cr
#'  `config` \tab The model configuration. \cr
#' 	`priors` \tab The priors that were used.\cr
#'	`posteriors` \tab A named list of posterior chains for each parameter.\cr
#'	`MH` \tab A `list` with information about the Metropolis-Hastings `tuning` values and `acceptance` rates for parameters with MH steps.\cr
#' 	`startingValues` \tab A named list of the starting values of the parameters. \cr
#' 	`constantValueOverrides` \tab A named list of constant parameter value overrides that were provided by the user. \cr
#' 	`pnums` \tab A vector of the participant numbers. \cr
#' 	`equalityConstraints` \tab Made based on `config$conditionEffects` and `config$factors`. A list mapping from the names of parameters to the name of the parameter that they obtain their value from, if any. \cr
#' 	`runConfig` \tab The run configuration. See `runConfig` argument.
#' }
#' 
#' @export
#' 
runParameterEstimation = function(data, config, iterations, parallel=NULL, mhOptim=FALSE, runConfig=NULL) {
  
  parallelConfig = NULL
  if (!is.null(parallel)) {
    if (is.list(parallel)) {
      parallelConfig = parallel
    } else if (is.numeric(parallel) && length(parallel) == 1 && parallel > 1) {
      parallelConfig = makeParallelConfig(parallel, outputFile = "RunProgress_")
    }
  }
  
  if (is.null(runConfig)) {
    runConfig = makeRunConfig(iterations)
  }
  
  rval = NULL
  
  if (!is.null(parallelConfig)) {
    rval = RPE_parallel(data=data, modCfg=config, runConfig=runConfig, parallelConfig=parallelConfig, mhOptim=mhOptim)
    
  } else {
    rval = RPE_mhOptim(data=data, modCfg=config, runConfig=runConfig, mhOptim=mhOptim)
  }
  
  rval
  
}

# Branches depending on whether mhOptim is enabled.
RPE_mhOptim = function(data, modCfg, runConfig, mhOptim) {
  
  mhOptimConfig = NULL
  if (is.logical(mhOptim) && mhOptim == TRUE) {
    mhOptimConfig = makeMHConfig()
  } else if (is.list(mhOptim)) {
    mhOptimConfig = mhOptim
  }
  
  doOptim = is.list(mhOptimConfig) && mhOptimConfig$optimSamples > 0 && mhOptimConfig$optimSampleIter > 0
  
  rval = NULL
  
  if (doOptim) {
    mho = optimizeMHTuning(data, modCfg, mhOptimConfig = mhOptimConfig)
    
    ## Treat iterations as after MHO
    if (runConfig$iterations > 0) {
      rval = continueSampling(mho$combinedResults, runConfig$iterations, combinedOnly = TRUE)
    } else {
      rval = mho$combinedResults
    }
    
  } else {
    rval = RPE_internal(data=data, modCfg=modCfg, runConfig=runConfig)
  }
  
  rval
}



# Final RPE function that calls into C++
RPE_internal = function(data, modCfg, runConfig) {
  
  # Data checks are in checkModelConfig
  modCfg = checkModelConfig(data, modCfg, immediateWarnings = TRUE, checkData = TRUE)
  
  # Equality constraints are not currently part of model config
  equalityConstraints = getConstrainedConditionEffects(modCfg)
  
  if (runConfig$verbose) {
    cat(paste0("Running ", runConfig$iterations, " iterations.\n"))
  }
  
  # Run estimation in C++
  results = CCM_CPP_runParameterEstimation(
    data = data,
    modelConfig = modCfg,
    runConfig = runConfig,
    equalityConstraints = equalityConstraints
  )
  
  results$data = data
  results$config = modCfg
  results$runConfig = runConfig
  
  results$runConfig$packageVersion = utils::packageVersion("CatContModel")
  
  results$equalityConstraints = equalityConstraints

  
  ####
  # Create list of actual starting values from posterior chains
  results$startingValues = list()
  
  for (n in names(results$posteriors)) {
    # Make it so that constant parameters have the same length posterior as the other kinds of parameters
    if (length(results$posteriors[[n]]) == 1) {
      results$posteriors[[n]] = rep(results$posteriors[[n]], results$runConfig$iterations + 1) # +1 because of the start value
    }
    
    # Store actual starting values rather than startingValueOverrides, which may be only a few of the parameters
    results$startingValues[[n]] = results$posteriors[[n]][1]
    
    results$posteriors[[n]] = results$posteriors[[n]][-1] # strip off the start value
  }
  
  # Set the class of results object to within-participants (WP)
  class(results) = c(class(results), "CCM_WP")
  
  results
}




#' Sample Additional Iterations of the Gibbs Sampler
#' 
#' After running the parameter estimation with [`runParameterEstimation`] and analyzing the results, 
#' you may decide that you want more iterations to be sampled. 
#' This function allows you to continue sampling iterations from where the parameter estimation left off.
#' 
#' @param results Results from [`runParameterEstimation`].
#' @param iterations The number of new iterations to sample (if `totalIterations == FALSE`).
#' @param totalIterations If `TRUE`, the `iterations` argument is treated as the total number of iterations that are desired. If the current number of samples is greater than or equal to `iterations`, this function returns immediately.
#' @param combinedOnly If `TRUE`, only the combined results (old combined with new) will be returned. If `FALSE`, the old, new, and combined results will be returned in a list.
#' 
#' 
#' @return A list with three elements: The `oldResults` (what you passed in), the `newResults` (the additional iterations), and the `combinedResults` (the `oldResults` and `newResults` combined using [`combineResults`]).
#' 
#' @family WP functions
#' 
#' @export
continueSampling = function(results, iterations, totalIterations=FALSE, combinedOnly=TRUE) {
  
  # TODO: Take parallelConfig as argument?
  
  # TODO: Allow character results? 
  #resultsFile = NULL
  #if (is.character(results)) {
    # Treat as a filename
  #  resultsFile = results
  #  if (!file.exists(resultsFile)) {
  #    stop(paste0("Cannot continueSampling: File not found: ", resultsFile))
  #  }
  #  results = readRDS(resultsFile)
  #}
  
  if (iterations <= 0) {
    warning("continueSampling can't be done if requested iterations <= 0. Returning the provided results.")
    return(results)
  }
  
  if (totalIterations) {
    completedIter = getCompletedIterations(results)
    if (completedIter >= iterations) {
      message("No need to continueSampling: Enough iterations have already been sampled. Returning the provided results.")
      return(results)
    }
    iterations = iterations - completedIter
  }
  
  if (resultIsType(results, "Parallel")) {
    res = continueSampling_parallel(results, iterations)
    return(res)
  } else if (resultIsType(results, "BP")) {
    stop("Cannot use continueSampling directly with BP results. Call continueSampling() for each element of res$groups.")
  }
  
  if (results$runConfig$iterations < 1) {
    stop("No existing iterations from which to continue.")
  }
  
  # The starting values are the last value in the previous posterior chain.
  startValues = list()
  for (parName in names(results$posteriors)) {
    startValues[[parName]] = results$posteriors[[parName]][results$runConfig$iterations]
  }
	
	# Copy the old config but with the new number of iterations.
	config = results$config
	config$startingParamValues = startValues
	
	runConfig = results$runConfig
	runConfig$iterations = iterations
	
	# TODO: Should actual values from last run be set as overrides? Or assume same defaults? 
	# Maybe this is an application for runConfig$packageVersion? If different version, do full overrides.
	# config$mhTuningOverrides = results$MH$tunings
	# config$priorOverrides = results$priors
	# No change: constantParamValues should be the same in config.
	# No change: startingParamValues already does change.
	
	newResults = RPE_internal(data=results$data, modCfg=config, runConfig=runConfig)

	rval = combineResults(results, newResults)
	
	# TODO: Auto-save over previous?
	#if (!is.null(resultsFile)) {
	#  saveRDS(rval, resultsFile)
	#}
	
	if (!combinedOnly) {
	  rval = list(oldResults=results, newResults=newResults, combinedResults=rval)
	}
	rval
}


#' Combine Multiple Posterior Parameter Chains
#' 
#' Can be used to combine multiple posterior chains that used the same data, model configuration, etc.
#' Can be used to combine parallel chains.
#' Note that you must remove burn-in iterations before merging results.
#' 
#' @param ... Multiple results objects from the [`runParameterEstimation`] function.
#' @param resList Instead of providing individual results objects as `...`, you can use `resList` to provide a list of results objects (numeric indices starting with 1).
#' @param doIntegrityChecks If `TRUE`, check that all of the results are comparable in terms of priors, constant parameter values, etc. The only time you should set this to `FALSE` is if you think there is a bug in the integrity checks.
#' @param compareExclude A vector of names of comparisons to NOT do. See `include` argument of [`compareResults`] for names.
#'
#' @return The combined chains in an updated results object.
#' 
#' @family WP functions
#'
#' @export
combineResults = function(..., resList=NULL, doIntegrityChecks=TRUE, compareExclude=NULL) {

	resultsList = list(...)
	if (is.null(resList)) {
	  # Do special for Parallel dots
		if (length(resultsList) == 1 && resultIsType(resultsList[[1]], "Parallel")) {
		  resultsList = resultsList[[1]]$chains
		  if (is.null(compareExclude)) {
		    compareExclude = "mhTuning"
		  }
		}
	} else {
	  resultsList = resList # TODO: Concatenate dots and resList?
	}
	
	if (length(resultsList) < 1) {
	  stop("No results provided.")
	} else if (length(resultsList) == 1) {
	  return(resultsList[[1]])
	}
	
	combined = resultsList[[1]]

	for (i in 2:length(resultsList)) {
		
		res = resultsList[[i]]
		
		# check that everything matches
		if (doIntegrityChecks) {
			compareResults(res, combined, exclude=compareExclude)
		}

		# Combine MH$acceptance rates
		combined$MH$acceptance$acceptanceRate = (combined$runConfig$iterations * combined$MH$acceptance$acceptanceRate) + 
		  (res$runConfig$iterations * res$MH$acceptance$acceptanceRate)
		combined$MH$acceptance$acceptanceRate = combined$MH$acceptance$acceptanceRate / (combined$runConfig$iterations + res$runConfig$iterations)
		
		
		# Tack on the posteriors
		for (n in names(combined$posteriors)) {
			combined$posteriors[[n]] = c(combined$posteriors[[n]], res$posteriors[[n]])
		}

		# Add the additional iterations
		combined$runConfig$iterations = combined$runConfig$iterations + res$runConfig$iterations
		
	}
	
	combined
}

#' Compare Two Results Objects for Configuration Consistency
#' 
#' The purpose of this function is to confirm that two sets of results have sampled from the same data and model.
#' The aspects to be compared are given by the `include` argument, excluding those aspects given by `exclude`.
#' This is a WP function only: No Parallel or BP objects allowed.
#' 
#' @param res1 The first results object from [`runParameterEstimation`].
#' @param res2 The second results object.
#' @param include A vector of comparisons to do.
#' @param exclude A vector of comparisons to NOT do. See names in `include`.
#' 
#' @return Nothing, but stops if the results don't match.
#' 
#' @family WP functions 
#' @export
compareResults = function(res1, res2, 
                          include = c("data", "modCfg", "priors", "constantParameters", "equalityConstraints", "mhTuning"), 
                          exclude=NULL)
{
	
  if (!(resultIsType(res1, "WP") && resultIsType(res2, "WP"))) {
    stop("Can only compare WP results.")
  }
  
  # Strip exclude from include
  include = include[ !(include %in% exclude) ]
  
  #errorFun = stop
  
  if ("data" %in% include) {
    if (any(res1$data != res2$data)) {
      stop("Mismatched data.")
    }
  }
	
	if ("modCfg" %in% include) {
	  
	  # These are simple items to compare
	  configTest = c("modelVariant", "dataType", "maxCategories", "minSD", 
	                 "catMuPriorApproximationPrecision", "factors", "calculateParticipantLikelihoods")
	  
	  if (res1$config$dataType == "linear") {
	    configTest = c(configTest, "responseRange", "catMuRange")
	  }
	  
	  for (n in configTest) {
      if (any(res1$config[[n]] != res2$config[[n]])) {
        stop(paste0("Mismatched config for setting ", n, "."))
      }
	  }
	  
	  # Check conditionEffects
		for (n in names(res1$config$conditionEffects)) {
			if (any(res1$config$conditionEffects[[n]] != res2$config$conditionEffects[[n]])) {
			  stop("Mismatched conditionEffects in config.")
			}
		}
	  
	  # Check constantParamValues outside of config
	  
	  # Don't compare starting param values (continueSampling).
	  
	  # Don't check priorOverrides because full priors are checked.
	  
	}
  
  if ("constantParameters" %in% include) {
    for (n in names(res1$config$constantParamValues)) {
      if (res1$config$constantParamValues[[n]] != res2$config$constantParamValues[[n]]) {
        stop(paste0("Mismatched constantParamValues for parameter ", n))
      }
    }
  }
  
  # Check all priors, not just the ones in config
  if ("priors" %in% include) {
    for (n in names(res1$priors)) {
      if (res1$priors[[n]] != res2$priors[[n]]) {
        stop(paste0("Mismatched priors for parameter ", n, "."))
      }
    }
  }
  
  # MH tuning values get warnings
  if ("mhTuning" %in% include) {
    mismatchedParams = c()
    
    for (n in names(res1$MH$tuning)) {
      if (res1$MH$tuning[[n]] != res2$MH$tuning[[n]]) {
        mismatchedParams = c(mismatchedParams, n)
      }
    }
    
    if (length(mismatchedParams) > 0) {
      warning(paste0("Mismatched MH tuning values for: ", paste(mismatchedParams, collapse=", "), "."))
    }
  }
  
  
  if ("equalityConstraints" %in% include) {
    if (length(res1$equalityConstraints) > 0 && length(res2$equalityConstraints) > 0) {
      if (any(sort(names(res1$equalityConstraints)) != sort(names(res2$equalityConstraints)))) {
        stop("Names of equality constraints do not match.")
      }
      for (n in names(res1$equalityConstraints)) {
        if (res1$equalityConstraints[[n]] != res2$equalityConstraints[[n]]) {
          stop(paste("Mismatched equality constraints for parameter ", n, ".", sep=""))
        }
      }
    } else if (length(res1$equalityConstraints) != length(res2$equalityConstraints)) {
      stop("The number of equality constraints do not match.")
    }
  }
	
}


