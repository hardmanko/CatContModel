

#' Verify Parameter Estimation Configuration Values
#' 
#' This function is used to verify the correctness of a configuration list and to add in default values for all missing elements. This function is used internally by, e.g., [`runParameterEstimation`].
#' 
#' @param config A configuration to be used as the `config` argument of [`runParameterEstimation`].
#' @param data The data that will be used as the `data` argument of [`runParameterEstimation`].
#' @param immediateWarnings If `TRUE`, warnings will be printed immediately. Regardless of the value, warnings will also be stored for later access with `warnings()`.
#' 
#' @return An updated configuration list, possibly with additional items added or existing items modified.
#' 
#' @export
#' 
#' @examples 
#' config = list(iterations = 1000, modelVariant = "betweenItem")
#' 
#' config$factors = data.frame(
#'   cond = c("a1", "a2", "a3", "b1", "b2", "b3"),
#'   letters = rep(c("a", "b"), each = 3),
#'   numbers = rep(c("1", "2", "3"), 2)
#' )
#' 
#' config$conditionEffects = list(pMem = c("letters", "numbers"),
#'   pContBetween = "all", #shortcut to listing them all
#'   contSD = "letters",
#'   pCatGuess = "numbers")
#' # All unspecified parameters will be set to "none".
#' 
#' # Make fake data
#' data = data.frame(pnum = 1, cond = config$factors$cond,
#'   study = 0, response = 0)
#' 
#' config = verifyConfigurationList(config, data)
verifyConfigurationList = function(config, data, immediateWarnings = FALSE) {
	
	usedConfigKeys = names(config)
	
	allAllowedConfigKeys = c("iterations", "modelVariant", "iterationsPerStatusUpdate", 
													 "cornerstoneConditionName", "maxCategories", "minSD", 
													 "catMuPriorApproximationPrecision",
													 "calculateParticipantLikelihoods", "conditionEffects",
													 "dataType", "responseRange", "catMuRange", "factors")
	
	disallowedConfigKeys = usedConfigKeys[ !(usedConfigKeys %in% allAllowedConfigKeys) ]
	
	if (length(disallowedConfigKeys) > 0) {
		msg = paste("The following configuration settings were used, but are not allowed (check their spelling): ", paste(disallowedConfigKeys, collapse=", "), sep="")
		stop(msg)
	}
	
	###############
	# iterations
	if (is.null(config$iterations)) {
		stop("config$iterations not set.")
	}
	
	##################################
	# iterationsPerStatusUpdate
	if (is.null(config$iterationsPerStatusUpdate)) {
		config$iterationsPerStatusUpdate = 10
		cat(paste("Note: config$iterationsPerStatusUpdate not set. Set to ", config$iterationsPerStatusUpdate, ".\n", sep=""))
	}
	
	####################
	# dataType
	
	if (is.null(config$dataType)) {
		config$dataType = "circular"
		cat(paste("Note: config$dataType not set. Set to ", config$dataType, ".\n", sep=""))
	} else if (!(config$dataType %in% c("circular", "linear"))) {
		stop("Invalid value for config$dataType. Choose from one of \"circular\" or \"linear\".")
	}
	
	if (config$dataType == "linear") {
		
		if (!is.null(config$responseRange) && length(config$responseRange) != 2) {
			stop("responseRange provided, but it is not valid. It must be a length 2 vector of numeric values.")
		}
		
		if (is.null(config$responseRange)) {
			config$responseRange = range(data$response)
			cat(paste("Note: config$responseRange not set. Set to (", paste(config$responseRange, collapse=", "), ").\n", sep=""))
		}

		if (config$responseRange[1] >= config$responseRange[2]) {
			stop("config$responseRange[1] >= config$responseRange[2]. The lower end of the response range must be lower!")
		}
		
		if (any(data$response < config$responseRange[1]) || any(data$response > config$responseRange[2])) {
			stop("Some responses are outside of responseRange. responseRange must contain all responses.")
		}
		
		if (is.null(config$catMuRange)) {
			config$catMuRange = config$responseRange #TODO: catMuRange is not currently documented anywhere. This is intentional, but you should still decide what to do with it.
		}
		
		if (config$catMuRange[1] >= config$catMuRange[2]) {
			stop("config$catMuRange[1] >= config$catMuRange[2]. The lower end of the catMu range must be lower!")
		}

	}
	
	########################
	# modelVariant
	possibleModelVariants = c("betweenAndWithin", "betweenItem", "withinItem", "ZL")
	visibleModelVariants = c("betweenItem", "withinItem", "ZL")
	if (is.null(config$modelVariant)) {
		stop(paste("config$modelVariant not set. Choose from one of: ", paste(visibleModelVariants, collapse = ", "), ".", sep="" ))
	}
	if (!(config$modelVariant %in% possibleModelVariants)) {
		stop(paste("Invalid model variant '", config$modelVariant, "' selected. Choose from one of: ", 
							 paste(visibleModelVariants, collapse = ", "), ".", sep="" ))
	}
	
	##################################
	# cornerstoneConditionName
	if (is.null(config$cornerstoneConditionName)) {
		#choose cornerstone condition based on the amount of data in the condition
		dataCounts = stats::aggregate(study ~ cond, data, length)
		config$cornerstoneConditionName = dataCounts$cond[which.max(dataCounts$study)]
		cat(paste("Note: config$cornerstoneConditionName not set. Set to ", config$cornerstoneConditionName, ".\n", sep=""))
	}
	if (!is.character(config$cornerstoneConditionName)) {
		config$cornerstoneConditionName = as.character(config$cornerstoneConditionName)
	}
	
	######################
	# maxCategories
	if (config$modelVariant == "ZL") {
		config$maxCategories = 0
		cat("Note: config$maxCategories has been set to 0 because you are using the ZL modelVariant.\n")
	}
	if (is.null(config$maxCategories)) {
		config$maxCategories = 16
		cat(paste("Note: config$maxCategories not set. Set to ", config$maxCategories, ".\n", sep=""))
	}
	
	####################
	# minSD
	if (is.null(config$minSD)) {
		config$minSD = 1
		cat(paste("Note: config$minSD not set. Set to ", config$minSD, " degree(s).\n", sep=""))
	}
	
	##########################################
	# catMuPriorApproximationPrecision
	if (is.null(config$catMuPriorApproximationPrecision)) {
		config$catMuPriorApproximationPrecision = 60
		cat(paste("Note: config$catMuPriorApproximationPrecision not set. Set to ", config$catMuPriorApproximationPrecision, " points at which the prior is evaluated.\n", sep=""))
	}
	
	###########################################
	# calculateParticipantLikelihoods
	if (is.null(config$calculateParticipantLikelihoods)) {
		config$calculateParticipantLikelihoods = FALSE
		cat(paste("Note: config$calculateParticipantLikelihoods not set. Set to ", config$calculateParticipantLikelihoods, ".\n", sep=""))
	}
	
	
	#######################
	# factors
	if (is.null(config$factors)) {
		#If no factors provided, assume one-factor design
		config$factors = data.frame(cond=unique(data$cond))
		config$factors$Factor = config$factors$cond
		cat("Note: config$factors not provided. A one-factor design is assumed.\n")
	} else {
		dataConds = unique(data$cond)
		allConds = union( dataConds, config$factors$cond )
		
		for (cond in allConds) {
			if (!(cond %in% config$factors$cond)) {
				stop(paste("The condition ", cond, " is in the data but not in config$factors.", sep=""))
			}
			if (!(cond %in% dataConds)) {
				stop(paste("The condition ", cond, " is in config$factors but not in the data.", sep=""))
			}
		}
	}
	for (n in names(config$factors)) {
		config$factors[ , n ] = as.character(config$factors[ , n ])
	}
	if (any(c("key", "group") %in% names(config$factors))) {
		stop('config$factors cannot have factors named either "key" or "group". These are reserved column names.')
	}
	for (n in names(config$factors)) {
		nameBad = grepl("[:.]+", n)
		levelsBad = any(grepl("[:.]+", config$factors[,n]))
		if (nameBad || levelsBad) {
			stop('Factor names and factor levels may not contain colon (":") or period (".").')
		}
	}
	
	######################
	# conditionEffects
	config$conditionEffects = verifyConditionEffects(config, immediateWarnings=immediateWarnings)
	
	config
}

# Checks the correctness of the condition effects in 
# config$conditionEffects, including possibly modifying them.
# This function is used by verifyConfigurationList
verifyConditionEffects = function(config, immediateWarnings = FALSE) {
	
	parametersWithPossibleConditionEffects = getAllParams(NULL, config$modelVariant, filter=TRUE)
	
	if (is.null(config$conditionEffects)) {
		
		pceToUse = getDefaultParametersWithConditionEffects(config$modelVariant)
		
		config$conditionEffects = list()
		for (param in pceToUse) {
			config$conditionEffects[[ param ]] = "all"
		}
		
		cat(paste("Note: config$conditionEffects not set. It was set to use all factors for ", paste(pceToUse, collapse = ", "), ".\n", sep=""))
		
	} else { #!is.null(config$conditionEffects)
		
		for (n in names(config$conditionEffects)) {
			if (!(n %in% parametersWithPossibleConditionEffects)) {
				config$conditionEffects[[n]] = NULL
				msg = paste0("In config$conditionEffects, parameter \"", n, "\" was included, but it is not a parameter that can have condition effects (possibly because it is not used by the current modelVariant). Its condition effects have been ignored.")
				if (immediateWarnings) {
					warning(msg, immediate. = TRUE)
				}
				warning(msg)
			}
		}
	}
	
	#set all unmentioned condition effects to "none"
	for (n in parametersWithPossibleConditionEffects) {
		if (!(n %in% names(config$conditionEffects))) {
			config$conditionEffects[[n]] = "none"
		}
	}
	
	#Double check that condition effect factor names are in factors
	factNames = getAllFactorNames(config$factors)
	for (n in names(config$conditionEffects)) {
		ce = config$conditionEffects[[n]]
		if (sum(ce %in% c("all", "none")) >= 2) {
			stop(paste0("config$conditionEffects$", n, " contains more than one instance of \"all\" or \"none\"."))
		}
		if (length(ce) > 1 || !all(ce %in% c("all", "none"))) {
			notIn = ce[ !(ce %in% factNames) ]
			if (length(notIn) > 0) {
				stop( paste0("config$conditionEffects$", n, " contains factor names not found in config$factors. The bad factor name(s): ", paste(notIn, collapse=", "), ".") )
			}
		}
	}
	
	config$conditionEffects
}


# This function is used by runParameterEstimation.
# It checks that config$conditionEffects are correct, given the 
# constant value overrides that are being used.
# It returns a, potentially modified, copy of config$conditionEffects.
checkConditionEffectsGivenConstantParameters = function(config, constantValueOverrides, immediateWarnings = FALSE) {
	
	for (param in names(config$conditionEffects)) {
		
		if (length(config$conditionEffects[[param]]) > 1 || config$conditionEffects[[param]] != "none") {
			thisCPN = paste(param, "_cond[", config$factors$cond, "]", sep="")
			
			inCVO = thisCPN %in% names(constantValueOverrides)
			
			if (all(inCVO)) {
				cat( paste0("Note: Parameter ", param, " has constant value overrides on all condition effect parameters. It has been noted to have no condition parameters.\n") )
				
				config$conditionEffects[[param]] = "none"
				
			} else if (any(inCVO)) {
				msg = paste0("Parameter ", param, " has constant value overrides on some condition effect parameters. It will have condition effects estimated, but some follow-up tests may not work correctly. After parameter estimation is complete, consider setting results$config$conditionEffects$", param, " to \"none\".")
				if (immediateWarnings) {
					warning( msg, immediate. = TRUE )
				}
				warning( msg )
			}
		}
		
	}
	
	config$conditionEffects
}

# This function is used by runParameterEstimation.
# It checks the starting value overrides for correctness.
# It returns a, potentially modified, copy of the starting value overrides.
checkStartingValueOverrides = function(config, svo) {
	
	for (n in names(svo)) {
		#Convert catMu starting values to radians
		if (config$dataType == "circular" && grepl("catMu", n, fixed=TRUE)) {
			svo[[n]] = CatContModel::d2r(svo[[n]])
		}
		
		if (grepl("catActive", n, fixed=TRUE)) {
			if (!(svo[[n]] %in% c(0, 1))) {
				stop("Starting values of catActive must be either 0 or 1.")
			}
		}
	}
	
	svo
}

# This function is used by runParameterEstimation.
# It checks the constant parameter value overrides for correctness.
# It returns a, potentially modified, copy of the starting value overrides.
checkConstantValueOverrides = function(config, cvo) {
	
	for (n in names(cvo)) {
		#Convert catMu constant values to radians
		if (config$dataType == "circular" && grepl("catMu", n, fixed=TRUE)) {
			cvo[[n]] = CatContModel::d2r(cvo[[n]])
		}
		
		if (grepl("catActive", n, fixed=TRUE)) {
			if (!(cvo[[n]] %in% c(0, 1))) {
				stop("Constant values of catActive must be either 0 or 1.")
			}
		}
	}
	
	cvo
}




#' Estimate Parameters of the Models
#' 
#' This function runs the Gibbs sampler for the selected delayed estimation model variant and data. In addition to this function documentation, you should read the manual and look at some of the examples (see the Examples section of the manual).
#'
#' @param config A list of configuration options. See details.
#' @param data The data to use in a data frame with 4 columns: `pnum`, `cond`, `study`, and `response`. `pnum` and `cond` may be strings. If `config$dataType == "circular"`, then `study` and `response` should be degrees in the interval [0, 360). If `config$dataType == "linear"`, then `study` and `response` can be in any units.
#' @param mhTuningOverrides A list of overrides of the default tuning parameters for the Metropolis-Hastings algorithm. See the Tuning Metropolis-Hastings Acceptance Rates section of the manual for more information.
#' @param priorOverrides A list of overrides of the default priors. See details and/or the Priors section of the manual.
#' @param startingValueOverrides A list of overrides of the default starting values. This is a list mapping from parameter name to starting values.
#' @param constantValueOverrides A list parameters giving overrides to the parameter values, which sets the parameters to constant values. This is a list mapping from parameter name to starting values. See [`setConstantParameterValue`] and [`setConstantCategoryParameters`].
#' 
#' 
#' @return A list containing the results of the Gibbs sampler, as follows:
#' \tabular{ll}{
#'  `config` \tab The configuration that was used. \cr
#'	`posteriors` \tab A named list of posterior distributions for each parameter.\cr
#' 	`mhAcceptance` \tab A `data.frame` containing information about the acceptance rates of parameters with Metropolis-Hastings steps.\cr
#' 	`mhTuning` \tab The Metropolis-Hastings tuning parameters that were used. \cr
#' 	`data` \tab The provided data.\cr
#' 	`priors` \tab The priors that were used.\cr
#' 	`startingValues` \tab A named list of the starting values of the parameters. \cr
#' 	`constantValueOverrides` \tab A named list of constant parameter value overrides that were provided by the user. \cr
#' 	`pnums` \tab A vector of the participant numbers. \cr
#' 	`equalityConstraints` \tab A list mapping from the names of parameters to the name of the parameter that they obtain their value from, if any.
#' }
#' 
#' @details
#' The elements of `config` include:
#' \tabular{ll}{
#' 	`iterations` \tab Required. The number of iterations of the Gibbs sampler to run. No default value.\cr
#' 	`modelVariant` \tab Required. Which variant of the model to use. Choose from "betweenItem", "withinItem", "betweenAndWithin", and "ZL". "betweenItem" and "withinItem" are the two model variants used by Hardman, Vergauwe, and Ricker.  "betweenAndWithin" is an unpublished model variant that has some problems. "ZL" is the Zhang and Luck (2008) model.\cr
#' 	`dataType` \tab The type of data you have, either `"circular"` or `"linear"`. Defaults to `"circular"`. \cr
#' 	`cornerstoneConditionName` \tab The name of the condition that will be used as the cornerstone condition. Defaults to the condition with the most observations in the data set. \cr
#' 	`maxCategories` \tab The maximum number of categories that each participant is allowed to have. Defaults to 16.\cr
#' 	`conditionEffects` \tab A list mapping from the name of a parameter to a character vector with the names of the factors that parameter will be allowed to vary by. To use all factors (or if you have only 1 factor), use the special value `"all"`. To prevent the parameter from varying by task condition, use the special value `"none"`. The default varies by model variant. Parameters that you do not include will be set to `"none"`. See the example. \cr
#' 	`factors` \tab A `data.frame` with a column named `cond` which contains the unique conditions in the data set. The other columns are factors, each containing factor levels corresponding to the conditions. If you have only 1 factor, you do not need to set this, but you may want to as the column names are used to identify the factors. See the example. \cr
#' 	`minSD` \tab The minimum standard deviation of the Von Mises or normal distributions. This affects the `contSD`, `catSD`, and `catSelectivity` parameters. Defaults to 1. \cr
#' 	`calculateParticipantLikelihoods` \tab Boolean. Whether (log) likelihoods should be calculated on each iteration for each participant. See [`calculateInappropriateFitStatistics`] for a poor reason to set this to `TRUE`.  Defaults to `FALSE`. \cr
#' 	`responseRange` \tab Only used if using the `"linear"` `dataType`, the possible range of response values as a length 2 vector, where the first value in the vector is the lower limit and the second the upper limit. By default, this is taken to be the range of the response values from the data set.
#' }
#' 
#' @section Prior Overrides:
#' Some things about the priors can be adjusted with the `priorOverrides` argument. The default priors and the prior distributions are specified in Hardman, Vergauwe, and Ricker (2017). Most participant level parameters have a hierarchical Normal prior with estimated mean, `mu`, and variance, `sigma2`. The prior on `mu` is Normal with fixed mean and variance. The prior on `sigma2` is Inverse Gamma with fixed alpha and beta. To set the priors, use the `priorOverrides` argument. For example, to set the priors on the distribution of `pMem`, you would set `pMem.mu.mu`, `pMem.mu.var`, `pMem.var.a`, and/or `pMem.var.b`. For all of the parameters, the priors are set in the latent space. See the Priors section of the manual for more information.
#' 
#' The condition effect parameters have Cauchy priors with constant location and scale. For example, for `pMem`, the location and scale priors are `pMem_cond.loc` and `pMem_cond.scale`. The locations all default to 0. The scales are different for different parameters. The locations should almost certaintly always be 0 unless you have a compelling argument to use a different value.
#' 
#' @export
#' 
runParameterEstimation = function(config, data, mhTuningOverrides=list(), 
																	priorOverrides=list(), startingValueOverrides=list(), 
																	constantValueOverrides=list()) 
{

	config = verifyConfigurationList(config, data, immediateWarnings = TRUE)

	startingValueOverrides = checkStartingValueOverrides(config, startingValueOverrides)
	
	constantValueOverrides = checkConstantValueOverrides(config, constantValueOverrides)
	
	# Double check that config$conditionEffects is reasonable given the constantValueOverrides
	config$conditionEffects = checkConditionEffectsGivenConstantParameters(config, constantValueOverrides, immediateWarnings = TRUE)
	
	equalityConstraints = getConstrainedConditionEffects(config)
	
	#Check the data
	requiredColumns = c("pnum", "cond", "study", "response")
	if (!all(requiredColumns %in% names(data))) {
		stop( paste("The required columns are not in the data set. These columns are ", paste(requiredColumns, collapse=", "), sep="") )
	}
	
	if (config$dataType == "circular" && (all(data$study < 10) || all(data$response < 10))) {
		msg = "Data appear to be in radians rather than degrees. The data should be in degrees."
		warning(msg, immediate. = TRUE)
		warning(msg)
	}


	results = CCM_CPP_runParameterEstimation(generalConfig = config, 
												data = data, 
												mhTuningOverrides = mhTuningOverrides, 
												priorOverrides = priorOverrides,
												startingValueOverrides = startingValueOverrides,
												constantValueOverrides = constantValueOverrides,
												equalityConstraints = equalityConstraints) 
	
	results$config = config
	results$data = data
	
	results$equalityConstraints = equalityConstraints

	#Convert catMu from radians to degrees
	if (results$config$dataType == "circular" && results$config$maxCategories > 0) {
		for (p in unique(results$data$pnum)) {
			for (i in 1:results$config$maxCategories) {
				cmName = paste("catMu[", p, ",", i-1, "]", sep="")
				results$posteriors[[ cmName ]] = CatContModel::r2d(results$posteriors[[ cmName ]])
			}
		}
	}
	
	
	results$startingValues = list()
	
	for (n in names(results$posteriors)) {
		#Make it so that constant parameters have the same length posterior as the other kinds of parameters
		if ( length(results$posteriors[[n]]) == 1) {
			results$posteriors[[n]] = rep(results$posteriors[[n]], config$iterations + 1) #+1 because of the start value
		}
		results$startingValues[[n]] = results$posteriors[[n]][1]
		
		results$posteriors[[n]] = results$posteriors[[n]][-1] #strip off the start value (mainly so that the number of iterations is correct)
	}
	
	class(results) = c(class(results), "CCM_WP")
	results
	
}

#' Sample Additional Iterations of the Gibbs Sampler
#' 
#' After running the parameter estimation with [`runParameterEstimation`] and analyzing the results, you may decide that you want more iterations to be sampled. This function allows you to continue sampling iterations from where the parameter estimation left off.
#' 
#' @param results The results from the [`runParameterEstimation`] function.
#' @param iterations The number of new iterations to sample.
#' 
#' @return A list with three elements: The `oldResults` (what you passed in), the `newResults` (the additional iterations), and the `combinedResults` (the `oldResults` and `newResults` merged together).
#' 
#' @family WP functions
#' 
#' @export
continueSampling = function(results, iterations) {
	lastIteration = removeBurnIn(results, results$config$iterations - 1)
	startValues = lastIteration$posteriors
	
	config = results$config
	config$iterations = iterations
	
	newResults = runParameterEstimation(config=config, data=results$data, 
																			mhTuningOverrides = results$mhTuning, priorOverrides = results$priors,
																			constantValueOverrides = results$constantValueOverrides,
																			startingValueOverrides=startValues)
	
	list(oldResults=results, newResults=newResults, combinedResults=mergeResults(results, newResults))
}


#' Merge the Results of Multiple Parameter Estimation Runs
#' 
#' Note that you must remove burn-in iterations before merging results.
#' 
#' @param ... Multiple results objects from the [`runParameterEstimation`] function.
#' @param doIntegrityChecks If `TRUE`, check that all of the results are comparable in terms of priors, Metropolis-Hastings tuning paramters, etc. The only time you should set this to FALSE is if you think there is a bug in the integrity checks.
#' @param rList Instead of providing individual results objects as `...`, you can use `rList` to provide a list of results objects.
#'
#' @return The merged results.
#' 
#' @family WP functions
#'
#' @export
mergeResults = function(..., doIntegrityChecks=TRUE, rList=NULL) {

	resultsList = list(...)
	if (!is.null(rList)) {
		resultsList = rList
	}

	merged = resultsList[[1]]

	for (i in 2:length(resultsList)) {
		
		res = resultsList[[i]]
		
		# check that everything matches
		if (doIntegrityChecks) {
			compareResults(res, merged, data=TRUE, mhTuning=TRUE, priors=TRUE, constantValueOverrides=TRUE, 
										 equalityConstraints=TRUE, config=TRUE, 
										 configIgnore = c("iterations", "iterationsPerStatusUpdate"))
		}

		# Combine mhAcceptance rates
		merged$mhAcceptance$acceptanceRate = (merged$config$iterations * merged$mhAcceptance$acceptanceRate) + (res$config$iterations * res$mhAcceptance$acceptanceRate)
		merged$mhAcceptance$acceptanceRate = merged$mhAcceptance$acceptanceRate / (merged$config$iterations + res$config$iterations)
		
		
		# Tack on the posteriors
		for (n in names(merged$posteriors)) {
			merged$posteriors[[n]] = c(merged$posteriors[[n]], res$posteriors[[n]])
		}

		# Add the additional iterations
		merged$config$iterations = merged$config$iterations + res$config$iterations
		
	}
	
	merged
}

#' Compare Two Results Objects for Configuration Consistency
#' 
#' By default, nothing is compared. The various comparisons must be enabled individually. This is a WP function only.
#' 
#' @param res1 The first results object from [`runParameterEstimation`].
#' @param res2 The second results object.
#' @param data Compare data sets.
#' @param constantValueOverrides Compare `constantValueOverrides`.
#' @param equalityConstraints Compare `equalityConstraints`.
#' @param priors Compare `priors`.
#' @param mhTuning Compare `mhTuning`.
#' @param config Compare `config`.
#' @param configIgnore Elements of the config list to not compare.
#' @param msgFun A function of one argument to pass messages to. Defaults to `stop`. Another good choice is `warning`.
#' 
#' @family WP functions
#' 
#' @export
compareResults = function(res1, res2, data = FALSE, constantValueOverrides = FALSE, equalityConstraints = FALSE, priors = FALSE, mhTuning = FALSE, config = FALSE, configIgnore = c("iterationsPerStatusUpdate"), msgFun = stop) {
	
	#TODO: work on this
	if (FALSE) {
		messages = NULL
		conditionalMessage = function(source, subsource, errMsg, errCond) {
			if (errCond) {
				msgFun(errMsg)
			}
			messages = rbind(messages, data.frame(source = source, subsource = subsource, errCond = errCond))
		}
		
		if (mhTuning) {
			for (n in names(res1$mhTuning)) {
				conditionalMessage("mhTuning", n, 
													 paste0("Mismatched MH tuning for parameter ", n, "."),
													 res1$mhTuning[[n]] != res2$mhTuning[[n]])
			}
		}
		
		if (priors) {
			for (n in names(res1$priors)) {
				conditionalMessage("priors", n, 
													 paste0("Mismatched priors for parameter ", n, "."),
													 res1$priors[[n]] != res2$priors[[n]])
			}
		}
	}
	
	
	if (mhTuning) {
		for (n in names(res1$mhTuning)) {
			if (res1$mhTuning[[n]] != res2$mhTuning[[n]]) {
				msgFun(paste("Mismatched MH tuning for parameter ", n, ".", sep=""))
			}
		}
	}
	
	if (priors) {
		for (n in names(res1$priors)) {
			if (res1$priors[[n]] != res2$priors[[n]]) {
				msgFun(paste("Mismatched priors for parameter ", n, ".", sep=""))
			}
		}
	}
	
	
	if (constantValueOverrides) {
		if (length(res1$constantValueOverrides) > 0 && length(res2$constantValueOverrides) > 0) {
			if (any(sort(names(res1$constantValueOverrides)) != sort(names(res2$constantValueOverrides)))) {
				msgFun("Names of constant value overrides do not match.")
			}
			for (n in names(res1$constantValueOverrides)) {
				if (res1$constantValueOverrides[[n]] != res2$constantValueOverrides[[n]]) {
					msgFun(paste("Mismatched constant value overrides for parameter ", n, ".", sep=""))
				}
			}
		} else if (length(res1$constantValueOverrides) != length(res2$constantValueOverrides)) {
			msgFun("The number of constant value overrides do not match.")
		}
	}
	
	if (equalityConstraints) {
		if (length(res1$equalityConstraints) > 0 && length(res2$equalityConstraints) > 0) {
			if (any(sort(names(res1$equalityConstraints)) != sort(names(res2$equalityConstraints)))) {
				msgFun("Names of equality constraints do not match.")
			}
			for (n in names(res1$equalityConstraints)) {
				if (res1$equalityConstraints[[n]] != res2$equalityConstraints[[n]]) {
					msgFun(paste("Mismatched equality constraints for parameter ", n, ".", sep=""))
				}
			}
		} else if (length(res1$equalityConstraints) != length(res2$equalityConstraints)) {
			msgFun("The number of equality constraints do not match.")
		}
	}
	
	if (config) {
		for (n in names(res1$config)) {
			if (!(n %in% c(configIgnore, "conditionEffects"))) {
				if (any(res1$config[[n]] != res2$config[[n]])) {
					msgFun(paste("Mismatched config for setting ", n, ".", sep=""))
				}
			}
		}
		if (!("conditionEffects" %in% configIgnore)) {
			for (n in names(res1$config$conditionEffects)) {
				if (any(res1$config$conditionEffects[[n]] != res2$config$conditionEffects[[n]])) {
					msgFun("Mismatched conditioEffects.")
				}
			}
		}
	}
	
	if (data) {
		if (any(res1$data != res2$data)) {
			msgFun("Mismatched data.")
		}
	}
	
}

#' Examine the Metropolis-Hastings Acceptance Rates
#' 
#' Most of the parameters in the model do not have conjugate priors, so the Metropolis-Hastings (MH) procedure is used to sample from their conditional posterior distribution. This function helps to examine the acceptance rate of the MH procedure for different parameter groups.
#' 
#' When a category is inactive, `catMu` is always accepted. Thus, when considering the acceptance rate of the `catMu` parameters, you should only consider the acceptances of active categories. This function gives two acceptance rates for `catMu`: Active only, which you should use, and the total, including inactive categories, which is basically irrelevant.
#' 
#' Note that the results of this function are not changed by removing burn-in iterations, except that the `catMu` acceptance rates are recalculated in a way that makes them meaningless once burn-in iterations have been removed. In short: Only use this function before burn-in iterations are removed.
#' 
#' @param results The results from the [`runParameterEstimation`] function.
#' 
#' @return A `data.frame` containing a summary of the MH acceptance rates for each group of parameters.
#' 
#' @export
examineMHAcceptance = function(results) {
	
	summarizeAcceptance = function(x) {
		qs = as.vector(stats::quantile(x, c(0, 0.025, 0.5, 0.975, 1)))
		qs = c(mean(x), qs)
		qs
	}
	
	tempMhAccept = results$mhAcceptance
	
	#separate out active catMu acceptances.
	for (pnum in results$pnums) {
		
		for (cat in 1:results$config$maxCategories) {
			
			cmn = paste("catMu[", pnum, ",", cat - 1, "]", sep="")
			can = paste("catActive[", pnum, ",", cat - 1, "]", sep="")
			
			if (cmn %in% names(results$constantValueOverrides) || !(cmn %in% names(results$posteriors))) {
				next #skip constant catMu and, if somehow the catMu is not in the posteriors, skip it as well
			}
			
			totalPossibleAcceptances = results$config$iterations
			
			nActive = sum(results$posteriors[[ can ]])
			nInactive = totalPossibleAcceptances - nActive
			
			totalAcceptances = results$mhAcceptance[ results$mhAcceptance$name == cmn, ]$acceptanceRate * totalPossibleAcceptances
			activeAcceptances = totalAcceptances - nInactive
			activeAcceptanceRate = activeAcceptances / nActive
			if (is.nan(activeAcceptanceRate) || is.na(activeAcceptanceRate) || activeAcceptanceRate < 0) {
				activeAcceptanceRate = 0
			}
			
			temp = data.frame(name=paste(cmn, "(Active)"), group="catMu (Active only)", acceptanceRate = activeAcceptanceRate)
			
			tempMhAccept = rbind(tempMhAccept, temp)

		}
	}
	
	
	summ = stats::aggregate(acceptanceRate ~ group, tempMhAccept, summarizeAcceptance)
	summ = do.call(data.frame, summ)
	
	names(summ) = c("paramGroup", "mean", "min", "2.5%", "median", "97.5%", "max")
	summ$paramGroup = as.character(summ$paramGroup)
	
	summ$paramGroup[ summ$paramGroup == "catMu" ] = "catMu (Total, ignore)"
	
	summ = summ[ order(summ$paramGroup), ]
	
	summ
}

