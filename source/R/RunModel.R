require(Rcpp)


getDefaultParametersWithConditionEffects = function(modelVariant) {
	pce = c("pMem", "contSD") #ZL and all other models
	
	if (modelVariant == "betweenAndWithin") {
		pce = c(pce, "pBetween", "pContBetween", "pContWithin")
		
	} else if (modelVariant == "betweenItem") {
		pce = c(pce, "pContBetween")
		
	} else if (modelVariant == "withinItem") {
		pce = c(pce, "pContWithin")
	}
	
	pce
}


#' Verify Parameter Estimation Configuration Values
#' 
#' @param config A configuration to be used as the \code{config} argument of \link{runParameterEstimation}.
#' @param data The data that will be used as the \code{data} argument of \link{runParameterEstimation}.
#' 
#' @return An updated configuration list.
#' 
#' @export
verifyConfigurationList = function(config, data) {
	
	allAllowedConfigKeys = c("iterations", "modelVariant", "iterationsPerStatusUpdate", 
													 "cornerstoneConditionName", "maxCategories", "minSD", 
													 "calculateParticipantLikelihoods", "conditionEffects",
													 "dataType", "responseRange", "catMuRange", "factors", "factorNames")
	disallowedConfigKeys = names(config)[ !(names(config) %in% allAllowedConfigKeys) ]
	if (length(disallowedConfigKeys) > 0) {
		msg = paste("The following configuration settings were used, but are not allowed (check their spelling): ", paste(disallowedConfigKeys, collapse=", "), sep="")
		stop(msg)
	}

	
	if (is.null(config$iterations)) {
		stop("config$iterations not set.")
	}
	
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
	

	possibleModelVariants = c("betweenAndWithin", "betweenItem", "withinItem", "ZL")
	visibleModelVariants = c("betweenItem", "withinItem", "ZL")
	if (is.null(config$modelVariant)) {
		stop(paste("config$modelVariant not set. Choose from one of: ", paste(visibleModelVariants, collapse = ", "), ".", sep="" ))
	}
	if (!(config$modelVariant %in% possibleModelVariants)) {
		stop(paste("Invalid model variant '", config$modelVariant, "' selected. Choose from one of: ", 
							 paste(visibleModelVariants, collapse = ", "), ".", sep="" ))
	}
	
	
	if (is.null(config$iterationsPerStatusUpdate)) {
		config$iterationsPerStatusUpdate = 10
		cat(paste("Note: config$iterationsPerStatusUpdate not set. Set to ", config$iterationsPerStatusUpdate, ".\n", sep=""))
	}
	
	
	if (is.null(config$cornerstoneConditionName)) {
		#choose cornerstone condition based on the amount of data in the condition
		dataCounts = stats::aggregate(study ~ cond, data, length)
		config$cornerstoneConditionName = dataCounts$cond[which.max(dataCounts$study)]
		cat(paste("Note: config$cornerstoneConditionName not set. Set to ", config$cornerstoneConditionName, ".\n", sep=""))
	}
	if (!is.character(config$cornerstoneConditionName)) {
		config$cornerstoneConditionName = as.character(config$cornerstoneConditionName)
	}
	
	
	if (is.null(config$maxCategories)) {
		config$maxCategories = 16
		if (config$modelVariant == "ZL") {
			config$maxCategories = 0
		}
		cat(paste("Note: config$maxCategories not set. Set to ", config$maxCategories, ".\n", sep=""))
	}
	

	
	if (is.null(config$minSD)) {
		config$minSD = 1
		cat(paste("Note: config$minSD not set. Set to ", config$minSD, " degree.\n", sep=""))
	}
	
	if (is.null(config$calculateParticipantLikelihoods)) {
		config$calculateParticipantLikelihoods = FALSE
		cat(paste("Note: config$calculateParticipantLikelihoods not set. Set to ", config$calculateParticipantLikelihoods, ".\n", sep=""))
	}
	
	
	parametersWithPossibleConditionEffects = c( getProbParams(NULL), getSdParams(NULL) )

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
				msg = paste("In config$conditionEffects, ", n, " was included, but it is not a parameter that can have condition effects.", sep="")
				cat( paste( msg, "\n", sep="") )
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
	
	if (is.null(config$factorNames)) {
		config$factorNames = guessFactorNames(config$factors)
	}
	
	config
}


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
#' This function runs the Gibbs sampler for the selected delayed estimation model.
#'
#' @param config A list of configuration options. See details.
#' @param data The data to use in a data frame with 4 columns: \code{pnum}, \code{cond}, \code{study}, and \code{response}. \code{pnum} and \code{cond} may be strings. \code{study} and \code{response} should be degrees in the interval [0, 360).
#' @param mhTuningOverrides A list of overrides of the default tuning parameters for the Metropolis-Hastings algorithm. For all parameters with Metropolis-Hastings steps, a normal candidate distribution centered on the current value is used. These tuning parameters are the standard deviations of the candidate distributions.
#' @param priorOverrides A list of overrides of the default priors. See details.
#' @param startingValueOverrides A list of overrides of the default starting values.
#' @param constantValueOverrides A list parameters giving overrides to the parameter values, which sets the parameters to constant values.
#' 
#' 
#' @return An list containing the results of the Gibbs sampler, as follows:
#' \tabular{ll}{
#'  \code{config} \tab The primary configuration that was used. \cr
#'	\code{posteriors} \tab A named list of posterior distributions for each parameter.\cr
#' 	\code{mhAcceptance} \tab A data frame containing information about the acceptance rates of parameters with Metropolis-Hastings steps.\cr
#' 	\code{mhTuning} \tab The Metropolis-Hastings tuning parameters that were used. \cr
#' 	\code{data} \tab The provided data.\cr
#' 	\code{priors} \tab The priors that were used.\cr
#' 	\code{startingValues} \tab A named list of the starting values of the parameters. \cr
#' 	\code{constantValueOverrides} \tab A named list of constant parameter value overrides (that were provided by the user). \cr
#' 	\code{pnums} \tab A vector of the participant numbers. \cr
#' 	\code{equalityConstraints} \tab A list mapping from the names of parameters to the name of the parameter that they obtain their value from, if any.
#' }
#' 
#' @details
#' The elements of \code{config} include:
#' \tabular{ll}{
#' 	\code{iterations} \tab Required. The number of iterations of the Gibbs sampler to run. No default value.\cr
#' 	\code{modelVariant} \tab Required. Which variant of the model to use. Choose from "betweenItem", "withinItem", and "ZL". "betweenItem" and "withinItem" are the two model variants used by Hardman, Vergauwe, and Ricker. "ZL" is the Zhang and Luck (2008) model. Defaults to "betweenItem".\cr
#' 	\code{dataType} \tab The type of data you have, either \code{"circular"} or \code{"linear"}. Defaults to \code{"circular"}. \cr
#' 	\code{cornerstoneConditionName} \tab The name of the condition that will be used as the cornerstone condition. Defaults to the condition with the most observations in the data set.\cr
#' 	\code{maxCategories} \tab The maximum number of categories that each participant is allowed to have. Defaults to 16.\cr
#' 	\code{conditionEffects} \tab Replaces parametersWithConditionEffects. A list mapping from the name of a parameter to a character vector with the names of the factors that parameter will be allowed to vary by. To use all factors (or if you have only 1 factor), use the special value \code{"all"}. To prevent the parameter from varying by task condition, use the special value \code{"none"}. The default varies by model variant. Parameters that you do not include will be set to \code{"none"}. See the example. \cr
#' 	\code{factors} \tab A \code{data.frame} with a column named "cond" which contains the unique conditions in the data set. The other columns are factors and factor levels corresponding to the conditions. If you have only 1 factor, you do not need to set this, but you may want to as the column names are used to identify the factors. See the example. \cr
#' 	\code{minSD} \tab The minimum standard deviation of the Von Mises or normal distributions. This affects the contSD, catSD, and catSelectivity parameters. Defaults to 1. \cr
#' 	\code{calculateParticipantLikelihoods} \tab Whether (log) likelihoods should be calculated on each iteration for each participant. See \code{\link{calculateInappropriateFitStatistics}}. \cr
#' 	\code{responseRange} \tab If using the "linear" dataType, the possible range of response values as a length 2 vector, where the first value in the vector is the lower limit. By default, this is inferred from the data.
#' }
#' 
#' Some things about the priors can be adjusted with the \code{priorOverrides} argument. The priors are specified in Hardman, Vergauwe, and Ricker. Most participant level parameters have a hierarchical Normal prior with estimated mean (mu) and variance (sigma2). The prior on mu is Normal with fixed mean and variance. The prior on sigma2 is Inverse Gamma with fixed alpha and beta. To set the priors, use the \code{priorOverrides} argument. For example, to set the priors on the distribution of \code{pMem}, you would set \code{pMem.mu.mu}, \code{pMem.mu.var}, \code{pMem.var.a}, and/or \code{pMem.var.b}. pMem.mu.var sets the prior variance on mu of the pMem parameters. For all of the parameters, the priors are set in the latent space.
#' 
#' The condition effect parameters have fixed Cauchy priors with a location and scale. For example, for pMem, the location and scale priors are \code{pMem_cond.loc} and \code{pMem_cond.scale}. The locations all default to 0. The scales are different for different parameters.
#' 
#' @export
#' 
#' @examples 
#' config = list(iterations = 1000, modelVariant = "betweenItem")
#' 
#' config$factors = data.frame(cond = c("a1", "a2", "a3", "b1", "b2", "b3"),
#' letters = rep(c("a", "b"), each = 3),
#' numbers = rep(c("1", "2", "3"), 2))
#' 
#' config$conditionEffects = list(pMem = c("letters", "numbers"),
#' pContBetween = "all", #shortcut to listing them all
#' contSD = "letters",
#' pCatGuess = "numbers")
#' #all not stated will be set to "none".
#' 
runParameterEstimation = function(config, data, mhTuningOverrides=list(), 
																	priorOverrides=list(), startingValueOverrides=list(), 
																	constantValueOverrides=list()) 
{

	config = verifyConfigurationList(config, data)
	

	startingValueOverrides = checkStartingValueOverrides(config, startingValueOverrides)
	
	constantValueOverrides = checkConstantValueOverrides(config, constantValueOverrides)
	
	# Double check that config$conditionEffects is reasonable given the constantValueOverrides
	# TODO: This should probably go into a function
	for (param in names(config$conditionEffects)) {
		
		if (config$conditionEffects[[param]] != "none") {
			thisCPN = paste(param, "_cond[", config$factors$cond, "]", sep="")
			
			inCVO = thisCPN %in% names(constantValueOverrides)
			
			if (all(inCVO)) {
				cat( paste0("Note: Parameter ", param, " has constant value overrides on all condition effect parameters. It has been noted to have no condition parameters.") )
				
				config$conditionEffects[[param]] = "none"
				
			} else if (any(inCVO)) {
				warning( paste0("Parameter ", param, " has constant value overrides on some condition effect parameters. It will have condition effects estimated, but some tests may not work correctly. After parameter estimation is complete, consider setting results$config$conditionEffects$", param, " to \"none\".") )
			}
		}

	}

	
	equalityConstraints = getConstrainedConditionEffects(config)
	
	#Check the data
	requiredColumns = c("pnum", "cond", "study", "response")
	if (!all(requiredColumns %in% names(data))) {
		stop( paste("The required columns are not in the data set. These columns are ", paste(requiredColumns, collapse=", "), sep="") )
	}
	
	if (config$dataType == "circular" && (all(data$study < 10) || all(data$response < 10))) {
		cat("Warning: Data appear to be in radians rather than degrees. The data should be in degrees.\n")
		warning("Data appear to be in radians rather than degrees. The data should be in degrees.")
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
	
	
	results
	
}

#' Sample Additional Iterations of the Gibbs Sampler
#' 
#' After running the parameter estimation with \code{\link{runParameterEstimation}} and analyzing the results, you may decide that you want more iterations to be sampled. This function allows you to continue sampling iterations from where the parameter estimation left off.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param iterations The number of new iterations to sample.
#' 
#' @return A list with three elements: The oldResults (what you passed in), the newResults (the additional iterations), and the combinedResults (the oldResults and newResults merged together).
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
#' @param ... Multiple results objects from the \code{\link{runParameterEstimation}} function.
#' @param doIntegrityChecks If TRUE, check that all of the results are comparable in terms of priors, Metropolis-Hastings tuning paramters, etc. The only time you should set this to FALSE is if you think there is a bug in the integrity checks.
#' @param rList Instead of providing individual results objects as \code{...}, you can use \code{rList} to provide a list of results objects.
#'
#' @return The merged results.
#'
#' @export
mergeResults = function(..., doIntegrityChecks=TRUE, rList=NULL) {

	resultsList = list(...)
	if (!is.null(rList)) {
		resultsList = rList
	}

	rval = resultsList[[1]]

	for (i in 2:length(resultsList)) {
		
		res = resultsList[[i]]
		
		# check that everything matches
		if (doIntegrityChecks) {
			
			for (n in names(res$mhTuning)) {
				if (res$mhTuning[[n]] != rval$mhTuning[[n]]) {
					stop(paste("Mismatched MH tuning for parameter ", n, ".", sep=""))
				}
			}

			
			for (n in names(res$priors)) {
				if (res$priors[[n]] != rval$priors[[n]]) {
					stop(paste("Mismatched priors for parameter ", n, ".", sep=""))
				}
			}
			
			
			if (length(res$constantValueOverrides) > 0 && length(rval$constantValueOverrides) > 0) {
				if (any(sort(names(res$constantValueOverrides)) != sort(names(rval$constantValueOverrides)))) {
					stop("Names of constant value overrides do not match.")
				}
				for (n in names(res$constantValueOverrides)) {
					if (res$constantValueOverrides[[n]] != rval$constantValueOverrides[[n]]) {
						stop(paste("Mismatched constant value overrides for parameter ", n, ".", sep=""))
					}
				}
			} else if (length(res$constantValueOverrides) != length(rval$constantValueOverrides)) {
				stop("The number of constant value overrides do not match.")
			}
			
			
			if (length(res$equalityConstraints) > 0 && length(rval$equalityConstraints) > 0) {
				if (any(sort(names(res$equalityConstraints)) != sort(names(rval$equalityConstraints)))) {
					stop("Names of equality constraints do not match.")
				}
				for (n in names(res$equalityConstraints)) {
					if (res$equalityConstraints[[n]] != rval$equalityConstraints[[n]]) {
						stop(paste("Mismatched equality constraints for parameter ", n, ".", sep=""))
					}
				}
			} else if (length(res$equalityConstraints) != length(rval$equalityConstraints)) {
				stop("The number of equality constraints do not match.")
			}
			
			
			for (n in names(res$config)) {
				if (!(n %in% c("iterations", "iterationsPerStatusUpdate", "conditionEffects"))) {
					if (any(res$config[[n]] != rval$config[[n]])) {
						stop(paste("Mismatched config for setting ", n, ".", sep=""))
					}
				}
			}
			
			
			if (any(res$data != rval$data)) {
				stop("Mismatched data.")
			}
		}
		

		# Combine mhAcceptance rates
		rval$mhAcceptance$acceptanceRate = (rval$config$iterations * rval$mhAcceptance$acceptanceRate) + (res$config$iterations * res$mhAcceptance$acceptanceRate)
		rval$mhAcceptance$acceptanceRate = rval$mhAcceptance$acceptanceRate / (rval$config$iterations + res$config$iterations)
		
		
		# Tack on the posteriors
		for (n in names(rval$posteriors)) {
			rval$posteriors[[n]] = c(rval$posteriors[[n]], res$posteriors[[n]])
		}

		# Add the additional iterations
		rval$config$iterations = rval$config$iterations + res$config$iterations
		
	}
	
	rval
}

#' Examine the Metropolis-Hastings Acceptance Rates
#' 
#' Most of the parameters in the model do not have conjugate priors, so the Metropolis-Hastings algorithm is used to sample from their conditional distribution. This function helps to examine the acceptance rate of the MH algorithm for different parameter groups.
#' 
#' When a category is inactive, catMu is always accepted. Thus. when considering the acceptance rate of the catMu parameters, you should only consider the acceptances of active categories. This function gives two acceptance rates for catMu: Active only, which you should use, and the total, including inactive categories, which is basically irrelevant.
#' 
#' Note that the results of this function are not changed by removing burn-in iterations, except that the catMu acceptance rates are recalculated in a way that makes them meaningless. In short: Only use this function before burn-in iterations are removed.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' 
#' @return A data frame containing a summary of the MH acceptance rates for each group of parameters.
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



#' Convert Participant Parameters to Matrices or Arrays
#' 
#' For standard parameters (for which each participant has only 1 parameter), the result is a matrix where each row of the matrix is an iteration and each column of the matrix is a participant.
#' 
#' For the category parameters (catMu and catActive), the result is an array where the 3 dimensions are participant, category, and iteration, in that order.
#' 
#' The column index for a pnum can be gotten with \code{which(results$pnums == pnum)}.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param A list of parameter names. If NULL, the default, all parameters are done.
#' 
#' @return A named list of matrices and/or arrays.
#' 
#' @export
convertPosteriorsToMatrices = function(results, param=NULL) {
	
	matrixParams = c("pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess", 
									 "contSD", "catSelectivity", "catSD")
	arrayParams = c("catActive", "catMu")
	
	if (!is.null(param)) {
		matrixParams = matrixParams[ matrixParams %in% param ]
		arrayParams = arrayParams[ arrayParams %in% param ]
	}

	post = list()
	for (mp in matrixParams) {
		post[[mp]] = matrix(0, nrow=results$config$iterations, ncol=length(results$pnums))
		colnames(post[[mp]]) = results$pnums
	}
	for (ap in arrayParams) {
		post[[ap]] = array(0, dim=c(length(results$pnums), results$config$maxCategories, results$config$iterations))
		dimnames(post[[ap]]) = list(results$pnums)
	}
	
	for (pInd in 1:length(results$pnums)) {
		
		pnum = results$pnums[pInd]
		
		istr = paste("[", pnum, "]", sep="")
		
		for (mp in matrixParams) {
			post[[mp]][,pInd] = results$posteriors[[paste(mp, istr, sep="")]]
		}
		
		if (results$config$maxCategories > 0) {
			for (cat in 1:results$config$maxCategories) {
				
				istr = paste("[", pnum, ",", cat - 1, "]", sep="")
				
				for (ap in arrayParams) {
					
					pname = paste(ap, istr, sep="")
					
					x = results$posteriors[[pname]]
					
					post[[ap]][pInd,cat,] = x
					
				}
				
			}
		}
	}
	
	post
}

#' Convert Posterior Distributions to a Single Matrix
#' 
#' This converts raw posteriors into a single matrix. This matrix can then be used with the boa or coda packages for assessing convergence. Some of the convergence diagnostics require you to make separate matrices for separate runs of the Gibbs sampler.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param stripConstantParameters Remove all parameters with a constant value. Constant parameters cannot converge.
#' @param stripCatActive Remove all of the cat active parameters. It is difficult to assess convergence for indicator parameters that are either 0 or 1.
#' @param stripCatMu Remove all of the category mean parameters. It is difficult to assess convergence for these parameters because they have multi-modal posterior distributions.
#' 
#' @return A matrix containing all of the posterior distributions for the selected parameters. Each column is one parameter.
#' 
#' @export
convertPosteriorsToMatrix = function(results, stripConstantParameters=TRUE, stripCatActive=TRUE, stripCatMu=TRUE) {

	rawPost = results$posteriors
		
	iterations = 0
	constantPar = NULL
	catActivePar = NULL
	catMuPar = NULL
	participantLLPar = NULL
	
	for (n in names(rawPost)) {
		iterations = max(c(iterations, length(rawPost[[n]])))
		
		if (all(rawPost[[n]] == rawPost[[n]][1])) {
			constantPar = c(constantPar, n)
		}
		
		if (grepl("catActive", n, fixed=TRUE)) {
			catActivePar = c(catActivePar, n)
		}
		
		if (grepl("catMu", n, fixed=TRUE)) {
			catMuPar = c(catMuPar, n)
		}
		
		if (grepl("participantLL", n, fixed=TRUE)) {
			participantLLPar = c(participantLLPar, n)
		}
		
	}
	
	excludedPar = participantLLPar
	if (stripConstantParameters) {
		excludedPar = c(excludedPar, constantPar)
	}
	if (stripCatActive) {
		excludedPar = c(excludedPar, catActivePar)
	}
	if (stripCatMu) {
		excludedPar = c(excludedPar, catMuPar)
	}
	
	usedParam = names(rawPost)
	usedParam = usedParam[ !(usedParam %in% excludedPar) ]
	
	np = length(usedParam)
	
	m = matrix(0, nrow=iterations, ncol=np)
	colnames(m) = usedParam
	
	for (i in 1:np) {
		
		temp = rawPost[[usedParam[i]]]
		
		if (length(temp) == 1) {
			temp = rep(temp, iterations)
		} else if (length(temp) != iterations) {
			warning(paste("Data vector for parameter ", usedParam[i], " of incorrect length.", sep=""))
		}
		
		m[,i] = temp
		
	}
	
	m
}

#' Remove Burn-In Iterations
#' 
#' Due to the large number of parameters in the model and the fact that many of the paramters have Metropolis-Hastings updating steps, convergence can be slow. This function removes burn-in iterations from the beginning of the chains so that only converged distributions can be analyzed.
#'  
#' Note that the Metropolis-Hastings acceptance rate given by \link{examineMHAcceptance} is unaffected by removing burn-in iterations.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param burnIn Integer. Number of burn-in iterations to remove or a vector with length greater than 1 giving the indices of the iterations to remove.
#' 
#' @return A new results object with burnIn iterations removed.
#' 
#' @export
removeBurnIn = function(results, burnIn) {
	
	if (length(burnIn) == 0) {
		warning("No burn-in iterations removed because length(burnIn) == 0.")
		return(results)
	}
	
	if (length(burnIn) == 1) {
		burnIn = 1:burnIn
	}
	
	if (burnIn[1] < 1) {
		stop("The lowest burn-in iteration is too low (less than 1).")
	}
	
	if (burnIn[length(burnIn)] > results$config$iterations) {
		stop("The highest burn-in iteration is too high.")
	}
	
	
	keep = 1:results$config$iterations
	keep = keep[-burnIn]
	
	
	for (n in names(results$posteriors)) {
		results$posteriors[[n]] = results$posteriors[[n]][keep]
	}
	
	results$config$iterations = length(keep)
	
	results
}


#' Set Parameters to Constant Values
#' 
#' This function helps with setting parameters to constant values. It returns a
#' list mapping from parameter name to parameter value. I sets all participant
#' parameters to the same value. It also sets the hierarchical, population
#' level parameters to constant values (this doesn't really have any effect:
#' if all of the participant level parameters are constant, the population level
#' parameters do nothing).
#' 
#' If \code{doConditionEffects} is TRUE, the condition effects are all set to 0, which
#' means that the participants will have the same parameter values in all conditions.
#' Note that this overrides the \code{conditionEffects} setting in the
#' configuration for \code{\link{runParameterEstimation}}. If you allow condition
#' effects with \code{conditionEffects} but set \code{doConditionEffects == TRUE}
#' you will get no condition effects. The reverse is not true.
#' 
#' If you call this function several times to get multiple lists of fixed parameter
#' values, know that you can combine the lists with the concatenate function, \code{c()}.
#' 
#' @param data The data you will use, in the same format as required by \code{\link{runParameterEstimation}}.
#' @param param The parameter to set to a constant value, e.g. "pMem".
#' @param value The constant value of the paramter. See also the \code{transformValueToLatent} argument.
#' @param doConditionEffects Whether condition effects should all be set to 0.
#' @param transformValueToLatent Whether \code{value} should be transformed to the latent space. If TRUE (the default) you should provide parameter values in the manifest space (e.g. probabilities should be between 0 and 1). If FALSE, you should provide parameter values in the latent space.
#' 
#' @return A list mapping from parameter name to parameter value.
#' 
#' @seealso \code{\link{setConstantCategoryParameters}} For setting catMu and catActive to constant values.
#' 
#' @examples
#' data = data.frame(pnum=rep(1:5, each=2), cond=rep(1:2, each=5))
#' constParam = setConstantParameterValue(data, param = "pContBetween", value = 0.5)
#' 
#' @export
setConstantParameterValue = function(data, param, value, doConditionEffects = TRUE, transformValueToLatent = TRUE) {
	
	if (transformValueToLatent) {
		trans = getParameterTransformation(param, NULL, inverse=TRUE)
		value = trans(value)
		if (param %in% getProbParams(NULL)) {
			value = min(max(value, -100), 100) #clamp probability parameters to be finite.
		}
	}
	
	pnums = sort(unique(data$pnum))
	conds = sort(unique(data$cond))
	
	rval = list()
	for (pnum in pnums) {
		partParam = paste( param, "[", pnum, "]", sep="")
		rval[[ partParam ]] = value
	}

	rval[[ paste( param, ".mu", sep="" ) ]] = 0
	rval[[ paste( param, ".var", sep="" ) ]] = 1
	
	if (doConditionEffects) {
		for (cond in conds) {
			condParam = paste( param, "_cond[", cond, "]", sep="")
			rval[[ condParam ]] = 0
		}
	}
	
	rval
}

#' Set catMu and catActive Parameters to Constant Values
#' 
#' 
#' @param data Your data set, in the same format as required by \code{\link{runParameterEstimation}}.
#' @param catParam A data.frame with constant values of catMu and catActive parameters. See details.
#' @param maxCategories The maximum number of categories that will be allowed when running parameter estimation.
#' @param activateConstantCats If TRUE, categories for which catMu values are given are forced to be active.
#' @param deactivateUnspecifiedCats If TRUE, catActive is set to 0 for all categories for which constant values of both catMu and catActive are not specified. If FALSE, catActive and catMu are left as free parameters for the extra categories.
#' 
#' @return A list mapping from parameter name to parameter value.
#' 
#' @details The catParam data.frame can have three columns: pnum, catMu, and catActive, where the catActive column is optional. Each row specifies the setting for a catMu and/or catActive for the given participant. Each participant can have multiple rows, each row specifying constant parameters for a different category. If you want to allow a parameter to be estimated, set the value for that parameter to NA. 
#' 
#' Because either or both of catMu and catActive can be specified, there are 4 possibilities per category: 
#' 1) both catMu and catActive set to constant values. This forces the model to use the category at a fixed location.
#' 2) catMu set to a constant value but catActive freely estimated. In this case, what catActive will tell you is mow much of the time a category at that fixed location was used. You might do this if you have candidates for category locations that should be used, but don't know how much each category might be used.
#' 3) catActive set to a constant value, but catMu freely estimated. You should not set catActive to 0 in this case, because inactive categories do nothing with the catMu parameters. Thus, by setting catActive to 1 but leaving catMu unspecified, it forces the model to use a category, but allows it to put that category anywhere.
#' 4) Neither catMu nor catActive specified. This is the default behavior of the model, where both are freely estimated.
#' 
#' Each of these 4 possibilities is set per category. Thus, you can do complex things like force each participant to have 3 active categories in fixed locations, 2 more active categories in freely-estimated locations, and some additional number of fully freely-estimated categories.
#' 
#' @examples 
#' catParam = data.frame(pnum = rep(1, 4),
#' 	catMu = c(60, 120, NA, NA),
#' 	catActive = c(1, NA, 1, NA))
#' data = data.frame(pnum=1) #This is just for testing
#' constCat = setConstantCategoryParameters(data, catParam, maxCategories = 5, 
#' 	activateConstantCats = TRUE, deactivateUnspecifiedCats = TRUE)
#' 
#' @export
setConstantCategoryParameters = function(data, catParam, maxCategories, activateConstantCats = TRUE, 
														deactivateUnspecifiedCats = TRUE) 
{

	allPnums = sort(unique(data$pnum))

	if (is.null(catParam$catActive)) {
		catParam$catActive = NA
	}
	
	if (activateConstantCats) {
		catParam[ !is.na(catParam$catMu), "catActive" ] = 1.0
	}
	
	if (!all(catParam$catActive %in% c(0.0, 1.0, NA))) {
		stop("All catActive parameters must be either 0 or 1, or NA if not constant.")
	}

	rval = list()
	for (pnum in allPnums) {
		
		#mus = catParam$catMu[ catParam$pnum == pnum ]
		catMu = catParam$catMu[ catParam$pnum == pnum ]
		catActive = catParam$catActive[ catParam$pnum == pnum ]
		
		if (length(catMu) > maxCategories) {
			warning( paste("More categories provided for participant ", pnum, " than the maximum number of categories.", sep="") )
		}
		
		for (i in 1:maxCategories) {
			
			catInd = i - 1 #zero-indexed
			
			catMuName = paste( "catMu[", pnum, ",", catInd, "]", sep="")
			catActiveName = paste( "catActive[", pnum, ",", catInd, "]", sep="")
			
			if (!is.na(catMu[i])) {
				rval[[ catMuName ]] = catMu[i]
			}
			
			if (!is.na(catActive[i])) {
				rval[[ catActiveName ]] = catActive[i]
			}
			
			if (deactivateUnspecifiedCats && is.na(catActive[i]) && is.na(catActive[i])) {
				rval[[ catMuName ]] = 0.0
				rval[[ catActiveName ]] = 0.0
			}

		}
	}

	rval
}

upgradeResultsList = function(results, from, to) {
	
	if (from == "0.6.1" && to == "0.7.0") {
		results$config$conditionEffects = list()
		for (pp in results$config$parametersWithConditionEffects) {
			results$config$conditionEffects[[pp]] = "all"
		}
		results$config$parametersWithConditionEffects = NULL
		
		results$conditions = NULL
		
		results$config = verifyConfigurationList(results$config, results$data)
	} else {
		stop("Unsupported upgrade path")
	}
	
	results
}

