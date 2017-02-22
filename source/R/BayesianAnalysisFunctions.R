
#' Posterior Means and Credible Intervals for Participant Parameters
#' 
#' The credible intervals for catActive should be interpreted in the context of them only taking on integer values.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param params A vector of parameter names. If \code{NULL}, all valid parameters are used.
#' @param doCatActive If \code{TRUE}, the posterior mean and credible intervals of the \code{catActive} parameters will be calculated.
#' @param credLevel The credibility level of the credible intervals. Defaults to 0.95.
#' @param fun A user provided function that will be passed a vector of a single participant by condition posterior distribution and that should return a single value.
#' 
#' @return A data frame with several columns:
#' \tabular{ll}{
#'	\code{pnum} \tab The participant number.\cr
#' 	\code{cond} \tab The zero-indexed condition index.\cr
#' 	\code{param} \tab The parameter name.\cr
#' 	\code{mean} \tab The posterior mean.\cr
#' 	\code{ciLower,ciUpper} \tab The lower and upper bounds of the credible interval.\cr
#' 	\code{fun} \tab If \code{fun} was provided, the results of that function call.
#' }
#' 
#' @export
participantPosteriorSummary = function(results, params=NULL, doCatActive=TRUE, credLevel = 0.95, fun=NULL) {
	
	resDf = NULL
	
	if (is.null(params)) {
		params = c(getProbParams(results, filter=TRUE), getSdParams(results, filter=TRUE))
	}
	
	halfa = (1 - credLevel) / 2
	
	for (param in params) {
		
		hasConditionEffect = (param %in% getParametersWithConditionEffects(results$config$conditionEffects))
		
		conditions = results$config$factors$cond
		if (!hasConditionEffect) {
			conditions = results$config$cornerstoneConditionName
		}
		
		for (cond in conditions) {
			
			for (pnum in results$pnums) {
				
				partParam = getParameterPosterior(results, param, pnum=pnum, cond=cond, manifest=TRUE)
				
				qs = as.numeric( stats::quantile(partParam, c(halfa, 1 - halfa)) )
				
				temp = data.frame(pnum=pnum, cond=cond, param=param, 
													mean = mean(partParam), ciLower = qs[1], ciUpper = qs[2], stringsAsFactors = FALSE)
				
				if (!is.null(fun)) {
					temp$fun = fun(partParam)
				}
				
				resDf = rbind(resDf, temp)
			}
		}
	}
	
	post = convertPosteriorsToMatrices(results, "catActive")
	
	for (pnum in results$pnums) {
		
		ca = post$catActive[ which(results$pnums == pnum), ,  ]
		
		itCount = apply(ca, 2, sum)
		
		qs = as.numeric( stats::quantile(itCount, c(halfa, 1 - halfa)) )
		
		temp = data.frame(pnum=pnum, cond=results$config$cornerstoneConditionName, param="catActive", 
											mean = mean(itCount), ciLower = qs[1], ciUpper = qs[2], stringsAsFactors = FALSE)
		
		if (!is.null(fun)) {
			temp$fun = fun(itCount)
		}
		
		resDf = rbind(resDf, temp)
		
	}
	
	resDf = resDf[ order(resDf$pnum, resDf$param, resDf$cond), ]
	
	resDf
}

#' Population/Condition Posterior Means and Credible Intervals
#' 
#' Calculates posterior means and credible intervals for the population means in each condition for the
#' given parameters. For each condition, condition effects are added to population means, the result is 
#' transformed to the manifest space, and the mean and credible interval for the manifest value is calculated.
#' Note that this is different from adding condition effects to participant-level parameters, tranforming 
#' the result, calculating on each iteration the mean of the transformed participant parameters, and 
#' calculating the posterior mean and credible interval of the iteration means. 
#' Using iteration means rather than population means will generally result in less than the true
#' amount of variability, which is why population means are used. Note, however, that this is a little strange,
#' because in the model, condition effects are not added to population means, but participant means.
#' 
#' In the return value of this function, the lower and upper columns give the endpoints of the credible interval.
#' 
#' Note that parameters without condition effects do not have values specific to conditions, so the cond will be NA for those parameters.
#'
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param params A vector of parameter names. If NULL, the default, all valid parameters are used.
#' @param credLevel The credibility level of the credible intervals. Defaults to 0.95.
#' 
#' @return A data.frame containing the results.
#' 
#' @export
posteriorMeansAndCredibleIntervals = function(results, params=NULL, credLevel=0.95) {
	
	if (is.null(params)) {
		params = c(getProbParams(results, filter=TRUE), getSdParams(results, filter=TRUE))
	}
	
	rval = NULL
	
	halfa = (1 - credLevel) / 2
	
	for (param in params) {
		
		trans = getParameterTransformation(param, results)
		
		mu = results$posteriors[[ paste(param, ".mu", sep="") ]]
		
		for (cond in results$config$factors$cond) {
			
			
			condEff = results$posteriors[[ paste(param, "_cond[", cond, "]", sep="") ]]
			
			combined = trans(mu + condEff)
			
			quants = as.vector(stats::quantile(combined, c(halfa, 1 - halfa)))
			
			temp = data.frame(param=param, cond=cond, mean=mean(combined), lower=quants[1], upper=quants[2], stringsAsFactors=FALSE)
			
			paramWithConditionEffects = getParametersWithConditionEffects(results$config$conditionEffects)
			if (!(param %in% paramWithConditionEffects)) {
				temp$cond = "NA"
			}
			
			#For parameters without condition effects, only include one condition 
			#(the cornerstone condition, but it doesn't matter)
			if (param %in% paramWithConditionEffects || 
					cond == results$config$cornerstoneConditionName) 
			{
				rval = rbind(rval, temp)
			}
			
		}
	}
	
	rval = rval[ order(rval$param, rval$cond), ]
	
	rval
}


#' Transformed Participant Level Parameters for a Participant, Condition, and Iteration
#' 
#' Note that inactive categories have already been removed, which is why the return value does not include the `catActive` parameters.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param pnum A single participant number.
#' @param cond A single condition name.
#' @param iteration The index of an iteration of the chain.
#' @param removeInactiveCategories Whether catMu parameters for inactive categories should be removed. If FALSE, both catMu and catActive parameters are provided.
#' 
#' @return A list containing transformed parameters.
#' 
#' @export
getTransformedParameters = function(results, pnum, cond, iteration, removeInactiveCategories = TRUE) {
	
	allParam = c(getProbParams(results, filter=FALSE), getSdParams(results, filter=FALSE))
	
	combinedParam = list()
	for (pp in allParam) {
		condParam = results$posteriors[[ paste(pp, "_cond[", cond, "]", sep="") ]][iteration]
		
		partParam = results$posteriors[[ paste(pp, "[", pnum, "]", sep="") ]][iteration]
		
		trans = getParameterTransformation(pp, results)
		
		combinedParam[[pp]] = trans(partParam + condParam)
	}
	
	#Do catMu
	ca = cm = rep(NA, results$config$maxCategories)
	for (i in 1:results$config$maxCategories) {
		istr = paste("[", pnum, ",", i - 1, "]", sep="")
		cm[i] = results$posteriors[[ paste("catMu", istr, sep="") ]][iteration]
		ca[i] = results$posteriors[[ paste("catActive", istr, sep="") ]][iteration]
	}
	
	if (removeInactiveCategories) {
		combinedParam$catMu = cm[ ca == 1 ]
	} else {
		combinedParam$catMu = cm
		combinedParam$catActive = ca
	}
	
	if (results$config$dataType == "circular") {
		combinedParam$catMu = combinedParam$catMu %% 360 #limit to the interval [0, 360)
	}
	
	combinedParam
}




#' Vectors of Posterior Parameter Chains
#' 
#' This function helps with retrieving posterior chains for parameters.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param The name of the parameter.
#' @param pnum The participant number. If NULL, the population mean of the parameter will be used.
#' @param cond The condition to use. If NULL, no condition effect will be applied.
#' @param manifest If TRUE, the default, the chain will be converted from the latent to the manifest space. If FALSE, it will be left in the latent space.
#' 
#' @return A vector containing the parameter chain.
#' 
#' @export
getParameterPosterior = function(results, param, pnum, cond, manifest=TRUE) {
	
	if (!is.null(pnum)) {
		param_base = results$posteriors[[ paste(param, "[", pnum, "]", sep="") ]]
	} else {
		param_base = results$posteriors[[ paste(param, ".mu", sep="") ]]
	}
	
	param_cond = rep(0, length(param_base))
	if (!is.null(cond)) {
		param_cond = results$posteriors[[ paste(param, "_cond[", cond, "]", sep="") ]]
	}
	
	param_latent = param_base + param_cond
	
	transformationFunction = getParameterTransformation(param, results)
	if (!manifest) {
		transformationFunction = function(x) { x }
	}
	transformationFunction(param_latent)
	
}



#' Sample Data from Model with Specific Parameter Values
#' 
#' Samples data from the model given the provided parameter values. This is useful for observing the patterns of data generated by the model and for sampling from the posterior predictive distribution of the data.
#' 
#' @param study A vector of study angles in degrees.
#' @param param A list of parameter values. The values to include are \code{pMem}, \code{pBetween}, \code{pContBetween}, \code{pContWithin}, \code{pCatGuess}, \code{contSD}, \code{catMu}, \code{catSelectivity}, \code{catSD}. \code{catMu} should be a vector. If there are no categories, \code{catMu} should be \code{NULL}.
#' @param modelVariant The model variant, as a string. Should be one of "betweenItem", "withinItem", and "ZL".
#' @param dataType One of \code{"circular"} or \code{"linear"}.
#' @param responseRange A length 2 vector giving the theoretical minimum and maximum values of a response. Should be provided if \code{dataType} is \code{"linear"}.
#' 
#' @return A data frame containing the \code{study} angles, the sampled \code{response} angles, the response \code{type} (e.g. continuous memory response), and, if the response was categorical in nature, the category it was from (\code{cat}).
#'  
#' @export
sampleDataFromModel = function(study, param, modelVariant, dataType = "circular", responseRange = NULL) {
	
	trials = length(study)
	
	if (modelVariant == "betweenItem") {
		param$pBetween = 1
	} else if (modelVariant == "withinItem") {
		param$pBetween = 0
	} else if (modelVariant == "ZL") {
		param$pBetween = 1
		param$pContBetween = 1
		param$pCatGuess = 0
		param$catMu = NULL
	}
	
	realzationFunction = NULL
	unifGuessFunction = NULL
	if (dataType == "circular") {
		
		realzationFunction = function(mu, sd) {
			CatContModel::rvmd(1, mu, sd)
		}
		unifGuessFunction = function() {
			stats::runif(1, 0, 360)
		}
		
	} else if (dataType == "linear") {
		
		realzationFunction = function(mu, sd) {
			msm::rtnorm(1, mu, sd, responseRange[1], responseRange[2])
		}
		unifGuessFunction = function() {
			stats::runif(1, responseRange[1], responseRange[2])
		}
	}
	
	
	catCount = length(param$catMu)

	
	response = rep(0, trials)
	
	cat = rep(0, trials)
	type = rep("none", trials)
	
	for (i in 1:trials) {
		
		inMemory = (stats::rbinom(1, 1, param$pMem) == 1)
		
		if (!inMemory) {
			
			isCatGuess = (stats::rbinom(1, 1, param$pCatGuess) == 1)
			
			if (isCatGuess && (catCount > 0)) {
				cat[i] = sample(1:catCount, size=1)
				
				response[i] = realzationFunction(param$catMu[cat[i]], param$catSD)
				type[i] = "catGuess"
			} else {
				response[i] = unifGuessFunction()
				type[i] = "unifGuess"
			}
			
		} else {
			
			
			#Pick an error for the continuous representation
			
			contLocation = realzationFunction(study[i], param$contSD)
			
			if (dataType == "circular") {
				contLocation = contLocation %% 360
			}

			
			if (catCount == 0) {
				response[i] = contLocation
				type[i] = "continuous"
				
			} else {
				
				#Pick a category for the color
				catWeights = categoryWeightsFunction(study[i], param$catMu, param$catSelectivity, dataType = dataType)
				catEgory = sample(1:catCount, size=1, prob=catWeights)
				cat[i] = catEgory
				
				
				# Get the categorical response location
				catLocation = realzationFunction(param$catMu[catEgory], param$catSD)
				if (dataType == "circular") {
					catLocation = catLocation %% 360
				}
				
				isBetween = (stats::rbinom(1, 1, param$pBetween) == 1)
				
				if (isBetween) {
					#This is a between response
					isContinuous = (stats::rbinom(1, 1, param$pContBetween) == 1)
					
					if (isContinuous) {
						response[i] = contLocation
						type[i] = "continuous"
					} else {
						response[i] = catLocation
						type[i] = "categorical"
					}
				} else {
					#This is a within response
					
					locMix = NA
					if (dataType == "circular") {
						locMix = circMean(c(contLocation, catLocation), c(param$pContWithin, 1 - param$pContWithin))
						locMix = locMix %% 360
					} else if (dataType == "linear") {
						locMix = param$pContWithin * contLocation + (1 - param$pContWithin) * catLocation
					}
					
					response[i] = locMix
					
					type[i] = "within"
				}
			}
		}
		
	}
	data.frame(study=study, response=response, type=type, cat=cat)
}


