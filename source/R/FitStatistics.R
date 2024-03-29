
#' Calculate WAIC from a Likelihood Matrix
#' 
#' This is a general-purpose WAIC calculation function that can be used with any appropriate model. 
#' It is not specific to CatContModel and can be used with other models.
#' For information about WAIC, see Bayesian Data Analysis, Third Edition by Gelman, Carlin, Stern, Dunson, Vehtari, and Rubin.
#' 
#' @param likeMat A matrix of likelihoods (not log). Each row corresponds to one observation in the data set (e.g. a study/response pair).
#' Each column corresponds to one sample from the posterior of the parameters (e.g. one iteration of a Gibbs sampler).
#' 
#' @return A list containing elements:
#' \tabular{ll}{
#' 	`WAIC_1`, `WAIC_2` \tab WAIC values that are calculated using two different methods of estimating the effective number of free parameters (`P_1` and `P_2`).\cr
#'  `P_1`, `P_2` \tab Estimated effective number of free parameters used to calculate `WAIC_1` and `WAIC_2`, respectively. \cr
#' 	`LPPD` \tab Log Posterior Predictive Density. Used to calculate both `WAIC_1` and `WAIC_2`.
#' }
#' 
#' @export
calculateWAIC_likelihoodMatrix = function(likeMat) {
  
  # Initialize accumulator variables that will be added to
  LPPD = 0
  P_1 = 0
  P_2 = 0
  
  for (obs in 1:nrow(likeMat)) {
    # A vector of likelihoods for a single observation, length equal to number of iterations
    L = likeMat[obs,]
    
    # Calculate LPPD for this observation
    LPPD = LPPD + log( mean(L) )
    
    # Calculate P_1
    a = log( mean(L) )
    b = mean( log(L) )
    
    P_1 = P_1 + 2 * (a - b)
    
    # Calculate P_2
    P_2 = P_2 + stats::var( log(L) )
    
  }
  
  # The two versions of WAIC differ in terms of which P is used
  WAIC_1 = -2 * (LPPD - P_1)
  WAIC_2 = -2 * (LPPD - P_2)
  
  # Return a list containing the results
  list(WAIC_1 = WAIC_1, WAIC_2 = WAIC_2, P_1 = P_1, P_2 = P_2, LPPD = LPPD)
  
}

calcWAIC_MF_singleSubsample = function(res) {
  
  postMat = convertPosteriorsToMatrix(res, stripConstantParameters = FALSE, stripCatActive = FALSE, stripCatMu = FALSE)
  
  likeMat = CCM_CPP_MF_batchLikelihood(res$config, res$data, postMat)

  # Calculate WAIC per participant
  allWAIC = NULL
  for (pnum in res$pnums) {
    rows = which(res$data$pnum == pnum)
    waic_part = calculateWAIC_likelihoodMatrix(likeMat[rows,])
    
    waic_part$pnum = pnum
    
    allWAIC = rbind(allWAIC, as.data.frame(waic_part))
  }
  
  # Sum across participants to get the total
  totals = apply(allWAIC[ , c("WAIC_1", "WAIC_2", "P_1", "P_2", "LPPD") ], 2, sum)
  totals = as.list(totals)
  totals$pnum = "Total"
  totals = as.data.frame(totals)
  
  #if (onlyTotal) {
  #  allWAIC = totals
  #} else {
  allWAIC = rbind(totals, allWAIC)
  #}
  
  allWAIC
  
}


#' Calculate Whole-Model WAIC
#' 
#' WAIC is a whole-model fit statistic, like AIC and BIC. WAIC cannot be directly compared with AIC or BIC, but it is conceptually very similar.
#' 
#' WAIC is an appropriate fit statistic for these models because
#' 1) it can be calculated without using posterior means, which some parameters do not have (specifically `catMu` and `catActive`). 
#' 2) it estimates the effective number of free parameters, which is wildly different from the actual number of free parameters.
#' See the documentation for [`calculateInappropriateFitStatistics`] for more information on the number of parameters.
#' 
#' There are two ways to estimate the effective number of free parameters for WAIC and results from both are reported as `WAIC_1` and `WAIC_2`.
#' The estimated number of free parameters are reported as `P_1` and `P_2`.
#' Both ways of calculating WAIC use the same `LPPD` or "log posterior predictive density".
#'  
#' WAIC is a whole-model fit statistic. These models have parameters that are shared by participants (the hierarchical 
#' population mean and variance parameters and the condition effect parameters). This means that the participants are 
#' not independent, so examining WAIC for individual participants is unprincipled and may give inaccurate results. 
#' Thus, you should leave `onlyTotal` at the default value of `TRUE`.
#' 
#' The main use for WAIC is to compare model variants to select an appropriate model for a data set 
#' (i.e. does the between-item or ZL model fit this data best?). 
#' WAIC can also be used to compare models that differ in other ways, such as which parameters have condition effects, 
#' the maximum number of categories, reduced models with some parameters set to constant values, more or less restrictive priors, etc.
#' 
#' You can estimate how much the variability in the posterior chains affects WAIC by using the `subsamples` and 
#' `subsampleProportion` arguments. If `subsamples` is greater than 1, multiple subsamples from the posterior chains 
#' will be taken and the standard deviation of WAIC (and its components) across the subsamples will be calculated. 
#' The number of iterations used in each subsample is a proportion of the total number of iterations and is set by 
#' `subsampleProportion`. Note that this is not the standard deviation of WAIC over repeated samples of data sets, 
#' so it tells you nothing about what would happen if you had different data. It essentially tells you whether or 
#' not you ran enough iterations to have a stable WAIC estimate. The closer `subsampleProportion` is to 1, the less 
#' independent the subsamples will be, so you should use a reasonably low value of `subsampleProportion`. The degree to 
#' which the subsamples are independent influences to what extent the standard deviation is underestimated: The less 
#' independent, the larger the underestimate will be. If you want fully independent subsamples, you can set 
#' `subsampleProportion` to `NULL`. However, this means that the number of subsamples and the proportion of iterations 
#' in each subsample to be inversely related, which means that you have to choose between a low number of subsamples or
#' a low number of iterations per subsample.
#' 
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param subsamples Number of subsamples of the posterior chains to take. If greater than 1, subsampleProportion should be set to a value between 0 and 1 (exclusive).
#' @param subsampleProportion The proportion of the total iterations to include in each subsample. This should probably only be less than 1 if `subsamples` is greater than 1. If `NULL`, `subsampleProportion` will be set to `1 / subsamples` and no iterations will be shared between subsamples (i.e. each subsample will be independent, except inasmuch as there is autocorrelation between iterations).
#' @param onlyTotal If `TRUE`, exclude participant-level WAIC values (which aren't really valid because of the fact that participants are not independent). I recommend you leave this at `TRUE`.
#' @param summarize Logical or `NULL`. Whether WAIC values from multiple subsamples should be summarized. If `NULL` (default), results will be summarized only if there are multiple subsamples.
#' @param progress Boolean. If `TRUE`, a text progress bar is shown.
#' @param impl Implementation of WAIC calculation. `PC`: Pure c++ (faster). `BL_LM`: Batch likelihood, then likelihood matrix (slower).
#' 
#' @return A data.frame containing WAIC values, estimates of the effective number of free parameters, and LPPD.
#' \tabular{ll}{
#'  `pnum` \tab If `onlyTotal == FALSE`, the participant related to the statistic. \cr
#'  `group` \tab If using a BP design, the group that the results are from. If `"all"`, that is the sum of the values across all groups. \cr
#' 	`stat` \tab The name of the fit statistic.\cr
#' 	`value` \tab Only if the results are not summarized. The value of the statistic. \cr
#' 	`mean` \tab This and the following columns are provided if the results are summarized. The mean statistic value. \cr
#' 	`sd` \tab The standard deviation of the statistic. \cr
#' 	`min, median, max` \tab The minimum, median, and maximum of the statistic. \cr
#' 	`p2.5, p97.5` \tab The 2.5 and 97.5 percentiles of the statistic. 	
#' }
#' 
#' @family generic functions
#' @seealso [calculateWAIC_likelihoodMatrix()]
#' 
#' @export
calculateWAIC = function(res, subsamples = 1, subsampleProportion = 1, onlyTotal = TRUE, summarize = NULL, 
												 progress = subsamples > 1, impl=c("PC", "BL_LM")) {
	
	summarize = valueIfNull(summarize, subsamples > 1)
	
	# catMu is always degrees in R now
	#res = calculateWAIC_convertPosteriorCatMu(res)
	
	allWAIC = calculateWAIC_subsamples(res, subsamples, subsampleProportion, progress = progress, impl = impl[1])
	agg = calculateWAIC_aggregate(res, allWAIC, onlyTotal = onlyTotal, summarize = summarize)
	agg
	
}



calculateWAIC_subsamples = function(res, subsamples, subsampleProportion, progress, impl) {
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(res$runConfig$iterations, subsamples, subsampleProportion)
	
	allWAIC = NULL
	
	if (subsamples > 1) {
	  pb = utils::txtProgressBar(0, 1, 0, style=3)
	}
	
	for (sub in 1:length(subsampleIterationsToRemove)) {
		if (length(subsampleIterationsToRemove[[sub]]) > 0) {
			subRes = removeBurnIn(res, subsampleIterationsToRemove[[sub]])
		} else {
			subRes = res
		}
		
		if (resultIsType(subRes, "WP")) {
			groups = list(default = subRes)
		} else if (resultIsType(subRes, "BP")) {
			groups = subRes$groups
		}
		
		for (grp in names(groups)) {
		  
		  thisSubsample = groups[[ grp ]]
		  
		  # Call the WAIC implementation
		  if (impl == "BL_LM") {
		    colWAIC = calcWAIC_MF_singleSubsample(thisSubsample)
		  } else if(impl == "PC") {
		    colWAIC = CCM_CPP_calculateWAIC(thisSubsample)
		  }
			
			statNames = names(colWAIC)
			statNames = statNames[ statNames != "pnum" ]
			
			for (sn in statNames) {
				temp = data.frame(subsample = sub, group = grp, pnum = colWAIC$pnum, stat = sn, value = colWAIC[ , sn ], stringsAsFactors = FALSE)
				allWAIC = rbind(allWAIC, temp)
			}
			
		}
		
	  if (subsamples > 1) {
		  utils::setTxtProgressBar(pb, sub / length(subsampleIterationsToRemove))
	  }
	}
	
	if (subsamples > 1) {
	  close(pb)
	}
	
	allWAIC
}


calculateWAIC_aggregate = function(res, allWAIC, onlyTotal, summarize = TRUE) {

	summaryFuns = list(mean = mean, 
										 sd = stats::sd, 
										 min = min, 
										 p2.5 = function(x) { stats::quantile(x, 0.025) },
										 median = stats::median, 
										 p97.5 = function(x) { stats::quantile(x, 0.975) },
										 max = max)
	if (!summarize) {
		summaryFuns = list(value = function(x) { x })
	}
	
	# For BP designs, sum total WAIC across groups to get the all-group total
	if (resultIsType(res, "BP")) {
		totalWAIC = allWAIC[ allWAIC$pnum == "Total", ]
		totalWAIC$group = "all"
		
		summedTotal = stats::aggregate(value ~ subsample * group * pnum * stat, totalWAIC, sum)
	
		allWAIC = rbind(allWAIC, summedTotal)
	}

	# Aggregate, collapsing across subsamples
	summaryValues = stats::aggregate(value ~ stat * group * pnum, allWAIC, function(x) { NA })
	summaryValues$value = NULL
	
	for (fn in names(summaryFuns)) {
		temp = stats::aggregate(value ~ stat * group * pnum, allWAIC, summaryFuns[[ fn ]])
		summaryValues[ , fn ] = temp$value
	}
	
	if (onlyTotal) {
		summaryValues = summaryValues[ summaryValues$pnum == "Total", ]
		summaryValues$pnum = NULL
	}
	
	if (resultIsType(res, "WP")) {
		summaryValues$group = NULL
	}
	
	summaryValues
}


# Depreciated
# Convert catMu from degrees to radians if circular.
calculateWAIC_convertPosteriorCatMu = function(res) {
  
  if (res$config$dataType != "circular" || res$config$maxCategories == 0) {
    return(res)
  }
  
  if (resultIsType(res, "WP")) {
    res = convertCatMuUnits(res, CatContModel::d2r)
    
  } else if (resultIsType(res, "BP")) {
    for (grp in names(res$groups)) {
      res$groups[[ grp ]] = convertCatMuUnits(res$groups[[ grp ]], CatContModel::d2r)
    }
  } else {
    stop("Invalid results type.")
  }
  
  res
}


#' Calculate Standard Fit Statistics that are Inappropriate for these Models
#' 
#' AIC and BIC are fit statistics that are commonly used to compare models, but they are 
#' inappropriate for the models fit by this package (with the possible exception of the 
#' ZL model). The reason is that both AIC and BIC use the number of free parameters in the 
#' model as a penalty term. However, the actual number of free parameters and the effective 
#' number of free parameters are very different for these models. Part of the reason is that 
#' the `catActive` parameters don't provide very much flexibility to the model but are counted 
#' as free parameters. In addition, `catMu` does
#' nothing if the category is inactive, which is common. Finally, the hierarchical nature
#' of the model means that the population level parameters actually constrain the participant
#' level parameters, reducing the effective number of free parameters, but AIC and BIC count 
#' the hierarchical parameters as free parameters. All of this combined, the effective number 
#' of free parameters is far less than the actual number of free parameters.
#' 
#' An appropriate fit statistic for these models is WAIC, which estimates the effective 
#' number of free parameters. See the [`calculateWAIC`] function of this package.
#' WAIC is relatively straightforward to calculate for models for which Bayesian parameter 
#' estimation was done.
#' 
#' Although AIC and BIC have an overly heavy penalty term, there can still be some value in using
#' those statistics to compare to other models in a back-of-the-envelope way.
#' The way that this function calculates the fit statistics is by evenly dividing the number
#' of free parameters in the whole between the participants, so each participant has the same
#' penalty term. For each iteration of the Gibbs sampler, the likelihoods of the model are
#' calculated for each participant. Thus, for each iteration, it is possible to calculate
#' AIC and BIC, which allows for an estimate of the uncertainty in the fit statistics.
#' For each fit statistic, both mean and standard deviation are reported.
#' 
#' @param results The results from the [`runParameterEstimation`] function. Note that you must set `config$calculateParticipantLikelihoods` to `TRUE` in order to use this function.
#' @param onlyTotal If `TRUE`, exclude participant-level values (which aren't valid because of the fact that participants are not independent).
#' 

#' @return A `data.frame` containing several fit statistics.
#' 
#' @seealso [`calculateWAIC`] for an appropriate fit statistic.
#' 
#' @family WP functions
#'
#' @export
calculateInappropriateFitStatistics = function(results, onlyTotal = TRUE) {
	
	if (!resultIsType(results, "WP")) {
		stop("calculateInappropriateFitStatistics only accepts WP results objects. See the Glossary (listed in the package functions index).")
	}
	
	if (!results$config$calculateParticipantLikelihoods) {
		stop("In order to calculate the inappropriate fit statistics, you need to have calculated the participant likelihoods. Set 'calculateParticipantLikelihoods' to TRUE in the config and rerun the parameter estimation.")
	}
	
	logWarning("AIC and BIC fit summaries are inappropriate for these models because the actual number of free parameters and the effective number of free parameters are very different. You should use WAIC instead. See the calculateWAIC function.")
	
	
	calcFits = function(pnum, ll, nParam, nObs) {
		aic_p = -2 * (ll - nParam)
		
		bic_p = -2 * ll + log(nObs) * nParam
		
		data.frame(pnum = pnum, nParam = nParam, nObs = nObs,
							 ll_mean = mean(ll), ll_sd = stats::sd(ll),
							 aic_mean = mean(aic_p), aic_sd = stats::sd(aic_p),
							 bic_mean = mean(bic_p), bic_sd = stats::sd(bic_p))
	}
	
	post = convertPosteriorsToMatrices(results, "catActive")
	
	#count free parameters
	anyCatActiveChanged = any( apply(post$catActive, c(1,2), function(x) {any(x != x[1])}) )
	nFreeParam = 0
	for (n in names(results$posteriors)) {
		x = results$posteriors[[n]]
		
		isConstant = all(x == x[1]) || all(is.na(x))
		
		isLL = startsWith(n, "participantLL")
		isCatActive = startsWith(n, "catActive")
		
		if (isCatActive) {
			isConstant = !anyCatActiveChanged
		}
		
		if (!isConstant && !isLL) {
			nFreeParam = nFreeParam + 1
		}
	}
	
	nParticipants = length(results$pnums)
	
	nFreeParamPerParticipant = nFreeParam / nParticipants #may be fractional
	
	totalLL = 0
	fitSummary = NULL
	
	for (n in results$pnums) {
		
		nObs = sum(results$data$pnum == n)
		
		ll = results$posteriors[[ paste("participantLL[", n, "]", sep="") ]]
		totalLL = totalLL + ll
		
		temp = calcFits(n, ll, nParam = nFreeParamPerParticipant, nObs=nObs)
		
		fitSummary = rbind(fitSummary, temp)
		
	}
	
	temp = calcFits("Total", totalLL, nParam = nFreeParam, nObs = nrow(results$data))
	
	fitSummary = rbind(fitSummary, temp)
	
	if (onlyTotal) {
		fitSummary = fitSummary[ fitSummary$pnum == "Total", ]
	}
	
	fitSummary
}



