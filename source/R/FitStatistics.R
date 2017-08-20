
#' Calculate Whole Model WAIC
#' 
#' WAIC is a whole-model fit statistic, like AIC. However, WAIC cannot be directly compared with AIC, but it is conceptually very similar.
#' 
#' WAIC is an appropriate fit statistic for these models because
#' 1) it can be calculated without using posterior means, which some parameters do not have (specifically catMu and catActive). 
#' 2) it estimates the effective number of free parameters, which is wildly different from the actual number of free parameters.
#' See the documentation for [`calculateInappropriateFitStatistics`] for more information on the number of parameters.
#' 
#' There are two ways to estimate the effective number of free parameters for WAIC and results from both are reported.
#'  
#' WAIC is a whole-model fit statistic. These models have parameters that are shared by participants (the hierarchical population mean and variance parameters and the condition effect parameters). This means that the participants are not independent, so examining WAIC for individual participants is unprincipled and may give inaccurate results. Thus, you should leave `onlyTotal` at the default value of `TRUE`.
#' 
#' Note that you can use WAIC to compare model variants, like the between-item and within-item variants. You can also compare models that differ in other ways, such as which parameters have condition effects, the maximum number of categories, reduced models with some parameters set to constant values, more or less restrictive priors, etc.
#' 
#' You can estimate how much the variability in the posterior chains affects WAIC by using the `subsamples` and `subsampleProportion` arguments. If `subsamples` is greater than 1, multiple subsamples from the posterior chains will be taken and the standard deviation of WAIC (et al.) across the subsamples will be calculated. The number of iterations used in each subsample is a proportion of the total number of iterations and is set by `subsampleProportion`. Note that this is not the standard deviation of WAIC over repeated samples of data sets, so it tells you nothing about what would happen if you had different data. It essentially tells you whether or not you ran enough iterations to have a stable WAIC estimate. The closer `subsampleProportion` is to 1, the less independent the subsamples will be, so you should use a reasonably low value of `subsampleProportion`. The degree to which the subsamples are independent influences to what extent the standard deviation is underestimated: The less independent, the larger the underestimate will be. If you want fully independent subsamples, you can set `subsampleProportion` to NULL. However, this means that the number of subsamples and the proportion of iterations in each subsample to be inversely related, which means that you have to choose between a low number of subsamples or a low number of iterations per subsample.
#' 
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param subsamples Number of subsamples of the posterior chains to take. If greater than 1, subsampleProportion should be set to a value between 0 and 1 (exclusive).
#' @param subsampleProportion The proportion of the total iterations to include in each subsample. This should probably only be less than 1 if `subsamples` is greater than 1. If `NULL`, `subsampleProportion` will be set to `1 / subsamples` and no iterations will be shared between subsamples (i.e. each subsample will be independent, except inasmuch as there is autocorrelation between iterations).
#' @param onlyTotal If `TRUE`, exclude participant-level WAIC values (which aren't really valid because of the fact that participants are not independent). I recommend you leave this at `TRUE`.
#' @param summarize Logical or `NULL`. Whether WAIC values from multiple subsamples should be summarized. If `NULL` (default), results will be summarized only if there are multiple subsamples.
#' 
#' @return A data.frame containing WAIC values, estimates of the effective number of free parameters, and LPPD.
#' \tabular{ll}{
#'  \code{pnum} \tab If \code{onlyTotal == FALSE}, the participant related to the statistic. \cr
#'  \code{group} \tab If using a BP design, the group that the results are from. If `"all"`, that is the sum of the values across all groups. \cr
#' 	\code{stat} \tab The name of the fit statistic.\cr
#' 	\code{value} \tab Only if the results are not summarized. The value of the statistic. \cr
#' 	\code{mean} \tab This and the following columns are provided if the results are summarized. The mean statistic value. \cr
#' 	\code{sd} \tab The standard deviation of the statistic. \cr
#' 	\code{min, median, max} \tab The minimum, median, and maximum of the statistic. \cr
#' 	\code{p2.5, p97.5} \tab The 2.5 and 97.5 percentiles of the statistic. 	
#' }
#' 
#' @md
#' @export
calculateWAIC = function(res, subsamples = 1, subsampleProportion = 1, onlyTotal = TRUE, summarize = NULL) {
	
	summarize = valueIfNull(summarize, subsamples > 1)
	
	res = CatContModel:::calculateWAIC_convertPosteriorCatMu(res)
	allWAIC = CatContModel:::calculateWAIC_subsamples(res, subsamples, subsampleProportion)
	agg = CatContModel:::calculateWAIC_aggregate(res, allWAIC, onlyTotal = onlyTotal, summarize = summarize)
	agg
	
}


# Convert catMu from degrees to radians if circular.
calculateWAIC_convertPosteriorCatMu = function(res) {
	
	if (res$config$dataType != "circular" || res$config$maxCategories == 0) {
		return(res)
	}
	
	if (resultIsType(res, "WP")) {
		res = calculateWAIC_convertPosteriorCatMu.WP(res)
	} else if (resultIsType(res, "BP")) {
		for (grp in names(res$groups)) {
			res$groups[[ grp ]] = calculateWAIC_convertPosteriorCatMu.WP(res$groups[[ grp ]])
		}
	}
	
	res
}

calculateWAIC_convertPosteriorCatMu.WP = function(results) {
	
	for (p in unique(results$data$pnum)) {
		for (i in 1:results$config$maxCategories) {
			cmName = paste("catMu[", p, ",", i-1, "]", sep="")
			results$posteriors[[ cmName ]] = CatContModel::d2r(results$posteriors[[ cmName ]])
		}
	}
	
	results
}


calculateWAIC_cppWrapper = function(results) {
	CCM_CPP_calculateWAIC(results)
}

calculateWAIC_subsamples = function(res, subsamples, subsampleProportion) {
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(res$config$iterations, subsamples, subsampleProportion)
	
	allWAIC = NULL
	
	pb = utils::txtProgressBar(0, 1, 0, style=3)
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
			colWAIC = calculateWAIC_cppWrapper(groups[[ grp ]])
			statNames = names(colWAIC)
			statNames = statNames[ statNames != "pnum" ]
			
			for (sn in statNames) {
				temp = data.frame(subsample = sub, group = grp, pnum = colWAIC$pnum, stat = sn, value = colWAIC[ , sn ], stringsAsFactors = FALSE)
				allWAIC = rbind(allWAIC, temp)
			}
			
		}
		
		utils::setTxtProgressBar(pb, sub / length(subsampleIterationsToRemove))
	}
	close(pb)
	
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
	
	# SUm total WAIC across groups to get the total total
	totalWAIC = allWAIC[ allWAIC$pnum == "Total", ]
	totalWAIC$group = "all"
	
	summedTotal = stats::aggregate(value ~ subsample * group * pnum * stat, totalWAIC, sum)

	allWAIC = rbind(allWAIC, summedTotal)

	# Aggregate, collapsing across subsamples
	summaryValues = stats::aggregate(value ~ stat * group * pnum, allWAIC, function(x) { NA })
	summaryValues$value = NULL
	
	for (fn in names(summaryFuns)) {
		summaryValues[ , fn ] = stats::aggregate(value ~ stat * group * pnum, allWAIC, summaryFuns[[ fn ]])$value
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




#' Calculate Standard Fit Statistics that are Inappropriate for these Models
#' 
#' AIC and BIC are fit statistics that are commonly used to compare models, but they are 
#' inappropriate for the models fit by this package (with the possible exception of the 
#' ZL model). The reason is that both AIC and BIC use the number of free parameters in the 
#' model as a penalty term. However, the actual number of free parameters and the effective 
#' number of free parameters are very different for these models. Part of the reason is that 
#' the `catActive` parameters don't provide very much flexibility. In addition, `catMu` does
#' nothing if the category is inactive, which is common. Finally, the hierarchical nature
#' of the model means that the population level parameters actually constrain the participant
#' level parameters, reducing the effective number of free parameters, but AIC and BIC count 
#' them as free parameters. All of this combined, the effective number of free parameters
#' is far less than the actual number of free parameters.
#' 
#' An appropriate fit statistc for these models is WAIC, which estimates the effective 
#' number of free parameters. See the [`calculateWAIC`] function of this package.
#' WAIC is relatively straightforward to calculate for models for which Bayesian parameter 
#' estimation was done.
#' 
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
#' @md
#' @export
calculateInappropriateFitStatistics = function(results, onlyTotal = TRUE) {
	
	if (!results$config$calculateParticipantLikelihoods) {
		stop("In order to calculate the inappropriate fit statistics, you need to have calculated the participant likelihoods. Set 'calculateParticipantLikelihoods' to TRUE in the config and rerun the parameter estimation.")
	}
	
	warning("AIC and BIC fit summaries are inappropriate for these models because the actual number of free parameters and the effective number of free parameters are very different. You should use WAIC instead. See the calculateWAIC function.")
	
	
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



