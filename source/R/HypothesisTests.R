
# Resulting parameters for a Cauchy RV that is a linear combination of other Cauchy RVs.
# The general form of the combinations is
# y = sum_i x_i * c_i
# where y is the resulting Cauchy random variable,
# x are the input CRVs, and c are coefficients.
# For example, if subtracting two Cauchys, c = c(1, -1) (note that order matters).
# For example, if taking the mean of many RVs, each weight is 1 / number of RVs.
cauchyRvLinearCombination = function(locations, scales, coefs) {
	
	location = sum(coefs * locations)
	
	scale = sum(abs(coefs) * scales)
	
	list(location=location, scale=scale)
}


#postVect: vector of the posterior distribution of the parameter you want to test
#h0_val: The value of the parameter tested in the null hypothesis.
#priorDensAtH0: The density of the prior distribution at h0_val
savageDickey = function(postVect, h0_val, priorDensAtH0) {
	
	success = tryCatch({
		ls = polspline::logspline(postVect)
		TRUE
	}, error = function(e) {
		print(e)
		return(FALSE)
	})
	
	if (success) {
		postDens = polspline::dlogspline(h0_val, ls)
		
		bf10 = priorDensAtH0 / postDens
	} else {
		bf10 = NA
	}
	
	list(bf01 = 1 / bf10, bf10 = bf10, success = success)
}


#' Test the Amount of Categorical Responding Present in the Data
#' 
#' Test the amount of categorical responding present in the data in each condition.
#' To do this test, you must have run the "betweenItem" model variant. Then, use this function on the results from that run to test whether categorical responding is or is not present in the data. A lack of categorical responding is the same as fully continuous responding, which happens when pContBetween = 1 and pCatGuess = 0. As explained in Hardman, Vergauwe, and Ricker (2017), testing 1 and 0 is not possible, so values near 1 and 0 should be tested instead. The test is performed in each condition separately.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param pContBetween_test Value of pContBetween to test.
#' @param pCatGuess_test Value of pCatGuess to test.
#' 
#' @return A data.frame with the results of the tests.
#' 
#' @seealso testMeanParameterValue for a more general way to do these kinds of hypothesis tests.
#' 
#' @export
testCategoricalResponding = function(results, pContBetween_test = 0.99, pCatGuess_test = 0.01) {

	H0_test = list(pContBetween = pContBetween_test, pCatGuess = pCatGuess_test)
	
	rval = NULL
	for (param in c("pContBetween", "pCatGuess")) {
		for (cond in results$config$factors$cond) {
			
			res = testMeanParameterValue(results, param = param, cond = cond, H0_value = H0_test[[param]])

			rval = rbind(rval, res)
			
		}
	}
	
	rval
}


convolveFuns = function(f, g, t, range=c(-Inf, Inf)) {
	
	#all of the ways to get to t.
	h = function(tau) {
		f(tau) * g(t - tau)
	}
	
	dens = stats::integrate(h, lower=range[1], upper=range[2])
	
	dens$value
}


#' Test the Mean Value for a Parameter in a Condition
#' 
#' Perform the test of whether the mean of the parameter is equal to some value. 
#' This test is performed in a condition.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param The name of the parameter to test.
#' @param cond The condition in which to test the value.
#' @param H0_value The null hypothesized value of the parameter. This should be provided in the manifest space (e.g., if testing a probability parameter, this should be between 0 and 1).
#' 
#' @export
#' 
#' @examples \dontrun{
#' testMeanParameterValue(results, param = "pMem", cond = "1", H0_value = 0.95)
#' }
#' 
testMeanParameterValue = function(results, param, cond, H0_value) {
	
	transformInverse = getParameterTransformation(results, param, inverse=TRUE)
	
	inverseH0 = transformInverse(H0_value)
	
	popMu0 = results$priors[[ paste(param, ".mu.mu", sep="") ]]
	popSd0 = sqrt( results$priors[[ paste(param, ".mu.var", sep="") ]] )
	
	condPrior = getConditionParameterPrior(results, param, cond)
	
	priorDensAtH0 = 0
	if (condPrior$scale == 0) {

		priorDensAtH0 = stats::dnorm(inverseH0, popMu0 + condPrior$location, popSd0)
		
	} else {
		
		fNorm = function(x) {
			stats::dnorm(x, popMu0, popSd0)
		}
		fCauchy = function(x) {
			stats::dcauchy(x, condPrior$location, condPrior$scale)
		}
		
		priorDensAtH0 = convolveFuns(fNorm, fCauchy, inverseH0)
	}
	

	#get the mean and condition effect posterior
	popMu = results$posteriors[[ paste(param, ".mu", sep="") ]]
	condEffect = results$posteriors[[ paste(param, "_cond[", cond, "]", sep="") ]]
	
	combinedPosterior = popMu + condEffect
	
	res = savageDickey(combinedPosterior, h0_val=inverseH0, priorDensAtH0 = priorDensAtH0)
	
	temp = data.frame(param = param, cond = cond, H0_value = H0_value, 
										bf01 = res$bf01, bf10 = res$bf10, success = res$success, stringsAsFactors = FALSE)
	
	temp
}


getSubsampleIterationsToRemove = function(totalIterations, subsamples, subsampleProportion) {
	if (subsamples < 1) {
		stop("You need to use at least one subsample.")
	}
	
	independentSubsamples = FALSE
	if (is.null(subsampleProportion)) {
		independentSubsamples = TRUE
		subsampleProportion = 1 / subsamples
	}
	
	subsampleProportion = min( max(subsampleProportion, 0), 1 )
	
	#if (subsamples > 1 && subsampleProportion == 1) {
	#	warning("The subsample proportion should be less than 1 if using multiple subsamples.")
	#}
	
	subsampleIterationsToRemove = list()
	
	if (independentSubsamples) {
		
		shuffledIterations = sample(1:totalIterations, totalIterations, replace=FALSE)
		
		for (sub in 1:subsamples) {
			
			iterationsToUse = floor(subsampleProportion * totalIterations)
			indicesToUse = ((sub - 1) * iterationsToUse + 1):(sub * iterationsToUse)
			subsampleIterationsToRemove[[sub]] = shuffledIterations[-indicesToUse]
			
		}
	} else {
		
		for (sub in 1:subsamples) {
			iterationsToRemove = round((1 - subsampleProportion) * totalIterations, 0)
			subsampleIterationsToRemove[[sub]] = sample(1:totalIterations, iterationsToRemove, replace = FALSE)
		}
		
	}
	
	subsampleIterationsToRemove
}






#' Test Differences Between Conditions
#' 
#' Tests the differences between the conditions in the experiment for different parameters. This test does all pairwise comparisons between all conditions. For omnibus tests, see \code{\link{testMainEffectsAndInteractions}}.
#' 
#' You can estimate how much the variability in the posterior chains affects the Bayes factors by using the \code{subsamples} and \code{subsampleProportion} arguments. If \code{subsamples} is greater than 1, multiple subsamples from the posterior chains will be taken and the standard deviation (and other measures) of the Bayes factors across the subsamples will be calculated. The number of iterations used in each subsample is a proportion of the total number of iterations and is set by \code{subsampleProportion}. Note that this is not the standard deviation of the Bayes factors over repeated samples of data sets, so it tells you nothing about what would happen if you had different data. It essentially tells you whether or not you ran enough iterations to have a stable Bayes factor estimate. The closer \code{subsampleProportion} is to 1, the less independent the subsamples will be, so you should use a reasonably low value of \code{subsampleProportion}. The degree to which the subsamples are independent influences to what extent the standard deviation is underestimated: The less independent, the larger the underestimate will be. If you want fully independent subsamples, you can set \code{subsampleProportion} to \code{NULL}. However, this means that the number of subsamples and the proportion of iterations in each subsample to be inversely related, which means that you have to choose between a low number of subsamples or a low number of iterations per subsample.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param A vector of basic parameter names for which to perform condition tests (e.g. "pMem"). If NULL (the default), tests are performed for all parameters with condition effects.
#' @param subsamples Number of subsamples of the posterior chains to take. If greater than 1, subsampleProportion should be set to a value between 0 and 1 (exclusive).
#' @param subsampleProportion The proportion of the total iterations to include in each subsample. This should probably only be less than 1 if \code{subsamples} is greater than 1. If \code{NULL}, \code{subsampleProportion} will be set to \code{1 / subsamples} and no iterations will be shared between subsamples (i.e. each subsample will be independent, except inasmuch as there is autocorrelation between iterations).
#' @param summarize Boolean. Should the results be summarized with \code{\link{summarizeSubsampleResults}}?
#' 
#' @return A data frame containing test results. It has the following columns:
#' \tabular{ll}{
#' 	\code{param} \tab The name of the parameter being tested.\cr
#' 	\code{cond} \tab The zero-indexed condition indices being compared.\cr
#' 	\code{bfType} \tab The type of Bayes factor in the "bf" column. "10" means that the Bayes factor is in favor of the alternative hypothesis that the two conditions differed. "01" means that the Bayes factor is in favor of the null hypothesis that the two conditions did not differ.\cr
#' 	\code{bf} \tab The Bayes factor. If using multiple subsamples, this is the mean Bayes factor. \cr
#' 	\code{sd} \tab Standard deviation of the Bayes factors. \cr
#' 	\code{min, median, max} \tab The minimum, median, and maximum of the Bayes factors. \cr
#' 	\code{p2.5, p97.5} \tab The 2.5 and 97.5 percentiles of the Bayes factors. 	
#' }
#'
#' @export
testConditionEffects = function(results, param = NULL, subsamples = 1, subsampleProportion = 1, summarize=TRUE) {

	#TODO: Rename to pairwiseComparisonsOfConditions?
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(results$config$iterations, subsamples, subsampleProportion)
	
	
	if (is.null(param)) {
		param = getParametersWithConditionEffects(results$config$conditionEffects)
	}
	
	condNames = results$config$factors$cond
	csName = results$config$cornerstoneConditionName
	nonCsNames = condNames[ condNames != csName ]

	
	pb = utils::txtProgressBar(0, 1, 0, style=3)
	currentStep = 1
	lastStep = length(subsampleIterationsToRemove) * length(param)
	
	
	allSubsamples = NULL
	for (sub in 1:length(subsampleIterationsToRemove)) {
		
		if (length(subsampleIterationsToRemove[[sub]]) > 0) {
			resultSubsample = removeBurnIn(results, subsampleIterationsToRemove[[sub]])
		} else {
			resultSubsample = results
		}
		
		for (p in param) {
		
			equalConds = getEqualConditionParameters(results, p)
			uniqueGroups = unique(equalConds$group)
			
			pairs = expand.grid(i = 1:length(uniqueGroups), j = 1:length(uniqueGroups) )
			pairs = pairs[ pairs$i < pairs$j, ]
			
			for (r in 1:nrow(pairs)) {
				i = pairs$i[r]
				j = pairs$j[r]
				
				g1conds = equalConds$cond[ equalConds$group == uniqueGroups[i] ]
				g2conds = equalConds$cond[ equalConds$group == uniqueGroups[j] ]
				
				prior1 = getConditionParameterPrior(results, p, g1conds[1])
				prior2 = getConditionParameterPrior(results, p, g2conds[1])
				
				post1 = resultSubsample$posteriors[[ paste(p, "_cond[", g1conds[1], "]", sep="") ]]
				post2 = resultSubsample$posteriors[[ paste(p, "_cond[", g2conds[1], "]", sep="") ]]
				
				priorDensAtH0 = stats::dcauchy(0, prior1$location - prior2$location, prior1$scale + prior2$scale)
				
				testRes = savageDickey(post1 - post2, h0_val=0, priorDensAtH0 = priorDensAtH0)
				
				condLabel = paste( paste(g1conds, collapse=","), " - ", paste(g2conds, collapse=","), sep="")
				
				temp = data.frame(param=p, cond=condLabel, 
													bf01=testRes$bf01, bf10=testRes$bf10, success = testRes$success, 
													stringsAsFactors=FALSE)
				
				allSubsamples = rbind(allSubsamples, temp)
			}
		
			utils::setTxtProgressBar(pb, value = currentStep / lastStep)
			currentStep = currentStep + 1
			
		}
		
	}
	
	close(pb)
	
	#TODO: Failure shouldn't be possible any more, so this could be removed.
	if (any(allSubsamples$success == FALSE)) {
		failed = allSubsamples[ allSubsamples$success == FALSE, ]
		failed = unique(failed)
		failed = failed[ , c("param", "cond") ]
		cat("Unable to test condition effects for the following parameters and conditions:")
		print(failed)
		cat("\n")
	}
	
	if (summarize) {
		rval = summarizeSubsampleResults(allSubsamples, aggregateBy = c("param", "cond"))
		rval = rval[ order(rval$param, rval$bfType, rval$cond), ]
	} else {
		rval = allSubsamples
		attr(rval, "aggregateColumns") = c("param", "cond")
	}
	
	rval
}



#' Calculate Whole Model WAIC
#' 
#' WAIC is a whole-model fit statistic, like AIC. However, WAIC cannot be directly compared with AIC, but it is conceptually very similar.
#' 
#' WAIC is an appropriate fit statistic for these models because
#' 1) it can be calculated without using posterior means, which some parameters do not have (specifically catMu and catActive). 
#' 2) it estimates the effective number of free parameters, which is wildly different from the actual number of free parameters.
#' See the documentation for \code{\link{calculateInappropriateFitStatistics}} for more information on the number of parameters.
#' 
#' There are two ways to estimate the effective number of free parameters for WAIC and results from both are reported.
#'  
#' WAIC is a whole-model fit statistic. These models have shared parameters that are shared by participants (the hierarchical parameters and condition effects). This means that the participants are not independent, so examining WAIC for individual participants is unprincipled and may give inaccurate results. Thus, you should leave \code{onlyTotal} at the default value of TRUE.
#' 
#' Note that you can use WAIC to compare model variants, like the between-item and within-item variants. You can also compare models that differ in other ways, such as which parameters have condition effects, the maximum number of categories, reduced models with some parameters set to constant values, more or less restrictive priors, etc.
#' 
#' You can estimate how much the variability in the posterior chains affects WAIC by using the \code{subsamples} and \code{subsampleProportion} arguments. If \code{subsamples} is greater than 1, multiple subsamples from the posterior chains will be taken and the standard deviation of WAIC (et al.) across the subsamples will be calculated. The number of iterations used in each subsample is a proportion of the total number of iterations and is set by \code{subsampleProportion}. Note that this is not the standard deviation of WAIC over repeated samples of data sets, so it tells you nothing about what would happen if you had different data. It essentially tells you whether or not you ran enough iterations to have a stable WAIC estimate. The closer \code{subsampleProportion} is to 1, the less independent the subsamples will be, so you should use a reasonably low value of \code{subsampleProportion}. The degree to which the subsamples are independent influences to what extent the standard deviation is underestimated: The less independent, the larger the underestimate will be. If you want fully independent subsamples, you can set \code{subsampleProportion} to NULL. However, this means that the number of subsamples and the proportion of iterations in each subsample to be inversely related, which means that you have to choose between a low number of subsamples or a low number of iterations per subsample.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param subsamples Number of subsamples of the posterior chains to take. If greater than 1, subsampleProportion should be set to a value between 0 and 1 (exclusive).
#' @param subsampleProportion The proportion of the total iterations to include in each subsample. This should probably only be less than 1 if \code{subsamples} is greater than 1. If \code{NULL}, \code{subsampleProportion} will be set to \code{1 / subsamples} and no iterations will be shared between subsamples (i.e. each subsample will be independent, except inasmuch as there is autocorrelation between iterations).
#' @param onlyTotal If \code{TRUE}, exclude participant-level WAIC values (which aren't really valid because of the fact that participants are not independent). I recommend you leave this at TRUE.
#' 
#' @return A data.frame containing WAIC values, estimates of the effective number of free parameters, and LPPD.
#' \tabular{ll}{
#'  \code{pnum} \tab If \code{onlyTotal == FALSE}, the participant related to the statistic. \cr
#' 	\code{stat} \tab The name of the fit statistic.\cr
#' 	\code{value} \tab The statistic value. If \code{subsamples > 1}, the mean statistic value across all subsamples. \cr
#' 	\code{sd} \tab The standard deviation of the statistic. \cr
#' 	\code{min, median, max} \tab The minimum, median, and maximum of the statistic. \cr
#' 	\code{p2.5, p97.5} \tab The 2.5 and 97.5 percentiles of the statistic. 	
#' }
#' 
#' @export
calculateWAIC = function(results, subsamples=1, subsampleProportion=1, onlyTotal=TRUE) {
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(results$config$iterations, subsamples, subsampleProportion)
	
	
	#??? This should not be necessary
	#results$config = verifyConfigurationList(config=results$config, data=results$data)
	
	#Convert catMu from degrees to radians if circular. This should really be done in C++
	if (results$config$dataType == "circular" && results$config$maxCategories > 0) {
		for (p in unique(results$data$pnum)) {
			for (i in 1:results$config$maxCategories) {
				cmName = paste("catMu[", p, ",", i-1, "]", sep="")
				results$posteriors[[ cmName ]] = CatContModel::d2r(results$posteriors[[ cmName ]])
			}
		}
	}
	
	allWAIC = NULL
	
	pb = utils::txtProgressBar(0, 1, 0, style=3)
	for (sub in 1:length(subsampleIterationsToRemove)) {
		if (length(subsampleIterationsToRemove[[sub]]) > 0) {
			noBurnIn = removeBurnIn(results, subsampleIterationsToRemove[[sub]])
		} else {
			noBurnIn = results
		}
		
		#if (length(subsampleIterationsToRemove) > 1) {
		#	cat(paste("Starting subsample ", sub, ".\n", sep=""))
		#}
		
		waic = CCM_CPP_calculateWAIC(noBurnIn)
		waic$subsample = sub
		
		allWAIC = rbind(allWAIC, waic)

		utils::setTxtProgressBar(pb, sub / length(subsampleIterationsToRemove))
	}
	close(pb)

	summaryWAIC = NULL
	
	for (p in unique(allWAIC$pnum)) {
		
		stats = c("WAIC_1", "WAIC_2", "P_1", "P_2", "LPPD")
		
		for (s in stats) {
			x = allWAIC[allWAIC$pnum == p, s]

			qs = as.numeric(stats::quantile(x, c(0, 0.025, 0.5, 0.975, 1)))
			
			if (subsamples > 1) {
				temp = data.frame(pnum = p, stat = s, value = mean(x), sd = stats::sd(x), 
													min=qs[1], p2.5=qs[2], median=qs[3], p97.5=qs[4], max=qs[5],
													stringsAsFactors=FALSE)
			} else {
				temp = data.frame(pnum = p, stat = s, value = mean(x), stringsAsFactors=FALSE)
			}
			
			summaryWAIC = rbind(summaryWAIC, temp)
			
		}
	}
	
	if (onlyTotal) {
		summaryWAIC = summaryWAIC[ summaryWAIC$pnum == "Total", ]
		summaryWAIC$pnum = NULL
	}
	
	summaryWAIC
}


#' Calculate Standard Fit Statistics that are Inappropriate for these Models
#' 
#' AIC and BIC are commonly used to compare models, but they are inappropriate for the
#' models fit by this package (with the possible exception of the ZL model). The reason
#' is that both AIC and BIC use the number of free parameters in the model as a penalty
#' term. However, the actual number of free parameters and the effective number of free
#' parameters are very different for these models. Part of the reason is that the 
#' catActive parameters don't provide very much flexibility. In addition, catMu does
#' nothing if the category is inactive, which is common. Finally, the hierarchical nature
#' of the model means that the population level parameters actually constrain the participant
#' level parameters, reducing the effective number of free parameters, but AIC and BIC count 
#' them as free parameters. All of this combined, the effective number of free parameters
#' is far less than the actual number of free parameters.
#' 
#' An appropriate fit statistc for these models is WAIC, which estimates the effective 
#' number of free parameters. See the \code{\link{calculateWAIC}} function of this package.
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
#' @param results The results from the \code{\link{runParameterEstimation}} function. Note that you must set \code{config$calculateParticipantLikelihoods} to \code{TRUE} in order to use this function.
#' @param onlyTotal If TRUE, exclude participant-level values (which aren't valid because of the fact that participants are not independent).
#' 
#' @return A data.frame containing several fit statistics.
#' 
#' @seealso \code{\link{calculateWAIC}} for an appropriate fit statistic.
#' 
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


