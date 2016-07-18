
#postVect: vector of the posterior distribution of the parameter you want to test
#h0_val: The value of the parameter tested in the null hypothesis.
#priorDensAtH0: The density of the prior distribution at h0_val
savageDickey = function(postVect, h0_val, priorDensAtH0) {

	ls = polspline::logspline(postVect)
	
	postDens = polspline::dlogspline(h0_val, ls)
	
	bf10 = priorDensAtH0 / postDens
	
	list(bf01=1/bf10, bf10=bf10)
}




#' Test the Amount of Categorical Responding Present in the Data
#' 
#' Test the amount of categorical responding present in the data in each condition.
#' To do this test, you must have run the "betweenItem" model variant. Then, use this function on the results from that run to test whether categorical responding is or is not present in the data. A lack of categorical responding is the same as fully continuous responding, which happens when pContBetween = 1 and pCatGuess = 0. As explained in the article, testing 1 and 0 is not possible, so values near 1 and 0 should be tested instead. The test is performed in each condition separately.
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
		for (cond in results$conditions$levels) {
			
			res = testMeanParameterValue(results, param = param, cond = cond, H0_value = H0_test[[param]])
			#temp = data.frame(param = param, cond = cond, H0_value = H0_test[[param]], bf01 = res$bf01, bf10 = res$bf10)
			
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
	
	transformInverse = getParameterTransformation(param, results, inverse=TRUE)
	
	inverseH0 = transformInverse(H0_value)
	
	popMu0 = results$priors[[ paste(param, ".mu.mu", sep="") ]]
	popSd0 = sqrt( results$priors[[ paste(param, ".mu.var", sep="") ]] )
	
	condLoc0 = results$priors[[ paste(param, "_cond.loc", sep="") ]]
	condScale0 = results$priors[[ paste(param, "_cond.scale", sep="") ]]
	
	priorDensAtH0 = 0
	if (cond == results$config$cornerstoneConditionName) {
		
		priorDensAtH0 = stats::dnorm(inverseH0, popMu0, popSd0)
		
	} else {
		
		fNorm = function(x) {
			stats::dnorm(x, popMu0, popSd0)
		}
		fCauchy = function(x) {
			stats::dcauchy(x, condLoc0, condScale0)
		}
		
		priorDensAtH0 = convolveFuns(fNorm, fCauchy, inverseH0)
	}
	

	#get the mean and condition effect
	popMu = results$posteriors[[ paste(param, ".mu", sep="") ]]
	condEffect = results$posteriors[[ paste(param, "_cond[", cond, "]", sep="") ]]
	
	combined = popMu + condEffect
	
	res = savageDickey(combined, h0_val=inverseH0, priorDensAtH0 = priorDensAtH0)
	
	temp = data.frame(param = param, cond = cond, H0_value = H0_value, bf01 = res$bf01, bf10 = res$bf10, stringsAsFactors = FALSE)
	
	temp
}



#' Test Condition Effects
#' 
#' Tests the differences between the conditions in the experiment for different parameters. This test does all pairwise comparisons between all conditions.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param A vector of basic parameter names for which to perform condition tests (e.g. "pMem"). If NULL (the default), tests are performed for all parameters with condition effects.
#' 
#' @return A data frame containing test results. It has the following columns:
#' \tabular{ll}{
#' 	\code{param} \tab The name of the parameter being tested.\cr
#' 	\code{cond} \tab The zero-indexed condition indices being compared.\cr
#' 	\code{bf01} \tab The Bayes factor in favor of the null hypothesis that the two conditions did not differ.\cr
#' 	\code{bf10} \tab The Bayes factor in favor of the alternative hypothesis that the two conditions differed.\cr
#' }
#'
#' @export
#'
testConditionEffects = function(results, param = NULL) {

	if (is.null(param)) {
		param = results$config$parametersWithConditionEffects
	}
	
	condNames = results$conditions$levels
	csName = results$config$cornerstoneConditionName
	nonCsNames = condNames[ condNames != csName ]
	
	rval = NULL
	
	for (p in param) {
		
		scale = results$priors[[ paste(p, "_cond.scale", sep="") ]]
		
		nonCsPrior = function(x) { stats::dcauchy(x, 0, scale) }
		difPrior = function(x) { stats::dcauchy(x, 0, scale * 2) }
		

		#do each non-cs condition on its own (i.e. comparisons to cs cond)
		for (i in 1:length(nonCsNames)) {
			x = results$posteriors[[ paste(p, "_cond[", nonCsNames[i], "]", sep="") ]]
			
			sdr = savageDickey(x, h0_val=0, priorDensAtH0 = nonCsPrior(0))
			
			temp = data.frame(param=p, cond=paste(csName, " - ", nonCsNames[i], sep=""),
												bf01=sdr$bf01, bf10=sdr$bf10)
			rval = rbind(rval, temp)
		}
		
		#do each pairwise comparison of non-cs conditions
		for(i in 1:length(nonCsNames)) {
			for(j in 1:length(nonCsNames)) {
				if (i < j) {
					
					xi = results$posteriors[[ paste(p, "_cond[", nonCsNames[i], "]", sep="") ]]
					xj = results$posteriors[[ paste(p, "_cond[", nonCsNames[j], "]", sep="") ]]
					
					sdr = savageDickey(xi - xj, h0_val=0, priorDensAtH0 = difPrior(0))
					
					temp = data.frame(param=p, cond=paste(nonCsNames[i], " - ", nonCsNames[j], sep=""), 
														bf01=sdr$bf01, bf10=sdr$bf10)
					rval = rbind(rval, temp)
				}
			}
		}
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
#' WAIC is a whole-model fit statistic. These models have shared parameters that are shared by participants (the hierarchical parameters and condition effects). This means that the participants are not independent, so examining WAIC for individual participants is unprincipled and may give inaccurate results. Thus, you should leave onlyTotal at the default value of TRUE.
#' 
#' Note that you can use WAIC to compare model variants, like the betweenItem and withinItem variants. You can also compare models that differ in other ways, such as which parameters have condition effects, the maximum number of categories, etc.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param onlyTotal If TRUE, exclude participant-level WAIC values (which aren't valid because of the fact that participants are not independent).
#' 
#' @return A data.frame containing WAIC values, estimates of the effective number of free parameters, and LPPD (basically log likelihood).
#' 
#' @export
calculateWAIC = function(results, onlyTotal=TRUE) {
	
	#Convert catMu from degrees to radians
	if (results$config$dataType == "circular" && results$config$maxCategories > 0) {
		for (p in unique(results$data$pnum)) {
			for (i in 1:results$config$maxCategories) {
				cmName = paste("catMu[", p, ",", i-1, "]", sep="")
				results$posteriors[[ cmName ]] = CatContModel::d2r(results$posteriors[[ cmName ]])
			}
		}
	}
	
	res = CCM_CPP_calculateWAIC(results)
	
	if (onlyTotal) {
		res = res[ res$pnum == "Total", ]
	}
	res
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
#' Thus, for each fit statistic, both mean and standard deviation are reported.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
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


