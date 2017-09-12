



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

convolveFuns = function(f, g, t, range=c(-Inf, Inf)) {
	
	#all of the ways to get to t.
	h = function(tau) {
		f(tau) * g(t - tau)
	}
	
	dens = stats::integrate(h, lower=range[1], upper=range[2])
	
	dens$value
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
#' @seealso [`testMeanParameterValue`] for a more general way to do these kinds of hypothesis tests.
#' 
#' @md
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
testConditionEffects.Old = function(results, param = NULL, subsamples = 1, subsampleProportion = 1, summarize=TRUE) {

	#TODO: Rename to pairwiseComparisonsOfConditions?
	#or conditionPairwiseComparisons?
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(results$config$iterations, subsamples, subsampleProportion)
	
	
	if (is.null(param)) {
		param = getParametersWithConditionEffects(results$config$conditionEffects)
	}
	
	pb = utils::txtProgressBar(0, 1, 0, style=3)
	currentStep = 1
	lastStep = length(subsampleIterationsToRemove) * length(param)
	
	
	allSubsamples = NULL
	for (sub in 1:length(subsampleIterationsToRemove)) {
		
		if (length(subsampleIterationsToRemove[[sub]]) > 0) {
			resSub = removeBurnIn(results, subsampleIterationsToRemove[[sub]])
		} else {
			resSub = results
		}
		
		for (p in param) {
		
			equalConds = getEqualConditionParameters(resSub, p)
			uniqueGroups = unique(equalConds$group)
			
			pairs = expand.grid(i = 1:length(uniqueGroups), j = 1:length(uniqueGroups) )
			pairs = pairs[ pairs$i < pairs$j, ]
			
			for (r in 1:nrow(pairs)) {

				g1conds = equalConds$cond[ equalConds$group == uniqueGroups[ pairs$i[r] ] ]
				g2conds = equalConds$cond[ equalConds$group == uniqueGroups[ pairs$j[r] ] ]
				
				prior1 = getConditionParameterPrior(resSub, p, g1conds[1])
				prior2 = getConditionParameterPrior(resSub, p, g2conds[1])
				
				post1 = resSub$posteriors[[ paste0(p, "_cond[", g1conds[1], "]") ]]
				post2 = resSub$posteriors[[ paste0(p, "_cond[", g2conds[1], "]") ]]
				
				priorDensAtH0 = stats::dcauchy(0, prior1$location - prior2$location, prior1$scale + prior2$scale)
				
				testRes = savageDickey(post1 - post2, h0_val=0, priorDensAtH0 = priorDensAtH0)
				
				condLabel = paste0( paste(g1conds, collapse=","), " - ", paste(g2conds, collapse=","))
				
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
	
	rval = cleanAndSummarizeMEIResults(allSubsamples, summarize = summarize, aggregateBy = c("param", "cond"))
	rval
}




