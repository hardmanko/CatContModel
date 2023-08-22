

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
#' To do this test, you must have run one of the "betweenItem" or "withinItem" model variants. 
#' Then, use this function on the results from that run to test whether categorical responding is or is not present in the data. 
#' A lack of categorical responding is the same as fully continuous responding, which happens when `pContBetween = 1` and `pCatGuess = 0`. 
#' As explained in Hardman, Vergauwe, and Ricker (2017), testing 1 and 0 is not possible, so values near 1 and 0 (0.99 or 0.01) should be tested instead. 
#' The test is performed in each condition separately.
#' 
#' @param results Results from [`runParameterEstimation`].
#' @param pCatGuess_test Values of pCatGuess to test.
#' @param pContBetween_test Values of pContBetween to test.
#' @param pContWithin_test Values of pContWithin to test.
#' 
#' @return A `data.frame` with the results of the tests. 
#' 
#' @seealso [`testMeanParameterValue`] for a more general way to do these kinds of hypothesis tests. [`posteriorMeansAndCredibleIntervals`] for means and credible intervals.
#' 
#' @family test functions
#' @family WP functions
#' 
#' @export
testCategoricalResponding = function(results, pCatGuess_test = c(0.01, 0.99), pContBetween_test = c(0.01, 0.99), pContWithin_test = c(0.01, 0.99)) {

	if (!resultIsType(results, "WP")) {
		stop("testCategoricalResponding only accepts WP results objects. See the Glossary (listed in the package functions index).")
	}
	
	H0_test = list(pCatGuess = pCatGuess_test)
	if (results$config$modelVariant == "betweenItem") {
	  H0_test$pContBetween = pContBetween_test
	  
	} else if (results$config$modelVariant == "withinItem") {
	  H0_test$pContWithin = pContWithin_test
	  
	} else if (results$config$modelVariant == "betweenAndWithin") {
	  H0_test$pContBetween = pContBetween_test
	  H0_test$pContWithin = pContWithin_test
	  
	} else {
	  stop("Invalid modelVariant.")
	}
	
	rval = NULL
	for (parName in names(H0_test)) {
	  for (testVal in H0_test[[parName]]) {
  		for (cond in results$config$factors$cond) {
  			
  			testRes = testMeanParameterValue(results, parName = parName, cond = cond, H0_value = testVal)
  
  			rval = rbind(rval, testRes)
  		}
	  }
	}
	
	rval
}

#' Test the Mean Value for a Parameter in a Condition
#' 
#' Perform the test of whether the mean of the parameter is equal to some value. 
#' This test is performed in a condition.
#' 
#' @param results Results from [`runParameterEstimation`].
#' @param parName The name of the parameter to test.
#' @param cond The condition in which to test the value.
#' @param H0_value The null hypothesized value of the parameter. This should be provided in the manifest space (e.g., if testing a probability parameter, this should be between 0 and 1).
#' 
#' @family test functions
#' @family WP functions
#' 
#' @export
#' 
#' @examples \dontrun{
#' testMeanParameterValue(results, parName = "pMem", cond = "1", H0_value = 0.95)
#' }
#' 
testMeanParameterValue = function(results, parName, cond, H0_value) {
	
	if (!resultIsType(results, "WP")) {
		stop("testMeanParameterValue only accepts WP results objects. See the Glossary (listed in the package index).")
	}
	
	transformInverse = getParameterTransformation(results, parName, inverse=TRUE)
	
	inverseH0 = transformInverse(H0_value)
	
	popMu0 = results$priors[[ paste(parName, "_part.mu.mu", sep="") ]]
	popSd0 = sqrt( results$priors[[ paste(parName, "_part.mu.var", sep="") ]] )
	
	condPrior = getConditionParameterPrior(results, parName, cond)
	
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
	popMu = results$posteriors[[ paste(parName, "_part.mu", sep="") ]]
	condEffect = results$posteriors[[ paste(parName, "_cond[", cond, "]", sep="") ]]
	
	combinedPosterior = popMu + condEffect
	
	res = savageDickey(combinedPosterior, h0_val=inverseH0, priorDensAtH0 = priorDensAtH0)
	
	# Option: Use a prior sampling approach instead?
	#CMBBHT::valueTest_SDDR()
	
	temp = data.frame(parName = parName, cond = cond, H0_value = H0_value, 
										bf01 = res$bf01, bf10 = res$bf10, success = res$success, stringsAsFactors = FALSE)
	
	temp
}


