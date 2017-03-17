



getPriorConditionEffects = function(results, param, priorSamples) {
	
	factors = results$config$factors
	
	priorDs = matrix(NA, nrow=priorSamples, ncol=nrow(factors))
	
	#first pass: Sample from free parameters and the cornerstone condition
	for (i in 1:nrow(factors)) {
		
		target = paste(param, "_cond[", factors$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, param, factors$cond[i])
		
		#If this is its own source (i.e. it is an unconstrained parameter), sample from it if it is free
		if (target == rootSource) {
			prior = getConditionParameterPrior(results, param, factors$cond[i])
			if (prior$scale == 0) {
				priorDs[,i] = prior$location
			} else {
				priorDs[,i] = stats::rcauchy(priorSamples, prior$location, prior$scale)
			}
		}
		
	}
	
	#second pass: Copy equality constrained parameters
	#you can't use different samples for equality constrained parameters
	#because then they would onlt be equal in distribution, not value.
	for (i in 1:nrow(factors)) {
		
		target = paste(param, "_cond[", factors$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, param, factors$cond[i])
		
		if (target != rootSource) {
			
			parts = getConditionParameterParts(rootSource)
			
			sourceCol = which(factors$cond == parts$cond)
			priorDs[,i] = priorDs[,sourceCol]
			
		}
	}
	
	colnames(priorDs) = factors$cond
	
	priorDs
	
}

getPosteriorConditionEffects = function(results, param) {
	
	factors = results$config$factors
	
	postDs = matrix(NA, nrow=results$config$iterations, ncol=nrow(factors))
	for (i in 1:nrow(factors)) {
		postDs[,i] = results$posteriors[[ paste(param, "_cond[", factors$cond[i], "]", sep="") ]]
	}
	
	colnames(postDs) = factors$cond
	
	postDs
}

#' Matrices of Prior and Posterior Condition Effects
#' 
#' Get matrices of prior and posterior condition effects for a given parameter. 
#' These can be used to perform tests with, e.g., the CMBBHT package.
#' 
#' @param results The results from the [`runParameterEstimation`] function.
#' @param param The parameter to get condition effects for.
#' @param priorSamples The number of samples to take from the priors on the condition effects.
#' 
#' @return A list with two matrices, `prior` and `post`. Each column is one condition effect and each row is one sample from the prior or posterior.
#' 
#' @md
#' @export
getConditionEffects = function(results, param, priorSamples=results$config$iterations) {
	prior = getPriorConditionEffects(results, param, priorSamples)
	post = getPosteriorConditionEffects(results, param)
	list(prior=prior, post=post)
}

#' Posterior Distributions of Main Effect and Interaction Parameters
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param The name of a parameter with condition effect.
#' @param testedFactors A character vector giving the names of factors for which a hypothesis test could be performed. If there is only one factor, the main effect will be used. If there is more than one factor, the interaction of the factors will be the effect that is used. 
#' @param dmFactors Character vector. The factors to use to construct the design matrix. For a fully-crossed (balanced) design, this can always be equal to \code{testFactors} (the default). For non-fully-crossed designs, you may sometimes want to create a design matrix using some factors, but perform a hypothesis test with only some of those factors (\code{testedFactors} must be a subset of \code{dmFactors}).
#' @param contrastType Character (or function). The contrast to use to create the design matrix. Can be any of the function names on the documentation page for \code{contr.sum}. For a non-fully-crossed (unbalanced) design, you should use either "contr.treatment" or "contr.SAS". For a balanced design, you can use anything, but psychologists are most used to "contr.sum", which uses sums-to-zero constraints.
#' 
#' @return A matrix with row being iterations and columns being effect parameters. The columns are named with the following scheme: "F1.L1:F2.L2" where "Fn" is the name of a factor and "Ln" is the level of that factor.
#' 
#' @export
getEffectParameterPosteriors = function(results, param, testedFactors, dmFactors = testedFactors, contrastType = NULL) {
	
	postCMs = getPosteriorConditionEffects(results, param)
	
	gmeihtf = results$config$factors
	gmeihtf$cond = NULL #NO EXTRA COLUMNS
	
	CMBBHT::getEffectParameters(cellMeans=postCMs, factors=gmeihtf, testedFactors = testedFactors, dmFactors = dmFactors, contrastType = contrastType)

}


#' Perform a Single Hypothesis Test of a Main Effect or Interaction
#' 
#' This function performs a single hypothesis test of a selected effect.
#' This is a lower-level function than \code{\link{testMainEffectsAndInteractions}} and it gives you more control than that function. For even more control, see [`getConditionEffects`] and the CMBBHT package (https://github.com/hardmanko/CMBBHT).
#' 
#' See the details of \code{\link{testMainEffectsAndInteractions}} for some discussion of fully-crossed vs non-fully-crossed designs.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param The name of the parameter for which to perform the test.
#' @param testedFactors Character vector. The factors for which to perform the hypothesis test as a vector of factor names. A single factor name results in the test of the main effect of the factor. Multiple factor names result in the test of the interaction of all of those factors.
#' @param dmFactors Character vector or formula. The factors to use to construct the design matrix. For a fully-crossed (balanced) design, this can always be equal to \code{testFactors} (the default). For non-fully-crossed designs, you may sometimes want to create a design matrix using some factors, but perform a hypothesis test with only some of those factors (\code{testedFactors} must be a subset of \code{dmFactors}). You may instead supply a \code{formula} like that taken by \code{\link{model.matrix}} which will be used to create the design matrix.
#' @param usedFactorLevels A \code{data.frame} with columns for each of the factors in \code{testedFactors}. Each row specifies factor levels that should be included in the test. This allows you to do things like pairwise comparisons of specific factor levels.
#' @param priorSamples Number of samples to take from the prior distribution of the effect parameters. You should not change this from the default unless you are using a custom \code{testFunction}, in which case you might want to use a different value.
#' @param testFunction See \code{\link[CMBBHT]{testHypothesis}}.
#' @param contrastType See \code{\link[CMBBHT]{testHypothesis}}.
#' 
#' @return Depends on the choice of \code{testFunction}. Typically a list with, at least, elements \code{bf10} and \code{bf01} which are the Bayes factors in favor of and against the tested effect, respectively.
#' 
#' @export
testSingleEffect = function(results, param, testedFactors, dmFactors = testedFactors, 
														usedFactorLevels = NULL, priorSamples = results$config$iterations, 
														testFunction = NULL, contrastType = NULL) 
{
	if (is.null(testFunction)) {
		testFunction = CMBBHT::testFunction_SDDR
	}
	
	gmeihtf = results$config$factors
	gmeihtf$cond = NULL
	
	ces = getConditionEffects(results, param, priorSamples)
	
	CMBBHT::testHypothesis(ces$prior, ces$post, gmeihtf, testedFactors, 
												 dmFactors=dmFactors, contrastType=contrastType,
												 testFunction = testFunction, usedFactorLevels = usedFactorLevels)
	
}

# You could add an argument: useFullDMForUnbalanced If TRUE and the design is unbalanced, the design matrix that is used for all tests will be the design matrix with all effects in it. This means that a main effect will not really be a marginal test, because interactions will be accounted for. This may or may not be what you want to do.
testMEI_singleParameter = function(results, param, priorSamples = NULL, doPairwise = FALSE, testFunction = NULL) {

	factorsToTest = getFactorsForConditionEffect(results$config, param)
	if (length(factorsToTest) == 0) {
		return(NULL)
	}

	#This is basically C/P from another function. Maybe wrap it up?
	combinationLayers = NULL
	
	for (layer in 1:length(factorsToTest)) {
		r = utils::combn(factorsToTest, layer)
		r = t(r)
		for (j in 1:nrow(r)) {
			x = sort(r[j,])
			y = paste(x, collapse=":")
			
			combinationLayers = rbind(combinationLayers, data.frame(combination=y, layer=length(x), stringsAsFactors=FALSE))
		}
	}
	
	allTests = NULL
	for (i in 1:nrow(combinationLayers)) {
		
		theseFactors = strsplit(combinationLayers$combination[i], ":", fixed=TRUE)[[1]]
		
		factorName = paste0(theseFactors, collapse=":")
		
		thisRes = testSingleEffect(results, param, testedFactors = theseFactors, 
															 priorSamples = priorSamples, 
															 testFunction = testFunction)
		
		omnibus = data.frame(param=param, factor=factorName, levels="Omnibus", 
												 bf10=thisRes$bf10, bf01=thisRes$bf01, success=thisRes$success)
		allTests = rbind(allTests, omnibus)
		
		
		if (doPairwise && combinationLayers$layer[i] == 1) {
			#do pairwise comps
			
			thisFactor = theseFactors #just to be clear that there is only 1 factor
			
			allLevels = unique(results$config$factors[ , thisFactor ])
			comb = t(utils::combn(allLevels, 2))
			for (j in 1:nrow(comb)) {
				
				uniqueFL = list()
				uniqueFL[[thisFactor]] = comb[j,]
				uniqueFL = as.data.frame(uniqueFL, stringsAsFactors = FALSE)
				
				thisRes = testSingleEffect(results, param, testedFactors = thisFactor, 
																	 usedFactorLevels = uniqueFL, priorSamples = priorSamples,
																	 testFunction = testFunction)
				
				levelNames = paste0(comb[j,], collapse=", ")
				
				pairwise = data.frame(param=param, factor=factorName, levels=levelNames, 
															bf10=thisRes$bf10, bf01=thisRes$bf01, success=thisRes$success)
				allTests = rbind(allTests, pairwise)
				
			} #j
		} #if doPairwise
	} #i
	
	allTests = allTests[ order(allTests$levels, allTests$factor), ]
	
	allTests
	
}



#' Test Main Effects and Interactions of Factors
#' 
#' Perform hypothesis tests of main effects and interactions for one-factor and multi-factor designs. If your design is fully crossed, your life is simple. If your design is not fully crossed, you may not be able to use this function for all of your tests (see Details).
#' 
#' You must provide a `data.frame` containing the mapping from conditions to factor levels. This should be provided in `results$config$factors`. See [`runParameterEstimation`] for more information about creating this `data.frame`. If you are using a one-factor design, this will have been created for you and you don't need to do anything. If using multiple factors, you should have given `config$factors` to `runParameterEstimation`.
#' 
#' This function uses kernel density estimation to estimate the densities of some relevant quantities. This procedure is somewhat noisy. As such, I recommend that you perform the procedure many times, the number of which can be configured with the `subsamples` argument. Then, aggregate results from the many repetitions of the procedure can be analyzed, which is done by default but can be changed by setting `summarize` to `FALSE`.
#' I recommend using many `subsamples` to see how noisy the estimation is. You can leave `subsampleProportion` at 1 or use a somewhat lower value. I would recommend against using a value of `subsampleProportion` that would result in fewer than 1,000 iterations being used per subsample.
#' 
#' For designs that are fully-crossed, sums-to-zero contrasts are used by default. Like any other kind of orthogonal contrast, sums-to-zero contrasts result in main effects and interactions that are independent of one another. Thus, for example, a main effect is the same regardless of whether you also included an interaction in the design or not. For designs that are not fully crossed, treatment contrasts are used by default. Treatment contrasts are non-orthogonal, which means that effects are not independent of one another. For example, a main effect might change depending on whether or not you include an interaction in the model. Thus, when working with non-fully-crossed designs, you must decide what effects you want to include in the model when you are testing an effect. 
#' This function does marginal tests: It only estimates what it needs to to do the test. Thus, if testing a main effect, only that main effect is estimated. If testing a two-factor interaction, only the two related main effects and that interaction are estimated.
#' You cannot do more this with this function, but see [`testSingleEffect`] for a function that gives you more control over the test that is performed. 
#' In addition, for designs that are not fully croseed, you can still use [`testConditionEffects`] to examine pairwise comparisons. See the introduction.pdf manual for more discussion of non-fully-crossed designs.
#' 
#' @param results The results from the [`runParameterEstimation`] function.
#' @param param Optional. Character vector of names of parameters to perform tests for. If NULL (default), is set to all parameters with condition effects.
#' @param summarize If TRUE (default), the results across subsamples will be summarized. If FALSE, the results from each of the subsamples will be returned. Those results can be later summarized with [`summarizeSubsampleResults`].
#' @param subsamples Number of subsamples of the posterior chains to take. If greater than 1, subsampleProportion should be set to a value between 0 and 1 (exclusive).
#' @param subsampleProportion The proportion of the total iterations to include in each subsample. This should probably only be less than 1 if `subsamples` is greater than 1. If `NULL`, `subsampleProportion` will be set to `1 / subsamples` and no iterations will be shared between subsamples (i.e. each subsample will be independent, except inasmuch as there is autocorrelation between iterations).
#' @param doPairwise Do pairwise tests of differences between levels of main effects (these are often called "post-hoc" tests).
#' 
#' @md
#' @export
testMainEffectsAndInteractions = function(results, param=NULL, 
																					subsamples = 50, subsampleProportion = 1, 
																					summarize=TRUE,	doPairwise = FALSE) 
{
	
	if (is.null(results$config$factors)) {
		stop('You must provide "results$config$factors".')
	}
	
	gmeihtf = results$config$factors
	gmeihtf$cond = NULL
	if (!CMBBHT::isDesignFullyCrossed(gmeihtf)) {
		warning("Design is not fully crossed, which means that main effects and interactions cannot be orthogonal. See the documentation of this function for more information.")
	}
	
	if (is.null(param)) {
		param = getParametersWithConditionEffects(results$config$conditionEffects)
	}
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(results$config$iterations, subsamples, subsampleProportion)
	
	BFs = NULL

	pb = utils::txtProgressBar(0, 1, 0, style=3)
	currentStep = 1
	lastStep = length(subsampleIterationsToRemove) * length(param)
	
	for (sub in 1:length(subsampleIterationsToRemove)) {
		
		if (length(subsampleIterationsToRemove[[sub]]) > 0) {
			resultSubsample = removeBurnIn(results, subsampleIterationsToRemove[[sub]])
		} else {
			resultSubsample = results
		}
		
		#Very important: The kernel density estimation procedure has a problem.
		#For the way in which it is used, the density depends on the number of
		#samples. More samples results in less density at the tested point (in a tail) 
		#of the kinds of distributions used. Thus, the prior and posterior sample counts must match.
		priorSamples = resultSubsample$config$iterations

		
		for (pInd in 1:length(param)) {
			
			result = testMEI_singleParameter(results=resultSubsample, param=param[pInd], 
															priorSamples=priorSamples, doPairwise=doPairwise)
			BFs = rbind(BFs, result)
			
			utils::setTxtProgressBar(pb, value = currentStep / lastStep)
			currentStep = currentStep + 1
			
		}
	}
	
	close(pb)
	
	if (any(BFs$success == FALSE)) {
		
		ff = stats::aggregate(success ~ param * factor * levels, BFs, function(x) { sum(!x) })
		ff$failures = ff$success
		ff$success = NULL

		warning( paste0("Bayes factor estimation failed for ", sum(BFs$success == FALSE), " subsamples. The cases with failures have been printed to the console. Failures have been stripped from the results.") )
		cat("\n\nFailures:\n")
		print( ff[ ff$failures > 0, ] )
		cat("\n\n")
		
		BFs = BFs[ BFs$success, ]
	}
	
	if (summarize) {
		rval = summarizeSubsampleResults(BFs, aggregateBy = c("param", "factor", "levels"))
	
		rval$notOmnibus = rval$levels != "Omnibus"
		rval = rval[ order(rval$notOmnibus, rval$param, rval$factor, rval$levels), ]
		
		rval$notOmnibus = NULL
	} else {
		rval = BFs
		attr(rval, "aggregateColumns") = c("param", "factor", "levels")
	}
	
	rval
	
}

#' Summarize Results from Multiple Subsamples
#' 
#' This function should be used with the values returned by \code{\link{testMainEffectsAndInteractions}} and \code{\link{testConditionEffects}} when the \code{summarize} argument is \code{FALSE}. It summarizes Bayes factors across many repeated estimates of those Bayes factors.
#' 
#' @param BFs A data.frame containing the individual Bayes factors. It should have a format like the result of \code{\link{testMainEffectsAndInteractions}} or \code{\link{testConditionEffects}}. 
#' @param proportioniles Percentiles divided by 100 to calculate.
#' @param geometricZs A numeric vector of Z-values. Geometric BF quantiles will be calculated based on the geometric mean and standard deviation for each of these provided Z-values.
#' @param consistencyCutoff A numeric vector of cutoffs. The proportion of Bayes factors above each cutoff is calculated. Defaults to \code{c(1, 3, 10)}.
#' @param logBF Summarize log Bayes factors? If \code{FALSE}, no logs will be taken. If \code{TRUE}, log base 10 BFs will be used. If a numeric value, that value will be used as the base for the logarithm.
#' @param aggregateBy Columns of \code{BFs} to aggregate by. Typically not required as it is read from an attribute of \code{BFs} called \code{aggregateColumns}.
#' 
#' @return A data frame containing summarized test results. It has the following columns:
#' \tabular{ll}{
#' 	\code{...} \tab The columns that were aggregated by. See the \code{aggregateBy} argument. \cr
#' 	\code{bfType} \tab The type of Bayes factor in this row. "10" means that the Bayes factor is in favor of the alternative hypothesis that there is an effect (e.g. two conditions differed; there is a main effect). "01" means that the Bayes factor is in favor of the null hypothesis that there is no effect.\cr
#' 	\code{bf} \tab The arithmetic mean Bayes factor. \cr
#' 	\code{sd} \tab Standard deviation of the Bayes factors. \cr
#' 	\code{geo.mean} \tab The geometric mean Bayes factor. Bayes factors estimated with the approach used in this package tend to vary exponentially, which makes the geometric mean a possibly better measure than the arithmetic mean. Note that the geometric mean and the median tend to be in closer agreement than the arithmetic mean and the median. \cr
#' 	\code{geo.sd} \tab Geometric standard deviation of the Bayes factors. Multiply the geo.mean by the geo.sd to go up one standard deviation. In general, \code{geo.mean * geo.sd^z} will give you the geometric value corresponding to the given z score. \cr
#' 	\code{geo.mean + n SD} \tab The geometric mean "plus" \code{n} standard deviations, where the \code{n} values are given by the \code{geometricZs} argument. \cr
#' 	`p(BF > n)` \tab The proportion of Bayes factors greater than `n`. See `consistencyCutoff` to set `n`. \cr
#' 	\code{Min, Median, Max} \tab The minimum, median, and maximum of the Bayes factors. \cr
#' 	\code{n\%} \tab Other percentiles, as given in the \code{proportioniles} argument.
#' }
#' 
#' @md
#' @export
summarizeSubsampleResults = function(BFs, proportioniles = c(0, 0.025, 0.5, 0.975, 1), 
															geometricZs = NULL, consistencyCutoff = c(1, 3, 10), logBF = FALSE, 
															aggregateBy = NULL) {
	

	if (logBF == TRUE) {
		logBF = 10
		if (consistencyCutoff == 1) {
			consistencyCutoff = 0
		}
	}
	
	if (is.null(aggregateBy)) {
		aggregateBy = attr(BFs, "aggregateColumns")
		if (is.null(aggregateBy)) {
			stop("No columns to aggregate by.")
		}
	}
	
	uniqueBF = unique(BFs[ , aggregateBy ])
	
	rval = NULL
	for (i in 1:nrow(uniqueBF)) {

		for (bf in c("bf01", "bf10")) {
			
			theseAgg = list()
			useRows = rep(TRUE, nrow(BFs))
			for (agg in aggregateBy) {
				theseAgg[[agg]] = uniqueBF[ i, agg ]
				
				useRows = useRows & (BFs[ , agg ] == uniqueBF[ i, agg ])
			}

			dfl = theseAgg
			dfl[["bfType"]] = substr(bf, 3, 4)
			
			x = BFs[ useRows, bf ]
			if (logBF != FALSE) {
				x = log(x, base=logBF)
			}
			
			if (length(x) > 1) {

				dfl[["bf"]] = mean(x)
				dfl[["sd"]] = stats::sd(x)
				

				dfl[["geo.mean"]] = geoMean(x)
				dfl[["geo.sd"]] = geoSD(x)
				if (logBF != FALSE) {
					dfl[["geo.mean"]] = NA
					dfl[["geo.sd"]] = NA
				}
				if (length(geometricZs) > 0) {
					for (j in 1:length(geometricZs)) {
						gq = geoQ(geometricZs[j], dfl[["geo.mean"]], dfl[["geo.sd"]])
						n = paste0("geo.mean ", ifelse( geometricZs[j] > 0, "+ ", " "), geometricZs[j], " sd")
						dfl[[ n ]] = gq
					}
				}
				
				for (j in 1:length(consistencyCutoff)) {
					name = paste0("p(BF > ", consistencyCutoff[j], ")")
					dfl[[name]] = mean(x > consistencyCutoff[j])
				}
				

				#percentiles
				if (length(proportioniles) > 0) {
					qs = stats::quantile(x, proportioniles)
					if (proportioniles[1] == 0) {
						names(qs) = c("Min", names(qs)[-1])
					}
					if (any(proportioniles == 0.5)) {
						ind = which(proportioniles == 0.5)
						n = names(qs)
						n[ind] = "Median"
						names(qs) = n
					}
					if (proportioniles[length(proportioniles)] == 1) {
						names(qs) = c(names(qs)[-length(proportioniles)], "Max")
					}
					for (j in 1:length(qs)) {
						dfl[[ names(qs)[j] ]] = as.numeric(qs[j])
					}
				}
				
			} else {
				
				dfl[["bf"]] = x

			}
			
			if (length(names(rval)) > length(names(dfl))) {
				for (n in names(rval)) {
					if (!(n %in% names(dfl))) {
						dfl[[n]] = NA
					}
				}
			}
			
			dfl = as.data.frame(dfl, stringsAsFactors=FALSE, check.names=FALSE)
			
			rval = rbind(rval, dfl)
		}
		
	}
	
	rval
}


#all(x > 0)
geoMean = function(x) {
	if (any(x <= 0)) {
		stop("All x must be > 0 to calculate the geometric mean.")
	}
	#prod(x)^(1 / length(x)) 
	#which translates to
	exp( 1/length(x) * sum(log(x)) )
}

geoSD = function(x) {
	m = geoMean(x)
	
	exp(sqrt(sum(log(x / m)^2) / length(x)))
}

geoZ = function(x, mu, sigma) {
	log(x / mu, base=sigma)
}

geoQ = function(z, mu, sigma) {
	mu * sigma ^ z
}


#' Marginal Priors on Main Effect and Interaction Parameters
#' 
#' Given the structure of the model, it is possible to specify the priors on condition effects.
#' The priors on the resulting main effect and interaction (MEI) effect parameters, however, are 
#' not directly specified, but can be calculated, which is what this function does.
#' 
#' The priors on condition effects are Cauchy distributions. The individual MEI effect parameters
#' are calculated by taking a linear combination (weighted sum) of the condition effect parameters.
#' 
#' In particular, let Y = sum(X * W), where Y is the resulting Cauchy distribution, X is the vector
#' of Cauchy distributions to be combined, and W are weights. A linear combination of Cauchy 
#' distributions is a Cauchy distribution with properties discussed below.
#' 
#' Let L and S be vectors of Locations and Scales of the Cauchy distributions in X. 
#' Then the location and scale parameters of Y are 
#' 
#' \code{L_Y = sum(L * W)}
#' \code{S_Y = sum(S * abs(W))}
#' 
#' I haven't found a citation for this anywhere, but have confirmed it to several significant digits in simulations.
#' 
#' The results of these calculations depend on a lot of information, which is most easily 
#' provided in the results of parameter estimation. To examine the effects of changing the
#' priors, you can test "new" priors with the \code{priorLoc} and \code{priorScale} arguments.
#' Note, however, that the typical proscription on nonzero prior locations holds here as well,
#' in the sense that it works mathematically, but nonzero prior locations are kind of bizarre.
#' 
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function. Note that you can run 1 iteration and still have all the information you need.
#' @param param The name of the parameter for which to calculate MEI effect parameter priors.
#' @param testedFactors Character vector. The factors for which to perform the hypothesis test as a vector of factor names. A single factor name results in the test of the main effect of the factor. Multiple factor names result in the test of the interaction of all of those factors.
#' @param dmFactors Character vector. The factors to use to construct the design matrix. For a fully-crossed (balanced) design, this can always be equal to \code{testFactors} (the default). For non-fully-crossed designs, you may sometimes want to create a design matrix using some factors, but perform a hypothesis test with only some of those factors (\code{testedFactors} must be a subset of \code{dmFactors}).
#' @param contrastType Character (or function). The contrast to use to create the design matrix. Can be any of the function names on the documentation page for \code{contr.sum}. For a non-fully-crossed (unbalanced) design, you should use either "contr.treatment" or "contr.SAS". For a balanced design, you can use anything, but psychologists are most used to "contr.sum", which uses sums-to-zero constraints.
#' @param priorLoc A new prior location to try, overriding the value in results$priors.
#' @param priorScale A new prior scale to try, overriding the value in results$priors.
#' 
#' @return A \code{data.frame} with four columns: the \code{factor} being used, the \code{effect} parameter, the prior \code{location}, and the prior \code{scale}.
#' 
#' @export
calculateMarginalEffectParameterPriors = function(results, param, testedFactors = NULL, dmFactors = NULL, contrastType = NULL, priorLoc = NULL, priorScale = NULL) {
	
	gmeihtf = results$config$factors
	gmeihtf$cond = NULL
	
	if (is.null(contrastType)) {
		if (CMBBHT::isDesignFullyCrossed(gmeihtf, FALSE)) {
			contrastType = "contr.sum"
		} else {
			contrastType = "contr.treatment"
		}
	}
	
	factors = results$config$factors
	
	testedFactorsSpecified = !is.null(testedFactors)
	if (is.null(testedFactors)) {
		testedFactors = getFactorsForConditionEffect(results$config, param)
		if (length(testedFactors) == 0) {
			return(NULL)
		}
	}
	
	calculateList = list()
	if (testedFactorsSpecified) {
		calculateList[[ paste0(fNames, collapse = ":") ]] = testedFactors
	} else {
		for (layer in 1:length(testedFactors)) {
			comb = utils::combn(testedFactors, layer)
			
			for (i in 1:ncol(comb)) {
				fNames = comb[,i]
				overallEffect = paste0(fNames, collapse = ":")
				
				calculateList[[ overallEffect ]] = fNames
			}
		}
	}
	
	
	priors = NULL
	
	for (n in names(calculateList)) {
		
		fNames = calculateList[[n]]
		
		uniqueFL = unique(subset(factors, select = fNames))
		
		locations = rep(NA, nrow(factors))
		scales = locations
		for (i in 1:nrow(factors)) {
			r = getConditionParameterPrior(results, param, factors$cond[i])
			locations[i] = r$location
			scales[i] = r$scale
			if (factors$cond[i] != results$config$cornerstoneConditionName) {
				if (!is.null(priorLoc)) {
					locations[i] = priorLoc
				}
				if (!is.null(priorScale)) {
					scales[i] = priorScale
				}
			}
		}
		
		if (is.null(dmFactors)) {
			thisDmFactors = fNames
		} else {
			thisDmFactors = dmFactors
		}

		S_s = CMBBHT::getPartialFilledS(gmeihtf, testedFactors = fNames, dmFactors = thisDmFactors, contrastType = contrastType)
		
		for (j in 1:ncol(S_s)) {
			res = cauchyRvLinearCombination(locations, scales, S_s[,j])
			
			temp = data.frame(factor = paste0(fNames, collapse = ":"), effect = colnames(S_s)[j], 
												location = res$location, scale = res$scale)
			priors = rbind(priors, temp)
		}
		
	}
	
	priors
}

