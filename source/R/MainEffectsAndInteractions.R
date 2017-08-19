
#' Posterior Distributions of Main Effect and Interaction Parameters
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param The name of a parameter with condition effect.
#' @param testedFactors A character vector giving the names of factors for which a hypothesis test could be performed. If there is only one factor, the main effect will be used. If there is more than one factor, the interaction of the factors will be the effect that is used. 
#' @param dmFactors See \code{\link[CMBBHT]{testHypothesis}}. Passed directly to \code{\link[CMBBHT]{getEffectParameters}}.
#' @param contrastType See \code{\link[CMBBHT]{testHypothesis}}. Passed directly to \code{\link[CMBBHT]{getEffectParameters}}.
#' 
#' @return A matrix with rows being iterations and columns being effect parameters. The columns are named with the following scheme: "F1.L1:F2.L2" where "Fn" is the name of a factor and "Ln" is the level of that factor and ":" indicates the combinations of factor levels.
#' 
#' @export
getEffectParameterPosteriors = function(results, param, testedFactors, dmFactors = testedFactors, contrastType = NULL) {
	
	postCMs = getPosteriorConditionEffects(results, param)
	
	gmeihtf = results$config$factors
	gmeihtf$cond = NULL #NO EXTRA COLUMNS
	
	CMBBHT::getEffectParameters(cellMeans=postCMs, factors=gmeihtf, testedFactors = testedFactors, dmFactors = dmFactors, contrastType = contrastType)
	
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
		testedFactors = getFactorsForConditionEffect(results, param)
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









# You could add an argument: useFullDMForUnbalanced If TRUE and the design is unbalanced, the design matrix that is used for all tests will be the design matrix with all effects in it. This means that a main effect will not really be a marginal test, because interactions will be accounted for. This may or may not be what you want to do.
testMEI_singleParameter.Old = function(results, param, priorSamples = NULL, doPairwise = FALSE, testFunction = NULL) {

	factorsToTest = getFactorsForConditionEffect(results, param)
	if (length(factorsToTest) == 0) {
		return(NULL)
	}

	combinationLayers = getMultilevelFactorCrossing(factorsToTest)
	
	allTests = NULL
	for (i in 1:nrow(combinationLayers)) {
		
		theseFactors = strsplit(combinationLayers$combination[i], ":", fixed=TRUE)[[1]]
		
		factorName = paste0(theseFactors, collapse=":")
		
		thisRes = testSingleEffect(results, param, testedFactors = theseFactors, 
															 priorSamples = priorSamples, 
															 testFunction = testFunction)
		if (!thisRes$success) {
			thisRes$bf10 = thisRes$bf01 = NA
		}
		
		omnibus = data.frame(param=param, factor=factorName, levels="Omnibus", 
												 bf10=thisRes$bf10, bf01=thisRes$bf01, success=thisRes$success)
		allTests = rbind(allTests, omnibus)
		
		
		if (doPairwise && combinationLayers$layer[i] == 1) {
			#do pairwise comps
			
			thisFactor = theseFactors #rename just to be clear that there is only 1 factor
			
			allLevels = unique(results$config$factors[ , thisFactor ])
			comb = t(utils::combn(allLevels, 2))
			for (j in 1:nrow(comb)) {
				
				usedFactorLevels = list()
				usedFactorLevels[[thisFactor]] = comb[j,]
				usedFactorLevels = as.data.frame(usedFactorLevels, stringsAsFactors = FALSE)
				
				thisRes = testSingleEffect(results, param, testedFactors = thisFactor, 
																	 usedFactorLevels = usedFactorLevels, priorSamples = priorSamples,
																	 testFunction = testFunction)
				if (!thisRes$success) {
					thisRes$bf10 = thisRes$bf01 = NA
				}
				
				levelNames = paste0(comb[j,], collapse=", ")
				
				pairwise = data.frame(param=param, factor=factorName, levels=levelNames, 
															bf10=thisRes$bf10, bf01=thisRes$bf01, success=thisRes$success)
				allTests = rbind(allTests, pairwise)
				
			} #j
		} #if doPairwise
	} #i

	allTests
	
}




testMainEffectsAndInteractions.Old = function(results, param=NULL, 
																					subsamples = 50, subsampleProportion = 1, 
																					summarize = TRUE, doPairwise = FALSE, 
																					testFunction = CMBBHT::testFunction_SDDR) 
{
	
	if (is.null(results$config$factors)) {
		stop('You must provide "results$config$factors".')
	}
	
	gmeihtf = results$config$factors
	gmeihtf$cond = NULL
	if (!CMBBHT::isDesignFullyCrossed(gmeihtf)) {
		warning("Design is not fully crossed, which means that main effects and interactions cannot be orthogonal. This complicates the meaning of tests. See the documentation of this function for more information.")
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
															priorSamples=priorSamples, doPairwise=doPairwise, 
															testFunction = testFunction)
			BFs = rbind(BFs, result)
			
			utils::setTxtProgressBar(pb, value = currentStep / lastStep)
			currentStep = currentStep + 1
			
		}
	}
	
	close(pb)
	
	rval = cleanAndSummarizeMEIResults(BFs, summarize=summarize, aggregateBy=c("param", "factor", "levels"))
	rval
	
}






