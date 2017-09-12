
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
#' `L_Y = sum(L * W)`
#' `S_Y = sum(S * abs(W))`
#' 
#' I haven't found a citation for this anywhere, but have confirmed it to several significant digits in simulations.
#' 
#' The results of these calculations depend on a lot of information, which is most easily 
#' provided in the results of parameter estimation. To examine the effects of changing the
#' priors, you can test "new" priors with the `priorLoc` and `priorScale` arguments.
#' Note, however, that the typical proscription on nonzero prior locations holds here as well,
#' in the sense that it works mathematically, but nonzero prior locations are kind of bizarre
#' for condition effect parameters.
#' 
#' 
#' @param results The results from the [`runParameterEstimation`] function. Note that you can run only 1 iteration and still have all the information you need.
#' @param param The name of the parameter for which to calculate MEI effect parameter priors.
#' @param testedFactors Character vector. The factors for which to perform the hypothesis test as a vector of factor names. A single factor name results in the test of the main effect of the factor. Multiple factor names result in the test of the interaction of all of those factors.
#' @param dmFactors Character vector. The factors to use to construct the design matrix. For a fully-crossed (balanced) design, this can always be equal to `testFactors` (the default). For non-fully-crossed designs, you may sometimes want to create a design matrix using some factors, but perform a hypothesis test with only some of those factors (`testedFactors` must be a subset of `dmFactors`).
#' @param contrastType Character (or function). The contrast to use to create the design matrix. Can be any of the function names on the documentation page for `contr.sum`. For a non-fully-crossed (unbalanced) design, you should use either "contr.treatment" or "contr.SAS". For a balanced design, you can use anything, but psychologists are most used to "contr.sum", which uses sums-to-zero constraints.
#' @param priorLoc A new prior location to try, overriding the value in results$priors.
#' @param priorScale A new prior scale to try, overriding the value in results$priors.
#' 
#' @return A `data.frame` with four columns: the `factor` being used, the `effect` parameter, the prior `location`, and the prior `scale`.
#' 
#' @family WP functions
#' @md
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




