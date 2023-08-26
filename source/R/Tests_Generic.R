# TODO Renaming:
# testSingleEffect -> HT_Factorial_Single
# testMainEffectsAndInteractions -> HT_Factorial
# testConditionEffects -> HT_ConditionPairs, HT_Pairwise
#
# Also rename this file to HypothesisTests.R
# And Tests_Generic_Utilities.R -> HypothesisTests_Utilities.R


#' Perform a Single Hypothesis Test of a Main Effect or Interaction
#' 
#' This function performs a single hypothesis test of a selected effect.
#' This is a lower-level function than [`testMainEffectsAndInteractions`] and it gives you more control than that function. 
#' For even more control, see [`getConditionEffects`] and the CMBBHT package (https://github.com/hardmanko/CMBBHT).
#' 
#' See the details of [`testMainEffectsAndInteractions`] for some discussion of fully-crossed vs non-fully-crossed designs.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param parName The name of the parameter for which to perform the test.
#' @param testedFactors Character vector. The factors for which to perform the hypothesis test as a vector of factor names. A single factor name results in the test of the main effect of the factor. Multiple factor names result in the test of the interaction of all of those factors.
#' @param dmFactors See \code{\link[CMBBHT]{testHypothesis}}.
#' @param usedFactorLevels See \code{\link[CMBBHT]{testHypothesis}}.
#' @param priorSamples Number of samples to take from the prior distribution of the effect parameters. You should not change this from the default unless you are non-default `testFunction` (like a curried version of \code{\link[CMBBHT]{testFunction_EPA}}), in which case you might want to use a different value.
#' @param testFunction See \code{\link[CMBBHT]{testHypothesis}}.
#' @param contrastType See \code{\link[CMBBHT]{testHypothesis}}.
#' @param addMu If `TRUE`, the grand mean of the parameter is added to the condition effects. If `FALSE`, nothing is added. If `NULL` (default), is equal to `manifest` for WP designs. `addMu` is forced to `TRUE` for BP designs to account for mean differences between groups.
#' @param manifest If `TRUE`, the resulting parameter will be in the manifest space. If `FALSE`, the resulting parameter will be in the latent space. (See [`Glossary`].) If `NULL`, it is set to `FALSE` for probability parameters and `TRUE` for SD parameters.
#' 
#' @return Depends on the choice of `testFunction`. Typically a list with, at least, elements `bf10` and `bf01` which are the Bayes factors in favor of and against the tested effect, respectively.
#' 
#' @family test functions
#' @family generic functions
#' 
#' @export
testSingleEffect = function(res, parName, testedFactors, dmFactors = testedFactors, 
														usedFactorLevels = NULL, priorSamples = res$runConfig$iterations, 
														testFunction = CMBBHT::testFunction_SDDR, contrastType = NULL,
														addMu = NULL, manifest = NULL) 
{

  if (resultIsType(res, "Parallel")) {
    stop("Cannot analyze Parallel results.")
  }
  
	priorSamples = valueIfNull(priorSamples, res$runConfig$iterations)
	
	
	if (is.null(manifest)) {
		if (parName %in% getParamNames(types="sd")) {
			manifest = TRUE
		} else {
			manifest = FALSE
		}
	}
	
	if (resultIsType(res, "WP")) {
		addMu = valueIfNull(addMu, manifest)
	} else if (resultIsType(res, "BP")) {
		addMu = TRUE
	}
	
	factors = normalizeFactors(res$config$factors)
	
	cefs = getConditionEffects(res, parName, priorSamples = priorSamples, addMu = addMu, manifest = manifest)
	
	if (!all(factors$key %in% cefs$colKeys$key)) {
		stop("factors contains a key not found in the condition effects.")
	}
	
	#Order correctly and select only needed cells
	cefs$prior = subset(cefs$prior, select=factors$key)
	cefs$post = subset(cefs$post, select=factors$key)
	
	#Select only real factors
	cmbbhtf = subset(factors, select = getAllFactorNames(factors))
	
	res = tryCatch({
		CMBBHT::testHypothesis(cefs$prior, cefs$post, cmbbhtf, testedFactors, 
													 dmFactors = dmFactors, contrastType = contrastType,
													 testFunction = testFunction, usedFactorLevels = usedFactorLevels)
	}, error = function(e) {
		logWarning(e$message)
		return(list(success=FALSE))
	})
	
	if (!res$success) {
		res$bf10 = res$bf01 = NA
	}
	
	res
}



testMEI_singleParameter = function(res, parName, priorSamples = res$runConfig$iterations, doPairwise = FALSE, 
																					 testFunction = CMBBHT::testFunction_SDDR, 
																					 addMu = NULL, manifest = NULL) 
{
	factors = normalizeFactors(res$config$factors)
	
	#Strip WP factors without condition effects, but use all others.
	factWithCondEff = getFactorsForConditionEffect(res, parName)
	FN = getFactorTypeToName(factors)
	FN$wp = FN$wp[ FN$wp %in% factWithCondEff ]
	
	factorsToTest = c(FN$bp, FN$wp)
	
	if (length(factorsToTest) == 0) {
		return(NULL)
	}
	
	combinationLayers = getMultilevelFactorCrossing(factorsToTest)
	
	allTests = NULL
	for (i in 1:nrow(combinationLayers)) {
		
		theseFactors = strsplit(combinationLayers$combination[i], ":", fixed=TRUE)[[1]]
		
		factorName = paste0(theseFactors, collapse=":")
		
		thisRes = testSingleEffect(res, parName, testedFactors = theseFactors, priorSamples = priorSamples, 
															 testFunction = testFunction, addMu = addMu, manifest = manifest)
		
		omnibus = data.frame(parName=parName, factor=factorName, levels="Omnibus", 
												 bf10=thisRes$bf10, bf01=thisRes$bf01, success=thisRes$success, 
												 stringsAsFactors = FALSE)
		allTests = rbind(allTests, omnibus)
		
		
		if (doPairwise && combinationLayers$layer[i] == 1) {
			
			thisFactor = theseFactors #rename just to be clear that there is only 1 factor
			
			allLevels = unique(factors[ , thisFactor ])
			comb = t(utils::combn(allLevels, 2))
			for (j in 1:nrow(comb)) {
				
				usedFactorLevels = list()
				usedFactorLevels[[thisFactor]] = comb[j,]
				usedFactorLevels = as.data.frame(usedFactorLevels, stringsAsFactors = FALSE)
				
				thisRes = testSingleEffect(res, parName, testedFactors = thisFactor, 
																	 usedFactorLevels = usedFactorLevels, priorSamples = priorSamples, 
																	 testFunction = testFunction, addMu = addMu, manifest = manifest)
				
				levelNames = paste0(comb[j,], collapse=", ")
				
				pairwise = data.frame(parName=parName, factor=factorName, levels=levelNames, 
															bf10=thisRes$bf10, bf01=thisRes$bf01, success=thisRes$success, 
															stringsAsFactors = FALSE)
				allTests = rbind(allTests, pairwise)
				
			} #j
		} #if doPairwise
	} #i
	
	allTests
}

#' Test Main Effects and Interactions of Factors
#' 
#' Perform hypothesis tests of main effects and interactions for one-factor and multi-factor designs. 
#' If your design is fully crossed, your life is simple. 
#' If your design is not fully crossed, you may not be able to use this function for all of your tests (see Details).
#' 
#' @details
#' You must provide a `data.frame` containing the mapping from conditions to factor levels. 
#' This should be provided in `res$config$factors`. See [`makeFactors`] for more information about creating this `data.frame`. 
#' If you are using a one-factor design, this will have been created for you and you don't need to do anything. 
#' If using multiple factors, you should have given `config$factors` to `runParameterEstimation`.
#' 
#' This function uses kernel density estimation to estimate the densities of some relevant quantities. 
#' This procedure is somewhat noisy. 
#' As such, I recommend that you perform the procedure many times, the number of which can be configured with the `subsamples` argument. 
#' Then, aggregate results from the many repetitions of the procedure can be analyzed, which is done by default but can be changed by setting `summarize` to `FALSE`.
#' I recommend using many `subsamples` to see how noisy the estimation is. 
#' You can leave `subsampleProportion` at 1 or use a somewhat lower value. 
#' I would recommend against using a value of `subsampleProportion` that would result in fewer than 1,000 iterations being used per subsample.
#' 
#' For designs that are fully-crossed, sums-to-zero contrasts are used by default. 
#' Like any other kind of orthogonal contrast, sums-to-zero contrasts result in main effects and interactions that are independent of one another. 
#' Thus, for example, a main effect is the same regardless of whether you also included an interaction in the design or not. 
#' For designs that are not fully crossed, treatment contrasts are used by default. 
#' Treatment contrasts are non-orthogonal, which means that effects are not independent of one another. 
#' For example, a main effect might change depending on whether or not you include an interaction in the model. 
#' Thus, when working with non-fully-crossed designs, you must decide what effects you want to include in the model when you are testing an effect. 
#' 
#' This function does marginal tests: It only estimates what it needs to do the test. 
#' Thus, if testing a main effect, only that main effect is modeled. If testing a two-factor interaction, 
#' only the two related main effects and that interaction are modeled. 
#' This is different than testing for the existence of, e.g., a main effect in the context of an interaction.
#' You cannot do more with this function, but see [`testSingleEffect`] for a function that gives you more control over the test that is performed. 
#' In addition, for designs that are not fully croseed, you can still use [`testConditionEffects`] to examine pairwise comparisons. 
#' See the introduction.pdf manual for more discussion of non-fully-crossed designs.
#' 
#' 
#' @section testMainEffectsAndInteractions.BP:
#' `testMainEffectsAndInteractions.BP` has the additional `manifest` argument. If `TRUE`, manifest
#' parameter values are used. If `FALSE`, latent parameter values are used. If `NULL` (default), latent
#' parameter values are used for the probability parameters and manifest parameter values are used for standard deviation parameters.
#' 
#' These tests are not totally ideal because tests are based on dummy parameters that are the sum of the 
#' grand mean of the participant-level parameters (the hierarchical mean parameter, e.g. `pMem.mu`) and the condition effect parameters (e.g. `pMem_cond[1]`). 
#' These dummy parameters are created for the purposes of the test, but are not to be found anywhere in the models. 
#' Thus, the values used by the test are a very good but slightly inexact approximation. 
#' For more discussion of this issue, see the the Binomial Tutorial part of the documentation for 
#' CMBBHT: https://github.com/hardmanko/CMBBHT/releases/download/v0.1.3/BinomialTutorial.pdf 
#' Focus on the Between-Participants Design section beginning on page 15, but earlier material provides important context.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param parNames Optional. Character vector of names of parameters to perform tests for. If `NULL`, is set to all parameters with condition effects.
#' @param summarize If `TRUE`, the results across subsamples will be summarized. If `FALSE`, the results from each of the subsamples will be returned. Those results can be later summarized with [`summarizeSubsampleResults`].
#' @param subsamples Number of subsamples of the posterior chains to take. If greater than 1, `subsampleProportion` should be set to a value between 0 and 1 (exclusive).
#' @param subsampleProportion The proportion of the total iterations to include in each subsample. This should probably only be less than 1 if `subsamples` is greater than 1. If `NULL`, `subsampleProportion` will be set to `1 / subsamples` and no iterations will be shared between subsamples (i.e. each subsample will be independent, except inasmuch as there is autocorrelation between iterations).
#' @param doPairwise Do pairwise tests of differences between levels of main effects (these are often called "post-hoc" tests).
#' @param testFunction See \code{\link[CMBBHT]{testHypothesis}} in the CMBBHT package.
#' @param addMu Passed to same argument of [`getConditionEffects`].
#' @param manifest Passed to same argument of [`getConditionEffects`].
#' @param progress Boolean. If `TRUE`, a text progress bar is shown.
#' 
#' @return Depends on the value of `summarize`. If `summarize == TRUE`, it will have the same return value as [`summarizeSubsampleResults`], so see that function. If `summarize == FALSE`, a `data.frame` with columns
#' * `parName`: The parameter name.
#' * `factor`: The factor name.
#' * `levels`: The levels of the factor. If "Omnibus", all levels are used in an omnibus test of the effect. If all of the tests are omnibus tests (which is true as long as `doPairwise == FALSE`), this column is omitted for simplicity.
#' * `bf10`: The Bayes factor in favor of there being a difference between factor levels.
#' * `bf01`: The Bayes factor against there being a difference between factor levels.
#' * `success`: If `TRUE`, the Bayes factor was estimated successfully. If `FALSE`, it was not.
#' 
#' 
#' @family test functions
#' @family generic functions
#' 
#' @export
testMainEffectsAndInteractions = function(res, parNames = NULL, 
																						 subsamples = 1, subsampleProportion = 1, 
																						 summarize = TRUE, doPairwise = FALSE, 
																						 testFunction = CMBBHT::testFunction_SDDR,
																						 addMu = NULL, manifest = NULL, progress = subsamples > 1) 
{
  if (resultIsType(res, "Parallel")) {
    stop("Cannot analyze Parallel results.")
  }
  
  factors = normalizeFactors(res$config$factors)
	
	cmbbhtf = subset(factors, select=getAllFactorNames(factors))
	if (!CMBBHT::isDesignFullyCrossed(cmbbhtf)) {
		logWarning("Design is not fully crossed, which means that main effects and interactions cannot be orthogonal. This complicates the meaning of tests. See the documentation of testMainEffectsAndInteractions for more information.")
	}
	
	if (is.null(parNames)) {
		#All parameters can vary because its BP. 
		parNames = getParamNames(res$config$modelVariant, types=c("prob", "sd"))
	  
	  # Param without condition effects get NULL rval from testMEI_singleParameter
		#param = getParametersWithConditionEffects(res$config$conditionEffects)
	}
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(res$runConfig$iterations, subsamples, subsampleProportion)
	#subsamples = length(subsampleIterationsToRemove)
	
	BFs = NULL
	
	if (progress) {
		pb = utils::txtProgressBar(0, 1, 0, style=3)
		
		currentStep = 1
		lastStep = subsamples * length(parNames)
	}

	for (sub in 1:subsamples) {
		
		if (length(subsampleIterationsToRemove[[sub]]) > 0) {
			resSub = removeBurnIn(res, subsampleIterationsToRemove[[sub]])
		} else {
			resSub = res
		}
		
		for (pn in parNames) {
			
			htr = testMEI_singleParameter(resSub, parName = pn,
																		priorSamples = resSub$runConfig$iterations, doPairwise = doPairwise, 
																		testFunction = testFunction, addMu = addMu, manifest = manifest)
			BFs = rbind(BFs, htr)
			
			if (progress) {
				utils::setTxtProgressBar(pb, value = currentStep / lastStep)
				currentStep = currentStep + 1
			}
			
		}
	}
	
	if (progress) {
		close(pb)
	}
	
	aggregateBy = c("parName", "factor", "levels")
	if (all(BFs$levels == "Omnibus")) {
		BFs$levels = NULL
		aggregateBy = aggregateBy[ aggregateBy != "levels" ]
	}
	
	rval = cleanAndSummarizeMEIResults(BFs, summarize=summarize, aggregateBy = aggregateBy)
	rval
}





#' Test Differences Between Conditions
#' 
#' Tests the differences between the conditions in the experiment for different parameters. 
#' This test does all pairwise comparisons between all conditions. For omnibus tests, see [`testMainEffectsAndInteractions`].
#' 
#' @details
#' 
#' You can estimate how much the variability in the posterior chains affects the Bayes factors by using the 
#' `subsamples` and `subsampleProportion` arguments. If `subsamples` is greater than 1, 
#' multiple subsamples from the posterior chains will be taken and the standard deviation (and other measures) 
#' of the Bayes factors across the subsamples will be calculated. The number of iterations used in each subsample 
#' is a proportion of the total number of iterations and is set by `subsampleProportion`. 
#' Note that this is not the standard deviation of the Bayes factors over repeated samples of data sets, 
#' so it tells you nothing about what would happen if you had different data. 
#' It essentially tells you whether or not you ran enough iterations to have a stable Bayes factor estimate. 
#' The closer `subsampleProportion` is to 1, the less independent the subsamples will be, 
#' so you should use a reasonably low value of `subsampleProportion`. 
#' The degree to which the subsamples are independent influences to what extent the standard deviation is underestimated: 
#' The less independent, the larger the underestimate will be. 
#' If you want fully independent subsamples, you can set `subsampleProportion` to `NULL`. 
#' However, this means that the number of subsamples and the proportion of iterations in each subsample to be 
#' inversely related, which means that you have to choose between a low number of subsamples or a low number of iterations per subsample.
#' 
#' Note that the return value is different if subsamples are used (`subsamples > 1`) and the results 
#' are summarized (`summarize = TRUE`). See [`summarizeSubsampleResults`] for the return value format in that case.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param parNames A vector of parameter names for which to perform condition tests (e.g. "pMem"). If `NULL` (the default), tests are performed for all parameters with condition effects.
#' @param addMu Passed to same argument of [`getConditionEffects`]. If using a Between-Participants design, `addMu` must be `TRUE`.
#' @param manifest Passed to same argument of [`getConditionEffects`].
#' @param credP Credible interval proportion for the posterior difference between the tested conditions. Ignored if `subsamples != 1`.
#' @param subsamples Number of subsamples of the posterior chains to take. See details. If greater than 1, subsampleProportion should be set to a value between 0 and 1 (exclusive).
#' @param subsampleProportion The proportion of the total iterations to include in each subsample. This should probably only be less than 1 if \code{subsamples} is greater than 1. If `NULL`, `subsampleProportion` will be set to `1 / subsamples` and no iterations will be shared between subsamples (i.e. each subsample will be independent, except inasmuch as there is autocorrelation between iterations).
#' @param summarize Boolean. Should the results be summarized with \code{\link{summarizeSubsampleResults}}?
#' @param progress Boolean. If `TRUE`, a text progress bar is shown.
#' 
#' @return A data frame containing test results with the following columns (see details for info about different return value if using subsamples and summarizing):
#' \tabular{ll}{
#' 	\code{parName} \tab The name of the parameter being tested.\cr
#' 	\code{key} \tab The conditions being compared, like cond1 - cond2, where "-" can be interpreted as a minus sign.\cr
#' 	\code{bf01} \tab Bayes factor in favor of the null hypothesis of no difference between the conditions.\cr
#' 	\code{bf10} \tab Bayes factor in favor of the alternative hypothesis of there being a difference between the conditions. Note: `bf10 * bf01 = 1`.\cr
#' 	\code{success} \tab `TRUE` if the Bayes factor was calculated successfully.\cr
#'  \code{difMean, difLower, difUpper} \tab Based on `credP` and the posterior distribution of the difference between conditions. The lower and upper ends of the credible interval (`ciLower` and `ciUpper`) and the mean difference (`postMean`).\cr
#' }
#'
#' @family test functions
#' @family generic functions
#'
#' @rdname testConditionEffects
#' @export
testConditionEffects = function(res, parNames = NULL, credP = 0.95, addMu = TRUE, manifest = TRUE, subsamples = 1, subsampleProportion = 1, summarize = FALSE, progress = subsamples > 1) 
{

  if (resultIsType(res, "BP") && !addMu) {
    stop("For Between-Participants designs, addMu must be TRUE.")
  } else if (resultIsType(res, "Parallel")) {
    stop("Cannot analyze Parallel results.")
  }
  
  if (!is.null(credP)) {
    if (credP <= 0 || credP >= 1) {
    	logWarning("credP is outside of the interval (0,1) and will be ignored.")
      credP = NULL
    }
  }
  
  if (subsamples > 1) {
    if (!is.null(credP)) {
    	logWarning("credP is ignored if subsamples != 1.")
      credP = NULL
    }
    if (!summarize) {
      logMsg("Note: if using subsamples, you probably want to set summarize to TRUE.")
    }
  } else if (subsamples == 1) {
    summarize = FALSE
  } else {
    stop("Invalid value for the subsamples argument.")
  }
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(getCompletedIterations(res), subsamples, subsampleProportion)
	
	
	if (is.null(parNames)) {
		#param = getParamNames(res$config$modelVariant, types=c("prob", "sd"))
		parNames = getParametersWithConditionEffects(res$config$conditionEffects)
	}
	
	if (progress) {
		pb = utils::txtProgressBar(0, 1, 0, style=3)
		currentStep = 1
		lastStep = length(subsampleIterationsToRemove) * length(parNames)
	}
	
	allSubsamples = NULL
	for (sub in 1:length(subsampleIterationsToRemove)) {
		
		if (length(subsampleIterationsToRemove[[sub]]) > 0) {
			resSub = removeBurnIn(res, subsampleIterationsToRemove[[sub]])
		} else {
			resSub = res
		}
		
		if (resultIsType(resSub, "WP")) {
			groups = list()
			groups[[ defaultGroupName() ]] = resSub
		} else if (resultIsType(resSub, "BP")) {
			groups = resSub$groups
		}
		
		for (pn in parNames) {
			
			equalConds = list()
			unique_groupEq = NULL
			
			for (grp in names(groups)) {
				
				eqCond = getEqualConditionParameters(groups[[grp]], pn)
				eqCond$group = grp
				
				eqCond$key = paste(grp, eqCond$cond, sep=":")
				
				equalConds[[ grp ]] = eqCond
				unique_eq = unique(eqCond$equalityGroup)
				
				unique_groupEq = rbind(unique_groupEq, data.frame(group = grp, equalityGroup = unique_eq, stringsAsFactors = FALSE))
			}
			
			pairs = expand.grid(i = 1:nrow(unique_groupEq), j = 1:nrow(unique_groupEq) )
			pairs = pairs[ pairs$i < pairs$j, ]
			
			if (nrow(pairs) >= 1) {
				
				condEff = getConditionEffects(resSub, parName = pn, priorSamples = resSub$runConfig$iterations, addMu = addMu, manifest = manifest)
				
				for (r in 1:nrow(pairs)) {
					
					i = pairs$i[r]
					j = pairs$j[r]
					
					eqCond_i = equalConds[[ unique_groupEq$group[ i ] ]]
					eqCond_j = equalConds[[ unique_groupEq$group[ j ] ]]
					
					eqKeys_i = eqCond_i[ eqCond_i$equalityGroup == unique_groupEq$equalityGroup[ i ], "key" ]
					eqKeys_j = eqCond_j[ eqCond_j$equalityGroup == unique_groupEq$equalityGroup[ j ], "key" ]
					
					# Use first key because all keys are in the same equality group
					difPrior = condEff$prior[ , eqKeys_i[1] ] - condEff$prior[ , eqKeys_j[1] ]
					difPost = condEff$post[ , eqKeys_i[1] ] - condEff$post[ , eqKeys_j[1] ]
					
					ht = CMBBHT::valueTest_SDDR(difPrior, difPost, testVal = 0)
					
					keyLabel = paste0( paste(eqKeys_i, collapse="/"), " - ", paste(eqKeys_j, collapse="/"))
					
					temp = data.frame(parName=pn, key=keyLabel, 
														bf01=ht$bf01, bf10=ht$bf10, success = ht$success, 
														stringsAsFactors=FALSE)
					
					if (!is.null(credP)) {
					  credLowerPercentile = (1 - credP) / 2
					  
					  difPostQs = stats::quantile(difPost, c(credLowerPercentile, 1 - credLowerPercentile))
					  names(difPostQs) = NULL
					  
					  temp$difMean = mean(difPost)
					  temp$difLower = difPostQs[1]
					  temp$difUpper = difPostQs[2]
					  
					  #temp$bfReject0 = temp$bf10 > 3
					  #temp$ciReject0 = temp$difLower > 0 | temp$difUpper < 0
					}
					
					allSubsamples = rbind(allSubsamples, temp)
					
				}
			}
			
			if (progress) {
				utils::setTxtProgressBar(pb, value = currentStep / lastStep)
				currentStep = currentStep + 1
			}
			
		}
	}
	
	if (progress) {
		close(pb)
	}
	
	# Remove default group name from key.
	#TODO: This is a bad hack
	if (resultIsType(res, "WP")) {
		allSubsamples$key = gsub(paste0(defaultGroupName(), ":"), "", allSubsamples$key)
	}
	
	if (summarize) {
	  rval = cleanAndSummarizeMEIResults(allSubsamples, summarize = summarize, aggregateBy = c("parName", "key"))
	} else {
	  rval = allSubsamples
	}
	
	rval
}

