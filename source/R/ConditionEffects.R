


#' Matrices of Prior and Posterior Condition Effects
#' 
#' Get matrices of prior and posterior condition effects for a given parameter (see [`Glossary`] for more information on what a condition effect is).
#' These can be used to perform tests with, e.g., the CMBBHT package.
#' 
#' Optionally, the grand mean parameter, mu, can be added to the condition effect by setting `addMu` to `TRUE` (default `FALSE`).
#' In addition, the resulting parameter values can be converted to the manifest parameter space by setting `manifest` to `TRUE` (default `FALSE`).
#' Converting a condition effect parameter to the manifest space without first adding in the grand mean would be nonsensical, so setting `addMu` to `FALSE` and `manifest` to `TRUE` at the same time is disallowed.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param parName The name of a single parameter (e.g. `"pMem"`).
#' @param priorSamples The number of samples to take from the priors on the condition effects. Defaults to the number of posterior iterations, which is usually appropriate.
#' @param addMu If `TRUE`, the grand mean of the parameter is added to the condition effects. If `FALSE`, nothing is added.
#' @param manifest If `TRUE`, the resulting parameter will be in the manifest space. If `FALSE`, the resulting parameter will be in the latent space. (See [`Glossary`].)
#' @param prior If `TRUE`, the prior condition effects are returned. Otherwise, no priors are returned.
#' @param posterior If `TRUE`, the posterior condition effects are returned. Otherwise, no posteriors are returned.
#' 
#' @return A list with two matrices, `prior` and `post`. Each column is one condition effect and each row is one sample from the prior or posterior.
#' 
#'
#' @family generic functions
#' 
#' @export
getConditionEffects = function(res, parName, priorSamples = res$runConfig$iterations, addMu = FALSE, manifest = FALSE, prior = TRUE, posterior = TRUE) 
{
	
	if (resultIsType(res, "WP")) {
		fun = getConditionEffects.WP
	} else if (resultIsType(res, "BP")) {
		fun = getConditionEffects.BP
	}
	
	fun(res, parName, priorSamples = priorSamples, addMu = addMu, manifest = manifest, prior = prior, posterior = posterior)
	
}


getDefaultParametersWithConditionEffects = function(modelVariant) {
  pce = c("pMem", "contSD") #ZL and all other models
  
  if (modelVariant == "betweenAndWithin") {
    pce = c(pce, "pBetween", "pContBetween", "pContWithin")
    
  } else if (modelVariant == "betweenItem") {
    pce = c(pce, "pContBetween")
    
  } else if (modelVariant == "withinItem") {
    pce = c(pce, "pContWithin")
  }
  
  pce
}


getPriorConditionEffects = function(results, parName, priorSamples, addMu = FALSE, manifest = FALSE) {
	
	if (!addMu && manifest) {
		stop("A manifest parameter is the transformation of the sum of the latent condition effect and the grand mean. Thus addMu == FALSE && manifest == TRUE is disallowed.")
	}
	
	factors = results$config$factors
	
	mu = 0
	if (addMu) {
		mu = stats::rnorm(priorSamples, 
							 results$priors[[ paste0(parName, "_part.mu.mu") ]], 
							 sqrt(results$priors[[ paste0(parName, "_part.mu.var") ]])
		)
	}
	
	priorDs = matrix(NA, nrow=priorSamples, ncol=nrow(factors))
	
	#first pass: Sample from free parameters and the cornerstone condition
	for (i in 1:nrow(factors)) {
		
		target = paste(parName, "_cond[", factors$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, parName, factors$cond[i])
		
		#If this is its own source (i.e. it is an unconstrained parameter), sample from it
		if (target == rootSource) {
			prior = getConditionParameterPrior(results, parName, factors$cond[i])
			
			# If scale == 0, this is just rep(prior$location, priorSamples)
			x = stats::rcauchy(priorSamples, prior$location, prior$scale)
			
			priorDs[,i] = x + mu
		}
		
	}
	
	# Second pass: Copy equality constrained parameters
	# you can't use different samples for equality constrained parameters
	# because then they would only be equal in distribution, not value.
	for (i in 1:nrow(factors)) {
		
		target = paste(parName, "_cond[", factors$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, parName, factors$cond[i])
		
		if (target != rootSource) {
			
			parts = getConditionParameterParts(rootSource)
			
			sourceCol = which(factors$cond == parts$cond)
			priorDs[,i] = priorDs[ , sourceCol ]
			
		}
	}
	
	if (manifest) {
		trans = getParameterTransformation(results, parName)
		priorDs = trans(priorDs)
	}
	
	colnames(priorDs) = factors$cond
	
	priorDs
	
}

getPosteriorConditionEffects = function(results, parName, addMu = FALSE, manifest = FALSE) {
	
	if (!addMu && manifest) {
		stop("A manifest parameter is the transformation of the sum of the latent condition effect and the grand mean. Thus addMu == FALSE && manifest == TRUE is disallowed.")
	}
	
	factors = results$config$factors
	
	mu = 0
	if (addMu) {
		mu = results$posteriors[[ paste0(parName, "_part.mu") ]]
	}
	
	postDs = matrix(NA, nrow=results$runConfig$iterations, ncol=nrow(factors))
	for (i in 1:nrow(factors)) {
		x = results$posteriors[[ paste0(parName, "_cond[", factors$cond[i], "]") ]]
		
		postDs[,i] = x + mu
	}
	
	if (manifest) {
		trans = getParameterTransformation(results, parName)
		postDs = trans(postDs)
	}
	
	colnames(postDs) = factors$cond
	
	postDs
}

getConditionEffects.WP = function(results, parName, priorSamples, addMu, manifest, prior = TRUE, posterior = TRUE) {
	
	colKeys = data.frame(cond = results$config$factors$cond)
	colKeys = normalizeFactors(colKeys)
	
	rval = list(colKeys = colKeys)
	
	if (prior) {
		rval$prior = getPriorConditionEffects(results, parName, priorSamples, addMu, manifest)
		colnames(rval$prior) = colKeys$key
	}
	
	if (posterior) {
		rval$post = getPosteriorConditionEffects(results, parName, addMu, manifest)
		colnames(rval$post) = colKeys$key
	}
	
	rval
}


getConditionEffects.BP = function(bpRes, parName, priorSamples = bpRes$runConfig$iterations, addMu = FALSE, manifest = FALSE, prior = TRUE, posterior = TRUE) {
	
	# TODO: Special case for catActive? (since it is BP, catActive can differ)
	
	cems = list()
	for (grp in names(bpRes$groups)) {
		
		ceff = getConditionEffects.WP(bpRes$groups[[grp]], parName=parName, priorSamples = priorSamples, addMu = addMu, manifest = manifest, prior = prior, posterior = posterior)
		
		tempCK = data.frame(group = grp, 
												cond = ceff$colKeys$cond, stringsAsFactors = FALSE)
		tempCK = normalizeFactors(tempCK) # to create key and order columns
		#tempCK$key = paste0(tempCK$group, ":", tempCK$cond)
		
		cems$colKeys = rbind(cems$colKeys, tempCK)
		
		ppNames = c("prior", "post")[ c(prior, posterior) ]
		for (pp in  ppNames) {
			colnames(ceff[[pp]]) = tempCK$key
			cems[[pp]] = cbind(cems[[pp]], ceff[[pp]])
		}
	}
	
	# Drop keys not in factors
	keep = cems$colKeys$key %in% bpRes$config$factors$key
	cems$prior[ , keep ]
	cems$post[ , keep ]
	cems$colKeys[ keep, ]
	
	cems
}


# This function collapses condition effects, keeping only unique
# combinations of used factors levels.
collapseConditionEffects = function(condEff, factors, usedFactors, uniqueFL = NULL, aggFun = mean) {

	if (is.null(uniqueFL)) {
		uniqueFL = unique(subset(factors, select = usedFactors))
	}
	
	if (nrow(uniqueFL) > 0) {
		keys = getMatchingKeysForUniqueFL(factors, uniqueFL)
	} else {
		keys = list(factors$key)
		uniqueFL = list()
		for (fn in getAllFactorNames(factors)) {

			flev = sort(unique(factors[, fn]))
			uniqueFL[[ fn ]] = paste(flev, collapse = "/")
		}
		uniqueFL = as.data.frame(uniqueFL)
	}
	
	ppNames = c("prior", "post")
	ppNames = ppNames[ ppNames %in% names(condEff) ]

	updatedCE = list()
	if ("prior" %in% ppNames) {
		updatedCE$prior = matrix(nrow = nrow(condEff$prior), ncol = length(keys))
	}
	if ("post" %in% ppNames) {
		updatedCE$post = matrix(nrow = nrow(condEff$post), ncol = length(keys))
	}

	for (i in 1:length(keys)) {
		for (ppi in ppNames) {
			
			pp = condEff[[ ppi ]]
			pp = subset(pp, select = keys[[ i ]])
			updatedCE[[ ppi ]][ , i ] = apply(pp, 1, aggFun)
			
		}
	}
	
	list(condEff = updatedCE, uniqueFL = uniqueFL, keys = keys)
	
}



#' Prior and/or Posterior Distributions of Main Effect and Interaction Parameters
#' 
#' This function provides matrices of prior and posterior main effect and interaction (MEI) parameters for given `testedFactors`. The results can be used with functions from the CMBBHT package (e.g. `summarizeEffectParameters` or `plotEffectParameterSummary`) or analyzed in other ways.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param parName The name of a single parameter (e.g. `"pMem"`).
#' @param testedFactors See \code{\link[CMBBHT]{testHypothesis}}. Passed directly to \code{\link[CMBBHT]{getEffectParameters}}. A character vector giving the names of factors for which a hypothesis test could be performed. If there is only one factor, the main effect will be used. If there is more than one factor, the interaction of the factors will be the effect that is used. 
#' @param dmFactors See \code{\link[CMBBHT]{testHypothesis}}. Passed directly to \code{\link[CMBBHT]{getEffectParameters}}.
#' @param contrastType See \code{\link[CMBBHT]{testHypothesis}}. Passed directly to \code{\link[CMBBHT]{getEffectParameters}}.
#' @param addMu Passed to same argument of [`getConditionEffects`].
#' @param manifest Passed to same argument of [`getConditionEffects`].
#' @param prior Passed to same argument of [`getConditionEffects`].
#' @param posterior Passed to same argument of [`getConditionEffects`].
#' 
#' @return A list of prior and posterior matrices (depending on the `prior` and `posterior` arguments). Each of the matrices have iterations in rows and effect parameters in columns. The columns are named with the following scheme: "F1.L1:F2.L2" where "Fn" is the name of a factor and "Ln" is the level of that factor and ":" indicates the combinations of factor levels.
#' 
#' @family generic functions
#'
#' @export
getMEIParameters = function(res, parName, testedFactors, dmFactors = testedFactors, contrastType = NULL, addMu = FALSE, manifest = FALSE, prior = TRUE, posterior = TRUE) {
	
	cef = getConditionEffects(res, parName, addMu = addMu, manifest = manifest, prior = prior, posterior = posterior)
	
	gmeihtf = res$config$factors
	gmeihtf[ , c("key", "group", "cond") ] = NULL #NO EXTRA COLUMNS
	
	ppNames = c("prior", "post")[ c(prior, posterior) ]
	efp = list()
	for (pp in  ppNames) {
		efp[[ pp ]] = CMBBHT::getEffectParameters(cellMeans = cef[[ pp ]], factors = gmeihtf, testedFactors = testedFactors, dmFactors = dmFactors, contrastType = contrastType)
	}
	
	efp
}


#' Credible Intervals for Manifest Condition Effect Priors
#' 
#' Which is to say, what prior beliefs do you have about the distribution of condition effects?
#' This function helps deal with the fact that the interpretation of priors on condition effects 
#' 1) depends on the value of the participant parameters and
#' 2) need to be translated from the latent space to the manifest space.
#' 
#' Given a vector of participant parameter values and a prior scale for the condition effect parameter,
#' this:
#' 1) Samples from the condition effect prior (0 centered).
#' 2) Calculates the manifest participant parameter values for each sample from the condition effect prior and
#' 3) Calculates the prior credible interval for the manifest parameter values.
#' This information can help you choose a prior scale value that captures your beliefs about how much the conditions differ from one another (or, more exactly, how much the non-cornerstone conditions differ from the cornerstone condition).
#' 
#' The left plot shows the prior median and the lower and upper credible interval bounds, giving an overall picture of the prior credible interval across the participant parameter values.
#' 
#' The right plot gives the total width of the credible interval and the directional widths from the median to the upper and lower bounds.
#'  
#' @param parName The name of a single parameter (e.g. `"pMem"`).
#' @param p_i A vector of manifest participant parameter values. Typically a series of numbers (used for x-axis in plotting). For probability parameters, use something like `seq(0, 1, 0.025)` for the whole range of probabilities. For SD parameters, use something like `seq(0, 40, 1)`.
#' @param ce_scale The scale of the Cauchy prior on the condition effect parameter.
#' @param cip Proportion of the prior inside of the credible interval.
#' @param n Number of samples to take from the prior on the credible interval. Use more for a more accurate approximation.
#' @param minSD Only used for standard deviation parameters. The minimum standard deviation (i.e. `config$minSD`).
#' @param plot Whether to make plots.
#' @param doMFRow If `TRUE`, uses `par(mfrow=c(1,2))` to make the two plot panels.
#' 
#' @return Invisibly, a `data.frame` containing columns:
#' * `p_i`: The participant parameter value (copied from the `p_i` argument).
#' * `lower`, `upper`: Lower and upper bounds of the credible interval.
#' * `median`: The median of the prior.
#' * `lowerW`, `upperW`: The distance from the median to the lower and upper credible intervals.
#' 
#' @examples \dontrun{
#' conditionEffectPriorCredibleInterval("pMem", seq(0, 1, 0.025), 0.3)
#' conditionEffectPriorCredibleInterval("contSD", seq(0, 40, 1), 3)
#' }
#' 
#' @export
conditionEffectPriorCredibleInterval = function(parName, p_i, ce_scale, cip = 0.95, n = 1e6, minSD = 1, plot = TRUE, doMFRow = TRUE) {
  
  qp = c((1 - cip) / 2, 0.5, (1 + cip) / 2)
  
  # Use same samples for each value of x
  ce = stats::rcauchy(n, 0, ce_scale)
  
  res = list(config = list(minSD = minSD)) # Fake res for getting transformations
  trans = getParameterTransformation(res, parName)
  inverse = getParameterTransformation(res, parName, inverse = TRUE)
  
  df = data.frame(p_i = p_i)
  
  for (i in 1:length(p_i)) {
    
    manifest = trans(inverse(p_i[i]) + ce)
    
    qs = stats::quantile(manifest, qp)
    df$lower[i] = qs[1]
    df$median[i] = qs[2]
    df$upper[i] = qs[3]
  }
  
  df$lowerW = df$median - df$lower
  df$upperW = df$upper - df$median
  df$totalW = df$upper - df$lower
  
  if (plot) {
    
    lowCol = "red"
    upCol = "blue"
    
    if (doMFRow) {
      graphics::par(mfrow=c(1, 2))
    }
    
    ylimL = range(df[ , c("lower", "median", "upper")])
    if (parName %in% getParamNames(types="prob")) {
      ylimL = c(0, 1)
    }
    
    plot(df$p_i, df$median, ylim=ylimL, type='l', xlab=parName, ylab="Manifest Parameter Value")
    graphics::lines(df$p_i, df$lower, col=lowCol)
    graphics::lines(df$p_i, df$upper, col=upCol)
    
    ylimR = range(df[ , c("lowerW", "totalW", "upperW")])
    
    plot(df$p_i, df$totalW, type = 'l', ylim=ylimR, xlab=parName, ylab="Credible Interval Width")
    graphics::lines(df$p_i, df$upperW, col=upCol)
    graphics::lines(df$p_i, df$lowerW, col=lowCol)
    graphics::legend("bottom", legend = c("upper", "total", "lower"), col=c(upCol, "black", lowCol), lty=1)
    
    if (doMFRow) {
      graphics::par(mfrow=c(1, 1))
    }
  }
  
  invisible(df)
}


#' Plot Manifest Condition Effect Prior Histogram
#' 
#' @param parName The name of a single parameter (e.g. `"pMem"`).
#' @param p_i A manifest participant parameter value. E.g., for probability parameters, give a probability.
#' @param ce_scale The scale of the Cauchy prior on the condition effect parameter.
#' @param sdCutoff For standard deviation parameters, extremely large sample values are common, so for plotting purposes the plot has to be cut off somewhere. `sdCutoff` sets the cutoff.
#' @param n Number of samples to take from the prior on the credible interval. Use more for a more accurate approximation.
#' @param minSD Only for standard deviation parameters. The minimum standard deviation (i.e. `config$minSD`).
#' 
#' @return Invisibly, the vector of manifest parameter samples from the condition effect prior that were used for plotting (so some samples are cut off for SD parameters; see the `sdCutoff` argument).
#' 
#' NOT EXPORTED
conditionEffectPriorHist = function(parName, p_i, ce_scale, sdCutoff=50, n=1e6, minSD=1) {
  
  if (length(p_i > 1)) {
    logWarning("Only scalar values of p_i are supported. Using the first value.")
    p_i = p_i[1]
  }
  
  ce = stats::rcauchy(n, 0, ce_scale)
  
  res = list(config = list(minSD = minSD))
  trans = getParameterTransformation(res, parName)
  inverse = getParameterTransformation(res, parName, inverse=TRUE)
  
  manifest = trans(inverse(p_i) + ce)
  if (parName %in% getParamNames(types="sd")) {
    manifest = manifest[manifest < sdCutoff]
  }
  
  main = paste0(parName, " = ", p_i, ", scale = ", ce_scale)
  graphics::hist(manifest, 
                 xlab=paste0("Manifest ", parName), 
                 ylab="Prior Density", 
                 main=main, prob=TRUE)
  
  invisible(manifest)
}

