


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
#' @param param The name of the parameter to get condition effects for.
#' @param priorSamples The number of samples to take from the priors on the condition effects. Defaults to the number of posterior iterations, which is usually appropriate.
#' @param addMu If `TRUE`, the grand mean of the parameter is added to the condition effects. If `FALSE`, nothing is added.
#' @param manifest If `TRUE`, the resulting parameter will be in the manifest space. If `FALSE`, the resulting parameter will be in the latent space. (See [`Glossary`].)
#' 
#' @return A list with two matrices, `prior` and `post`. Each column is one condition effect and each row is one sample from the prior or posterior.
#' 
#' @md
#' @family generic functions
#' @export
getConditionEffects = function(res, param, priorSamples = res$config$iterations, addMu = FALSE, manifest = FALSE) 
{
	
	if (resultIsType(res, "WP")) {
		fun = getConditionEffects.WP
	} else if (resultIsType(res, "BP")) {
		fun = getConditionEffects.BP
	}
	
	fun(res, param, priorSamples = priorSamples, addMu = addMu, manifest = manifest)
	
}


getPriorConditionEffects = function(results, param, priorSamples, addMu = FALSE, manifest = FALSE) {
	
	if (!addMu && manifest) {
		stop("A manifest parameter is the transformation of the sum of the latent condition effect and the grand mean. Thus addMu == FALSE && manifest == TRUE is disallowed.")
	}
	
	
	factors = results$config$factors
	
	mu = 0
	if (addMu) {
		mu = rnorm(priorSamples, 
							 results$prior[[ paste0(param, ".mu.mu") ]], 
							 sqrt(results$prior[[ paste0(param, ".mu.var") ]])
		)
	}
	
	priorDs = matrix(NA, nrow=priorSamples, ncol=nrow(factors))
	
	#first pass: Sample from free parameters and the cornerstone condition
	for (i in 1:nrow(factors)) {
		
		target = paste(param, "_cond[", factors$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, param, factors$cond[i])
		
		#If this is its own source (i.e. it is an unconstrained parameter), sample from it
		if (target == rootSource) {
			prior = getConditionParameterPrior(results, param, factors$cond[i])
			
			# If scale == 0, this is just rep(prior$location, priorSamples)
			x = stats::rcauchy(priorSamples, prior$location, prior$scale)
			
			priorDs[,i] = x + mu
		}
		
	}
	
	# Second pass: Copy equality constrained parameters
	# you can't use different samples for equality constrained parameters
	# because then they would only be equal in distribution, not value.
	for (i in 1:nrow(factors)) {
		
		target = paste(param, "_cond[", factors$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, param, factors$cond[i])
		
		if (target != rootSource) {
			
			parts = getConditionParameterParts(rootSource)
			
			sourceCol = which(factors$cond == parts$cond)
			priorDs[,i] = priorDs[ , sourceCol ]
			
		}
	}
	
	if (manifest) {
		trans = getParameterTransformation(results, param)
		priorDs = trans(priorDs)
	}
	
	colnames(priorDs) = factors$cond
	
	priorDs
	
}

getPosteriorConditionEffects = function(results, param, addMu = FALSE, manifest = FALSE) {
	
	if (!addMu && manifest) {
		stop("A manifest parameter is the transformation of the sum of the latent condition effect and the grand mean. Thus addMu == FALSE && manifest == TRUE is disallowed.")
	}
	
	
	factors = results$config$factors
	
	mu = 0
	if (addMu) {
		mu = results$posteriors[[ paste0(param, ".mu") ]]
	}
	
	postDs = matrix(NA, nrow=results$config$iterations, ncol=nrow(factors))
	for (i in 1:nrow(factors)) {
		x = results$posteriors[[ paste0(param, "_cond[", factors$cond[i], "]") ]]
		
		postDs[,i] = x + mu
	}
	
	if (manifest) {
		trans = getParameterTransformation(results, param)
		postDs = trans(postDs)
	}
	
	colnames(postDs) = factors$cond
	
	postDs
}



getConditionEffects.WP = function(results, param, priorSamples = results$config$iterations, addMu = FALSE, manifest = FALSE) {
	prior = getPriorConditionEffects(results, param, priorSamples, addMu, manifest)
	post = getPosteriorConditionEffects(results, param, addMu, manifest)
	
	colKeys = data.frame(cond = results$config$factors$cond)
	colKeys = normalizeFactors(colKeys)
	
	colnames(prior) = colKeys$key
	colnames(post) = colKeys$key
	
	list(prior=prior, post=post, colKeys = colKeys)
}


getConditionEffects.BP = function(bpRes, param, priorSamples = bpRes$config$iterations, addMu = FALSE, manifest = FALSE) {
	
	# TODO: Special case for catActive? (since it is BP, catActive can differ)
	
	cems = list(prior = NULL, post = NULL, colKeys = NULL)
	for (grp in names(bpRes$groups)) {
		
		ceff = CatContModel:::getConditionEffects.WP(bpRes$groups[[grp]], param=param, priorSamples = priorSamples, addMu = addMu, manifest = manifest)
		
		tempCK = data.frame(group = grp, 
												cond = ceff$colKeys$cond)
		tempCK = normalizeFactors(tempCK) # to create key and order columns
		#tempCK$key = paste0(tempCK$group, ":", tempCK$cond)
		
		cems$colKeys = rbind(cems$colKeys, tempCK)
		
		for (pp in c("prior", "post")) {
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
#' @export
collapseConditionEffects = function(condEff, factors, usedFactors, uniqueFL = NULL, aggFun = mean) {

	if (is.null(uniqueFL)) {
		uniqueFL = unique(subset(factors, select = usedFactors))
	} else {
		usedFactors = names(uniqueFL)
	}
	
	keys = CatContModel:::getMatchingKeysForUniqueFL(factors, uniqueFL)

	updatedCE = list(prior = matrix(nrow = nrow(condEff$prior), ncol = length(keys)), 
									 post = matrix(nrow = nrow(condEff$post), ncol = length(keys))
									 )
	
	for (i in 1:length(keys)) {
		for (ppi in c("prior", "post")) {
			
			pp = condEff[[ ppi ]]
			pp = subset(pp, select = keys[[ i ]])
			updatedCE[[ ppi ]][ , i ] = apply(pp, 1, aggFun)
			
		}
	}
	
	list(condEff = updatedCE, uniqueFL = uniqueFL, keys = keys)
	
}
