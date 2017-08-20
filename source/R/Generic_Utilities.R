

#' Glossary
#' 
#' @section Types of Result Object:
#' + *Generic results object*: A results object that can be either a WP or a BP results object. If a function takes a generic results object, the name of the results object argument will be `res`.
#' + *WP results object*: A within-participants results object. The return value of [`runParameterEstimation`]. If a function takes a WP results object, the name of the results object argument will be `results`.
#' + BP results object: A between-participants results object. The return value of [`mergeGroupResults.BP`]. If a function takes a BP results object, the name of the results object argument will be `bpRes`.
#' 
#' In this package, most functions can take either a WP or a BP results object. I use stong naming conventions for clarity, so pay attention to the name of the results object argument when calling a function to verify that you are providing the correct type of results object.
#' 
#' @section Latent vs Manifest Parameter Spaces:
#' *Latent parameters* exist in an unexpected/unnatural space: Probability parameters go from `-Inf` to `Inf` and standard deviation parameters can be negative. When the parameters are sampled (with [`runParameterEstimation`]), they exist in the latent space. The priors on the parameters are also with respect to the latest space. Thus, the latent space is the "true" space for the parameters, in the sense of statistical theory and given the specification of the model.
#' 
#' *Manifest parameters* exist in the expected/natural space: Probability parameters are between 0 and 1 and standard deviation parameters are positive, being forced to be greater than or equal to `results$config$minSD`.
#' 
#' @section Condition Effects:
#' In terms of the model specifications, a *condition effect* is a parameter that accounts for differences between the different within-participants conditions in the design. This is explained in more detail in the Condition Effects section in the Appendix of Hardman, Vergauwe, and Ricker (2017). In the equation on page 22, `P_j` is the condition effect parameter that is added to the participant parameter, `P_i`.
#' 
#' @md
#' @name Glossary
NULL


###############################################################################

#' Remove Burn-In Iterations
#' 
#' Due to the large number of parameters in the model and the fact that many of the paramters have Metropolis-Hastings updating steps, convergence can be slow. This function removes burn-in iterations from the beginning of the chains so that only converged posterior distributions will be analyzed.
#'  
#' Note that the Metropolis-Hastings acceptance rate given by [`examineMHAcceptance`] is unaffected by removing burn-in iterations.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param burnIn Integer scalar or vector. Number of burn-in iterations to remove or a vector with length greater than 1 giving the indices of the iterations to remove.
#' 
#' @return A new results object with burnIn iterations removed.
#' 
#' @family generic functions
#' @md
#' @export
removeBurnIn = function(res, burnIn) {
	
	rval = NULL
	if (resultIsType(res, "WP")) {
		rval = removeBurnIn.WP(res, burnIn)
	} else if (resultIsType(res, "BP")) {
		rval = removeBurnIn.BP(res, burnIn)
	}
	rval
	
}

removeBurnIn.BP = function(bpRes, burnIn) {
	
	for (g in names(bpRes$groups)) {
		bpRes$groups[[g]] = removeBurnIn(bpRes$groups[[g]], burnIn)
	}
	
	bpRes$config$iterations = bpRes$groups[[1]]$config$iterations
	
	bpRes
}

removeBurnIn.WP = function(results, burnIn) {
	
	if (length(burnIn) == 0) {
		warning("No burn-in iterations removed because length(burnIn) == 0.")
		return(results)
	}
	
	if (length(burnIn) == 1) {
		burnIn = 1:burnIn
	}
	
	if (burnIn[1] < 1) {
		stop("The lowest burn-in iteration is too low (less than 1).")
	}
	
	if (burnIn[length(burnIn)] > results$config$iterations) {
		stop("The highest burn-in iteration is too high.")
	}
	
	
	keep = 1:results$config$iterations
	keep = keep[-burnIn]
	
	
	for (n in names(results$posteriors)) {
		results$posteriors[[n]] = results$posteriors[[n]][keep]
	}
	
	results$config$iterations = length(keep)
	
	results
}

###############################################################################

# Check Type of Results Object
#
# @param res A generic results object (see [`Glossary`]).
# @param type One of `"WP"` (within-participants) or `"BP"` (between-participants).
# @family generic functions
# @export
resultIsType = function(res, type) {
	conversion = list(WP = "CCM_WP", BP = "CCM_BP")
	
	if (!any(conversion %in% class(res))) {
		stop("Results object is not of any valid type.")
	}
	
	conversion[[ type ]] %in% class(res)
}


###############################################################################

#' Get Parameter Transformation Function
#' 
#' Returns a function that will take parameter values in the latent space and convert them to the manifest space, or vice versa if `inverse` is `TRUE`. For probability parameters, the transformation is the inverse logit transformation. For standard deviation parameters, the transformation forces the parameter to be greater than some value, given by `results$config$minSD`.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param param Name of a parameter, e.g. `"pMem"`.
#' @param inverse If `TRUE`, the inverse transformation is returned, if possible. Some transformations do not have an inverse.
#' 
#' @return A function of one vector, matrix, or array argument that transforms the argument while maintaining the dimensionality of the argument.
#' 
#' @md
#' @family generic functions
#' @export
getParameterTransformation = function(res, param, inverse=FALSE) {

	if (param %in% getProbParams(NULL)) {
		
		if (inverse) {
			transformation = logit
		} else {
			transformation = logitInverse
		}
		
	} else if (param %in% getSdParams(NULL)) {
		
		if (inverse) {
			transformation = function(x) { x } #Give warning if x == minSd? Say that you don't know the inverse?
		} else {
			minSd = res$config$minSD
			
			transformation = local( function(x) { pmax(x, minSd); })
		}
		
	} else if (param %in% getCategoryParams(NULL)) {
		
		transformation = function(x) { x }
		
	}	else {
		stop( paste("No transformation found for parameter: ", param, ".", sep="") )
	}
	
	transformation
}

###############################################################################

#' Get Vectors of Parameter Names
#' 
#' `getProbParams` gets the names of "probability parameters" (e.g. `pMem`). `getSdParams` gets the names of "SD parameters" (e.g. `contSD`).
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param modelVariant Rather than providing `res`, you may provide this.
#' @param filter Filter out parameters not in the current model variant. If `FALSE`, the `res` and `modelVariant` arguments can be NULL.
#' 
#' @return A character vector of the names of parameters used by the model variant.
#' 
#' @md
#' @export
getProbParams = function(res, modelVariant = res$config$modelVariant, filter = FALSE) {
	pp = c("pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess")
	if (filter) {
		if (modelVariant == "betweenItem") {
			pp = c("pMem", "pContBetween", "pCatGuess")
		} else if (modelVariant == "withinItem") {
			pp = c("pMem", "pContWithin", "pCatGuess")
		} else if (modelVariant == "ZL") {
			pp = c("pMem")
		}
	}
	pp
}

#' @rdname getProbParams
#' @export
getSdParams = function(res, modelVariant = res$config$modelVariant, filter = FALSE) {
	pp = c("contSD", "catSelectivity", "catSD")
	if (filter && modelVariant == "ZL") {
		pp = c("contSD")
	}
	pp
}

#' @rdname getProbParams
#' @export
getAllParams = function(res, modelVariant = res$config$modelVariant, filter = FALSE) {
	c(getProbParams(modelVariant = modelVariant, filter = filter),
		getSdParams(modelVariant = modelVariant, filter = filter))
}


getCategoryParams = function(res, modelVariant = res$config$modelVariant, filter = FALSE) {
	param = c("catMu", "catActive")
	if (filter && modelVariant == "ZL") {
		param = character(0)
	}
	param
}



###############################################################################

#' Names of Factors with which a Parameter Varies
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param param A parameter name.
#' 
#' @return A character vector of factor names.
#' 
#' @md
#' @export
getFactorsForConditionEffect = function(res, param) {
	
	rval = NULL
	if (resultIsType(res, "WP")) {
		rval = getFactorsForConditionEffect.WP(res$config, param)
	} else if (resultIsType(res, "BP")) {
		rval = getFactorsForConditionEffect.BP(res, param)
	}

	if (is.null(rval)) {
		stop("Invalid arguments.")
	}
	
	rval
}

getFactorsForConditionEffect.WP = function(config, param) {
	
	allFactorNames = getAllFactorNames(config$factors)
	
	thisFactors = config$conditionEffects[[param]]
	if (length(thisFactors) == 1) {
		if (thisFactors == "all") {
			thisFactors = allFactorNames
		} else if (thisFactors == "none") {
			thisFactors = character(0)
		}
	}
	
	thisFactors = thisFactors[ thisFactors %in% allFactorNames ]
	
	thisFactors
}

getFactorsForConditionEffect.BP = function(bpRes, param) {
	ce = bpRes$config$conditionEffects[[ param ]]
	bpFactors = getFactorTypeToName(bpRes$config$factors)$bp
	c(ce, bpFactors)
}

###############################################################################

#' Update Factors while accounting for Used Condition Effects
#' 
#' Makes a copy of `res$config$factors` that has been updated to account for the fact that
#' if condition effects are disabled for a particular factor, `factors` should not show 
#' that factor as varying. Thus, if condition effects are disabled for a given factor,
#' the levels of that factor are all set to the level of that factor corresponding to 
#' the cornerstone condition.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param param A parameter name.
#' @param removeConstant Whether constant factors should be removed.
#' 
#' @export
updateFactorsForConditionEffects = function(res, param, removeConstant = FALSE) {
	
	if (resultIsType(res, "WP")) {
		fun = updateFactorsForConditionEffects.WP
	} else if (resultIsType(res, "BP")) {
		fun = updateFactorsForConditionEffects.BP
	}
	
	factors = fun(res, param)
	
	if (removeConstant) {
		factors = removeConstantFactors(factors, warnOnRemoval = FALSE)
	}
	
	factors
}

updateFactorsForConditionEffects.BP = function(bpRes, param) {
	
	factors = normalizeFactors(bpRes$config$factors)
	
	wpFactNames = getFactorTypeToName(factors)$wp
	
	for (grp in names(bpRes$groups)) {
		
		results = bpRes$groups[[ grp ]]
		cef = getFactorsForConditionEffect.WP(results$config, param)
		
		removedFactors = wpFactNames[ !(wpFactNames %in% cef) ]
		
		wpFact = normalizeFactors(results$config$factors)

		for (rf in removedFactors) {
			
			flev = sort(unique(wpFact[ , rf ]))
			
			#flev = wpFact[ wpFact$cond == results$config$cornerstoneConditionName, rf ]
			
			factors[ factors$group == grp, rf ] = paste(flev, collapse="/")
		}

	}
	
	factors
}

updateFactorsForConditionEffects.WP = function(results, param) {
	
	factors = normalizeFactors(results$config$factors, removeConstant = FALSE)
	
	factorNames = getAllFactorNames(factors, removeConstant = FALSE)
	
	condtionEffectFactors = getFactorsForConditionEffect.WP(results$config, param)
	
	for (fn in factorNames) {
		
		if (!(fn %in% condtionEffectFactors)) {
			
			flev = sort(unique(factors[ , fn ]))
			
			factors[ , fn ] = paste(flev, collapse="/")
		}

	}
	
	factors
}

###############################################################################

# @export
#getvaryingFactorNames = function(res, param) {
#	factors = updateFactorsForConditionEffects(res, param, removeConstant = TRUE)
#	getAllFactorNames(factors, removeConstant = FALSE)
#}


###############################################################################

#' @export
getAllFactorNames = function(factors, removeConstant = FALSE) {
	if (removeConstant) {
		factors = removeConstantFactors(factors, warnOnRemoval = FALSE)
	}
	
	bpCols = c("key", "group", "cond")
	allNames = names(factors)
	ns = allNames[ !(allNames %in% bpCols) ]
	ns
}


#' @export
getFactorNameToType = function(factors) {
	factors = normalizeFactors(factors, removeConstant = TRUE)
	
	ns = names(factors)
	
	ns = ns[ !(ns %in% c("cond", "group", "key")) ]
	
	nameToType = list()
	
	for (i in 1:length(ns)) {
		
		form = formula( paste0(ns[i], " ~ group") )
		nunique = function(x) { length(unique(x)) }
		agg = aggregate(form, factors, nunique)
		
		isWP = any(agg[ , ns[i] ] > 1)
		
		if (isWP) {
			nameToType[[ ns[i] ]] = "wp"
		} else {
			nameToType[[ ns[i] ]] = "bp"
		}
		
	}
	
	nameToType
}

#' Map from Factor Type to Factor Name
#' 
#' @param factors A WP/BP factors object.
#' 
#' @return A list with three elements. `all`: all factor names. `bp`: between-participant factor names. `wp`: within-participant factor names.
#' 
#' @family generic functions
#' @export
getFactorTypeToName = function(factors) {
	nameToType = getFactorNameToType(factors)
	
	typeToName = list(wp = character(0), bp = character(0))
	for (n in names(nameToType)) {
		typeToName[[ nameToType[[n]] ]] = c(typeToName[[ nameToType[[n]] ]], n)
		typeToName[[ "all" ]] = c(typeToName[[ "all" ]], n)
	}
	
	typeToName
}



#' @family generic functions
#' @export
removeConstantFactors = function(factors, warnOnRemoval = TRUE) {
	
	allFN = getAllFactorNames(factors)
	removed = NULL
	for (n in allFN) {
		if (length(unique(factors[, n])) == 1) {
			removed = c(removed, n)

			factors[,n] = NULL
		}
	}
	if (warnOnRemoval && !is.null(removed)) {
		warning( paste0("The following factors are constant and have been removed: ", paste(removed, collapse = ", "), ".") )
	}
	
	factors
}

#' @family generic functions
#' @export
normalizeFactors = function(factors, removeConstant = FALSE, warnOnRemoval = TRUE) {
	
	colNames = c("key", "group", "cond")
	ns = names(factors)
	ns = ns[ !(ns %in% colNames) ]

	if (!("group" %in% names(factors))) {
		factors$group = defaultGroupName()
	}
	if (!("key" %in% names(factors))) {
		factors$key = paste0(factors$group, ":", factors$cond)
	}
	
	# Remove any R factors
	for (n in names(factors)) {
		factors[ , n ] = as.character(factors[ , n ])
	}
	
	factors = factors[ , c(colNames, ns) ]
	
	if (removeConstant) {
		factors = removeConstantFactors(factors, warnOnRemoval)
	}
	
	factors
}


# DEPRECIATED
getFactorTypes = function(factors) {
	
	factors = normalizeFactors(factors)
	
	ns = names(factors)
	
	ns = ns[ !(ns %in% c("cond", "group", "key")) ]
	
	nameToType = list()
	
	for (i in 1:length(ns)) {
		
		form = formula( paste0(ns[i], " ~ group") )
		nunique = function(x) { length(unique(x)) }
		agg = aggregate(form, factors, nunique)
		
		isWP = any(agg[ , ns[i] ] > 1)
		
		if (isWP) {
			nameToType[[ ns[i] ]] = "wp"
		} else {
			nameToType[[ ns[i] ]] = "bp"
		}
		
	}
	
	typeToName = list(wp = character(0), bp = character(0))
	for (n in names(nameToType)) {
		typeToName[[ nameToType[[n]] ]] = c(typeToName[[ nameToType[[n]] ]], n)
		typeToName[[ "all" ]] = c(typeToName[[ "all" ]], n)
	}
	
	list(nameToType = nameToType, typeToName = typeToName)
}

###############################################################################
