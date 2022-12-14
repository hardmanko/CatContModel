

#' Glossary of Common Terms
#' 
#' This glossary covers some of the most common jargon used in this package.
#' 
#' @section Types of Result Object:
#' + *Generic results object*: A results object that can be either a WP or a BP results object. If a function takes a generic results object, the name of the results object argument will be `res`. Many functions take generic results objects.
#' + *WP results object*: A within-participants results object. The return value of [`runParameterEstimation`]. If a function takes a WP results object, the name of the results object argument will be `results`.
#' + *BP results object*: A between-participants results object. The return value of [`combineGroupResults.BP`]. If a function takes a BP results object, the name of the results object argument will be `bpRes`.
#' + *Parallel results object*: The return value of [`runParameterEstimation`] if the `parallel` argument is provided.
#' 
#' In this package, most functions can take either a WP or a BP results object. Only some functions can take Parallel objects.
#' Pay attention to the name of the results object argument when calling a function to verify that you are providing the correct type of results object.
#' The type of a results object can be checked with [`resultIsType`] or the R function `class`.
#' 
#' @section Latent vs Manifest Parameter Spaces:
#' *Latent parameters* exist in an unexpected/unnatural space: 
#' Probability parameters go from `-Inf` to `Inf` and standard deviation parameters can be negative. 
#' When the parameters are sampled (with [`runParameterEstimation`]), they exist in the latent space. 
#' The priors on the parameters are also with respect to the latest space. 
#' Thus, the latent space is the "true" space for the parameters, in the sense of statistical theory and given the specification of the model.
#' 
#' *Manifest parameters* exist in the expected/natural space: 
#' Probability parameters are between 0 and 1 and standard deviation parameters are positive, 
#' being forced to be greater than or equal to `config$minSD`.
#' 
#' Functions to transform between spaces are returned by [`getParameterTransformation`].
#' 
#' @section Condition Effects:
#' In terms of the model specifications, a *condition effect* is a parameter that accounts for differences between 
#' within-participants conditions in the design. This is explained in more detail in the Condition Effects section 
#' in the Appendix of Hardman, Vergauwe, and Ricker (2017). 
#' In the equation on page 22, `P_j` is the condition effect parameter that is added to the participant parameter, `P_i`.
#' 
#' @name Glossary
NULL

###############################################################################

#' Check Type of Results Object
#'
#' @param res A generic results object (see [`Glossary`]).
#' @param type One of `"WP"` (within-participants), `"BP"` (between-participants), or `"Parallel"` (WP estimated in parallel).
#' @family generic functions
#' @export
resultIsType = function(res, type) {
  conversion = list(WP = "CCM_WP", BP = "CCM_BP", Parallel = "CCM_Parallel")
  
  if (!any(conversion %in% class(res))) {
    stop("Results object is not of any valid type.")
  }
  
  nMatch = 0
  for (i in 1:length(type)) {
    nMatch = nMatch + as.numeric(conversion[[ type[i] ]] %in% class(res))
  }
  
  nMatch >= 1
  
  #conversion[[ type ]] %in% class(res)
}

##################################

getCompletedIterations = function(res, check=TRUE) {
  rval = NA
  if (resultIsType(res, "WP")) {
    rval = res$runConfig$iterations
    
  } else if (resultIsType(res, "BP")) {
    
    groupIter = rep(0, length(res$groups))
    for (i in 1:length(res$groups)) {
      groupIter[i] = res$groups[[i]]$runConfig$iterations
    }
    
    if (check && any(groupIter != groupIter[1])) {
      warning("Different between-participants groups have different numbers of iterations.")
    }
    
    rval = groupIter[1]
    
  } else if (resultIsType(res, "Parallel")) {
    
    chainIter = rep(0, length(res$chains))
    for (i in 1:length(res$chains)) {
      chainIter[i] = res$chains[[i]]$runConfig$iterations
    }
    
    if (check && any(chainIter != chainIter[1])) {
      warning("Different parallel chains have different numbers of iterations.")
    }
    
    rval = chainIter[1]
    #rval = res$chains[[1]]$runConfig$iterations
  }
  rval
}

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
#'
#' @export
removeBurnIn = function(res, burnIn) {
	
	rval = NULL
	if (resultIsType(res, "WP")) {
		rval = removeBurnIn.WP(res, burnIn)
	} else if (resultIsType(res, "BP")) {
		rval = removeBurnIn.BP(res, burnIn)
	} else if (resultIsType(res, "Parallel")) {
	  rval = removeBurnIn.Parallel(res, burnIn)
	} else {
	  stop("Invalid object type.")
	}
	rval
	
}

removeBurnIn.WP = function(results, burnIn) {
	
	if (length(burnIn) == 0) {
		warning("No burn-in iterations removed because length(burnIn) == 0.")
		return(results)
	}
  
	if (length(burnIn) == 1) {
	  if (burnIn <= 0) {
	    warning("No burn-in iterations removed because burnIn <= 0.")
	    return(results)
	  }
		burnIn = 1:burnIn
	}
	
	if (burnIn[1] < 1) {
		stop("The lowest burn-in iteration is too low (less than 1).")
	}
	
	if (burnIn[length(burnIn)] > results$runConfig$iterations) {
		stop("The highest burn-in iteration is too high.")
	}
	
	
	keep = 1:results$runConfig$iterations
	keep = keep[-burnIn]
	
	
	for (n in names(results$posteriors)) {
		results$posteriors[[n]] = results$posteriors[[n]][keep]
	}
	
	results$runConfig$iterations = length(keep)
	
	results
}

removeBurnIn.BP = function(bpRes, burnIn) {
  
  for (g in names(bpRes$groups)) {
    bpRes$groups[[g]] = removeBurnIn(bpRes$groups[[g]], burnIn)
  }
  
  bpRes$runConfig$iterations = bpRes$groups[[1]]$runConfig$iterations
  
  bpRes
}

removeBurnIn.Parallel = function(parRes, burnIn) {
  
  for (i in 1:length(parRes$chains)) {
    parRes$chains[[i]] = removeBurnIn.WP(parRes$chains[[i]], burnIn)
  }
  
  parRes
}


# Or prepareResults
cleanResults = function(res, burnIn) {
  
  if (resultIsType(res, "WP")) {
    rval = removeBurnIn(res, burnIn)
    
  } else if (resultIsType(res, "BP")) {
    rval = removeBurnIn(res, burnIn)
    
  } else if (resultIsType(res, "Parallel")) {
    res = removeBurnIn(res, burnIn)
    rval = combineResults(resList=res$chains)
    
  } else if (resultIsType(res, "CrossValidation")) {
    rval = removeBurnIn(res, burnIn)
    
  }
  
  rval
}

###############################################################################

#' Get Parameter Transformation Function
#' 
#' Returns a function that will take parameter values in the latent space and convert them to the manifest space, or vice versa if `inverse` is `TRUE`. 
#' For probability parameters, the transformation is the inverse logit transformation. 
#' For standard deviation parameters, the transformation forces the parameter to be greater than some value, given by `results$config$minSD`.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param param Name of a parameter, e.g. `"pMem"`.
#' @param inverse If `TRUE`, the inverse transformation is returned, if possible. Some transformations do not have an inverse.
#' @param minSD Minimum value for standard deviation parameters. If `minSD` is provided, `res` may be `NULL`.
#' 
#' @return A function of one vector, matrix, or array argument that transforms the argument while maintaining the dimensionality of the argument.
#' 
#' @family generic functions
#' 
#' @export
getParameterTransformation = function(res, param, inverse=FALSE, minSD=NULL) {

  paramParts = splitParamName(param)
  
  if ("probParam" %in% paramParts$types) {
    if (inverse) {
      transformation = logit
    } else {
      transformation = logitInverse
    }    
  } else if ("sdParam" %in% paramParts$types) {
    
    # Function does not need res, it can just take minSD
    if (is.null(minSD)) {
      minSD = res$config$minSD
    }
    if (is.null(minSD)) {
      stop("No minSD provided")
    }
    
    if (inverse) {
      transformation = local(function(x) {
        if (x < minSD) {
          warning(paste0("In SD parameter inverse transformation function: x < minSD (", x, " < ", minSD, ")"))
        }
        x
      })
    } else {

      transformation = local( function(x) { pmax(x, minSD); } )
    }
  } else if ("catParam" %in% paramParts$types) {
    transformation = function(x) { x }
    
  }	else {
    stop( paste("No transformation found for parameter: ", param, ".", sep="") )
  }
	
	transformation
}



###############################################################################
# Parameter name functions

#' Get Names of Model Parameters
#' 
#' 
#' @param modelVariant If provided, results will only include the names of parameters used by the model variant.
#' @param types A vector of names of parameter types.
#' @param base If `TRUE`, base parameter names are returned.
#' @param cond If `TRUE`, condition effects are added (e.g. "pMem_cond").
#' @param hyper If `TRUE`, hyperpriors are added (e.g. "pMem.mu" and "pMem.var").
#' @param asList If `TRUE`, returned names are formatted in a list instead of vector.
#' 
#' @return A vector of parameter names (unless `asList = TRUE`).
#' 
#' @export
getParamNames = function(modelVariant=NULL, types=c("prob", "sd", "cat"), base=TRUE, cond=FALSE, hyper=FALSE, asList=FALSE) {
  
  probParam = c("pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess")
  sdParam = c("contSD", "catSelectivity", "catSD")
  catParam = c("catMu", "catActive")
  
  if (!is.null(modelVariant)) {
    if (modelVariant == "betweenItem") {
      probParam = c("pMem", "pContBetween", "pCatGuess")
    } else if (modelVariant == "withinItem") {
      probParam = c("pMem", "pContWithin", "pCatGuess")
    } else if (modelVariant == "ZL") {
      probParam = "pMem"
      sdParam = "contSD"
      catParam = c()
    }
  }
  
  if (!("prob" %in% types)) {
    probParam = c()
  }
  if (!("sd" %in% types)) {
    sdParam = c()
  }
  if (!("cat" %in% types)) {
    catParam = c()
  }
  
  if (asList) {
    rval = list(prob=list(base=probParam), sd=list(base=sdParam), cat=list(base=catParam))
    if (cond) {
      rval$prob$cond = paste0(probParam, "_cond")
      rval$sd$cond = paste0(sdParam, "_cond")
    }
    if (hyper) {
      for (suffix in c(".mu", ".var")) {
        rval$prob$hyper = c(rval$prob$hyper, paste0(probParam, suffix))
        rval$sd$hyper = c(rval$sd$hyper, paste0(sdParam, suffix))
      }
    }
  } else {
    # as vector
    rval = c(probParam, sdParam, catParam)
    if (cond) {
      rval = c(rval, paste0(c(probParam, sdParam), "_cond"))
    }
    if (hyper) {
      muPar = paste0(c(probParam, sdParam), ".mu")
      varPar = paste0(c(probParam, sdParam), ".var")
      rval = c(rval, muPar, varPar)
    }
  }
  
  rval
}


# This could also have a list return option.
# rval = list(pMem = list(part=vector("char"), cond=vector("char"), hyper=vector("char")), contSD=list(...))
#
# TODO: This is not used anywhere
getIndexedParamNames = function(data, maxCategories, modelVariant=NULL, types=c("participant", "cond", "hyper")) {
  
  standardPN = getParamNames(modelVariant=modelVariant, types=c("prob", "sd"))
  catPN = getParamNames(modelVariant=modelVariant, types=c("cat"))
  
  indexedNames = vector("character", 0)
  
  pnums = unique(data$pnum)
  conds = unique(data$cond)
  
  for (param in standardPN) {
    if ("participant" %in% types) {
      indexedNames = c(indexedNames, paste0(param, "[", pnums, "]"))
    }
    if ("hyper" %in% types) {
      indexedNames = c(indexedNames, paste0(param, c(".mu", ".var")))
    }
    if ("cond" %in% types) {
      indexedNames = c(indexedNames, paste0(param, "_cond[", conds, "]"))
    }
  }
  
  if (is.null(maxCategories)) {
    catInd = ":" # matlab style
  } else {
    catInd = 1:maxCategories
  }
  
  for (cp in catPN) {
    for (pnum in pnums) {
      indexedNames = c(indexedNames, paste0(cp, "[", pnum, ",", catInd, "]"))
    }
  }
  
  indexedNames
}


#' Split Parameter Name
#' 
#' Utility that splits indexed parameter names (like `pMem[123]`) into parts.
#' 
#' @param pname Indexed parameter name.
#' 
#' @return A list with information about the parameter name.
#' 
#' @export
splitParamName = function(pname) {
  
  rval = list(baseName="invalid", types="invalid")
  
  if (grepl(".mu", pname, fixed=TRUE)) {
    
    # Hyperprior mean
    
    rval$baseName = stringr::str_remove(pname, ".mu")
    rval$types = c("mu", "hyper")
    
  } else if (grepl(".var", pname, fixed=TRUE)) {
    
    # Hyperprior variance
    
    rval$baseName = stringr::str_remove(pname, ".var")
    rval$types = c("var", "hyper")
    
  } else if (grepl("_cond", pname, fixed=TRUE)) {
    
    # Condition effect
    
    parts = stringr::str_split_fixed(pname, "_cond", n=Inf)
    
    rval$index = stringr::str_remove_all(parts[2], stringr::fixed("["))
    rval$index = stringr::str_remove_all(rval$index, stringr::fixed("]"))
    
    rval$baseName = parts[1]
    rval$types = "cond"
    
  } else if (grepl("[", pname, fixed=TRUE)) {
    
    # Participant parameter
    
    parts = stringr::str_split_fixed(pname, stringr::fixed("["), n=Inf)
    
    rval$index = stringr::str_remove_all(parts[2], stringr::fixed("]"))
    
    if (grepl(",", rval$index, fixed=TRUE)) {
      arrayParts = stringr::str_split_fixed(rval$index, stringr::fixed(","), n=Inf)
      rval$index = as.vector(arrayParts)
    }
    
    rval$baseName = parts[1]
    rval$types = "participant"
    
  } else {
    rval$baseName = pname
    rval$types = "generic"
  }
  
  if (!(rval$baseName %in% getParamNames())) {
    rval$types = c(rval$types, "invalid")
  }
  
  if (rval$baseName %in% getParamNames(types="prob")) {
    rval$types = c(rval$types, "probParam")
  } else if (rval$baseName %in% getParamNames(types="sd")) {
    rval$types = c(rval$types, "sdParam")
  } else if (rval$baseName %in% getParamNames(types="cat")) {
    rval$types = c(rval$types, "catParam")
  }
  
  rval
}



###############################################################################

#' Names of Factors with which a Parameter Varies
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param param A parameter name.
#' 
#' @return A character vector of factor names.
#' 
#' @family generic functions
#' 
#' @export
getFactorsForConditionEffect = function(res, param) {
	
	rval = NULL
	if (resultIsType(res, "WP")) {
		rval = getFactorsForConditionEffect.WP(res$config, param)
	} else if (resultIsType(res, "BP")) {
		rval = getFactorsForConditionEffect.BP(res, param)
	} else if (resultIsType(res, "Parallel")) {
	  stop("This function does not support Parallel results.")
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
	} else if (resultIsType(res, "Parallel")) {
	  stop("This function does not support Parallel results.")
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

#' Names of all Factors
#' 
#' @param factors A factors `data.frame`, e.g. `results$config$factors`.
#' @param removeConstant Whether constant (non-varying) factors should be included in the returned factor names.
#' 
#' @return A character vector of factor names.
#' 
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

#' Map from Factor Name to Factor Type
#' 
#' There are two types of factors, those that vary within groups (WP factors)
#' and those that vary between groups (BP factors). This function returns a `list`
#' that maps from the names of factors to the types of those factors.
#' 
#' @param factors A factors `data.frame`, e.g. `results$config$factors`.
#' 
#' @return A list with factor names as keys and either the value of `"wp"` or `"bp"` for each factor name.
#' 
#' @export
getFactorNameToType = function(factors) {
  
	factors = normalizeFactors(factors, removeConstant = TRUE, warnOnRemoval = FALSE)
	
	# ns is names of factors, not the base names like "cond"
	ns = names(factors)
	
	ns = ns[ !(ns %in% c("cond", "group", "key")) ]
	
	nameToType = list()
	
	if (length(ns) > 0) {
  	for (i in 1:length(ns)) {
  		
  		form = stats::formula( paste0(ns[i], " ~ group") )
  		nunique = function(x) { length(unique(x)) }
  		agg = stats::aggregate(form, factors, nunique)
  		
  		isWP = any(agg[ , ns[i] ] > 1)
  		
  		if (isWP) {
  			nameToType[[ ns[i] ]] = "wp"
  		} else {
  			nameToType[[ ns[i] ]] = "bp"
  		}
  		
  	}
	}
	
	nameToType
}

#' Map from Factor Type to Factor Name
#' 
#' @param factors A factors `data.frame`, e.g. `results$config$factors`.
#' 
#' @return A list with three elements. `all`: all factor names. `bp`: between-participant factor names. `wp`: within-participant factor names.
#' 
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


#' Remove Constant (Non-Varying) Factors
#' 
#' @param factors A factors `data.frame`, e.g. `results$config$factors`.
#' @param warnOnRemoval Whether a warning should be emitted when factors are removed.
#' 
#' @return The argument `factors` with constant factors removed.
#' 
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

#' Normalize the Format of Factors
#' 
#' Normalizing a factors `data.frame` involves
#' 1. Creating the `group` and `key` columns if those columns did not exist in the original.
#' 2. Converting R factors (i.e. those created with \code{\link[base]{factor}}) to character.
#' 3. Optionally, removing constant factors.
#' 
#' @param factors A factors `data.frame`, e.g. `results$config$factors`.
#' @param removeConstant Whether constant (non-varying) factors should be included in the returned factor names.
#' @param warnOnRemoval Whether a warning should be emitted when factors are removed.
#' 
#' @return The argument `factors` in a normalized format.
#' 
#' @export
normalizeFactors = function(factors, removeConstant = FALSE, warnOnRemoval = TRUE) {
	
	colNames = c("key", "group", "cond")
	ns = names(factors)
	ns = ns[ !(ns %in% colNames) ]

	# Create group and key columns if they did not exist.
	if (!("group" %in% names(factors))) {
		factors$group = defaultGroupName()
	}
	if (!("key" %in% names(factors))) {
		factors$key = paste0(factors$group, ":", factors$cond)
	}
	
	# Remove any R factors
	for (n in names(factors)) {
		if (is.factor(factors[ , n ])) {
			factors[ , n ] = as.character(factors[ , n ])
		}
	}
	
	# Reorder columns
	factors = factors[ , c(colNames, ns) ]
	
	if (removeConstant) {
		factors = removeConstantFactors(factors, warnOnRemoval)
	}
	
	factors
}



###############################################################################



#' Population/Condition Posterior Means and Credible Intervals
#' 
#' Calculates posterior means and credible intervals for the population means in each combination of relevant factor levels for the
#' given parameters. Irrelevant factors (i.e. ones that do not vary between the of cells of the design or for which no condition effects were estimated) are collapsed across. 
#' 
#' For each condition, condition effects are added to population means (if `addMu == TRUE`), the result is 
#' transformed to the manifest space (if `manifest == TRUE`), and the mean and credible interval for the manifest value is calculated.
#' Note that this is different from adding condition effects to participant-level parameters, tranforming 
#' the result, calculating on each iteration the mean of the transformed participant parameters, and 
#' calculating the posterior mean and credible interval of the iteration means. 
#' Using iteration means rather than population means will generally result in less than the true
#' amount of variability, which is why population means are used. Note, however, that this is a little strange,
#' because in the model, condition effects are not added to population means, but participant means.
#'
#' @param res A generic results object (see [`Glossary`]).
#' @param params A vector of parameter names. If `NULL`, the default, all valid parameters are used. Unlike most functions, `"catActive"` can be used as a parameter.
#' @param cip The credible interval proportion. Defaults to 95% credible intervals.
#' @param addMu See [`getConditionEffects`].
#' @param manifest See [`getConditionEffects`].
#' 
#' @return A `data.frame` containing the results. The factors of the design will each have their own column. For each parameter and meaningful combination of factor levels, the mean and credible interval are given.
#' 
#' @family generic functions
#'
#' @export
posteriorMeansAndCredibleIntervals = function(res, params=NULL, cip=0.95, addMu=TRUE, manifest=TRUE) {
	
  if (resultIsType(res, "Parallel")) {
    stop("This function does not support Parallel results.")
  }
  
	aggFuns = list(mean = mean,
								 lower = function(x) { stats::quantile(x, (1 - cip) / 2) },
								 upper = function(x) { stats::quantile(x, (1 + cip) / 2) })
	
	if (is.null(params)) {
		#params = getAllParams(res, filter=TRUE)
		params = getParamNames(res$config$modelVariant, types=c("prob", "sd"))
		if (res$config$modelVariant != "ZL") {
		  params = c(params, "catActive") # PMCI is able to summarize catActive
		}
	}

	
	factors = normalizeFactors(res$config$factors)
	allFN = getAllFactorNames(factors)
	
	fullFormula = stats::formula( paste( "x ~ ", paste(allFN, collapse = " * ") ) )
	
	allAgg = NULL
	
	# catActive special case
	if ("catActive" %in% params) {
		
		params = params[ params != "catActive" ]
		
		ica = getIterationCatActive(res)
		
		theseAgg = stats::aggregate(fullFormula, ica, function(x) { NA })
		theseAgg$x = NULL
		theseAgg$param = "catActive"
		
		for (fn in names(aggFuns)) {
			agg = stats::aggregate(fullFormula, ica, aggFuns[[ fn ]])
			theseAgg[ , fn ] = agg$x
		}
		
		allAgg = rbind(allAgg, theseAgg)
		
	}
	
	# The rest of the parameters
	for (param in params) {
		
		condEff = getConditionEffects(res, param, prior = FALSE, posterior = TRUE, 
																	addMu = addMu, manifest = manifest)
		
		ceFactors = getFactorsForConditionEffect(res, param)
		cce = collapseConditionEffects(condEff, factors, usedFactors = ceFactors)
		
		# If any factors are used, copy over missing factors
		if (length(ceFactors) > 0) {
		  baseUniqueFL = cce$uniqueFL # Make a copy so that you can modify the original
		  unusedFactors = allFN[ !(allFN %in% ceFactors) ]
		  for (fn in unusedFactors) {
		    
	      matchingFL = getMatchingFactorLevels(factors, baseUniqueFL, fn, collapse="/")
	      for (i in 1:length(matchingFL)) {
	        cce$uniqueFL[i,fn] = matchingFL[i]
	      }
		  }
		}
		# If no factors are used, cce$uniqueFL will already have collapsed factor levels in it
		
		df = reshapeMatrixToDF(cce$condEff$post, cce$uniqueFL)
		
		theseAgg = stats::aggregate(fullFormula, df, function(x) { NA })
		theseAgg$x = NULL
		theseAgg$param = param
		
		for (fn in names(aggFuns)) {
			agg = stats::aggregate(fullFormula, df, aggFuns[[ fn ]])
			theseAgg[ , fn ] = agg$x
		}
		
		allAgg = rbind(allAgg, theseAgg)
		
	}
	
	colNames = names(allAgg)
	allAgg = allAgg[ , c("param", colNames[ colNames != "param" ]) ]
	allAgg
}


###############################################################################

# collapses across columns in "drop" within iterations
getIterationCatActive = function(res, aggFun = mean, drop = "pnum") {
	
	df = getCatActiveDataFrame(res)
	
	ns = names(df)
	usedCols = ns[ !(ns %in% c("x", drop)) ]
	
	# Average across drop within iterations
	df$iteration = 1:res$runConfig$iterations
	form = paste0("x ~ iteration * ", paste(usedCols, collapse=" * "))
	partMean = stats::aggregate(stats::formula(form), df, aggFun)
	partMean$iteration = NULL
	
	partMean
}

###############################################################################

getCatActiveDataFrame = function(res) {
	
	rval = NULL
	
	if (resultIsType(res, "WP")) {
		rval = getCatActiveDataFrame.WP(res)
	} else if (resultIsType(res, "BP")) {
		rval = getCatActiveDataFrame.BP(res)
	} else if (resultIsType(res, "Parallel")) {
	  stop("This function does not support Parallel results.")
	}
	
	rval
}

getCatActiveDataFrame.BP = function(bpRes) {
	
	type2name = getFactorTypeToName(bpRes$config$factors)
	
	post = convertPosteriorsToMatrices(bpRes, param = "catActive")
	
	# Keep iteration and pnum, but sum across catActive parameters
	caSum = apply(post$catActive, c(1, 3), sum)
	caSum = t(caSum)
	
	group_part = dimnames(caSum)[[2]]
	gpParts = strsplit(group_part, split = ":", fixed=TRUE)
	
	design = NULL
	for (i in 1:length(gpParts)) {
		temp = data.frame(group = gpParts[[i]][1], pnum = gpParts[[i]][2], stringsAsFactors = FALSE)
		design = rbind(design, temp)
	}
	
	# Copy in factor levels
	gfact = bpRes$config$factors
	gfact = gfact[ , c("group", type2name$bp), drop=FALSE ]
	gfact = unique(gfact)
	
	bpFact = gfact[ , type2name$bp, drop=FALSE ]
	
	if (ncol(bpFact) > 0) {
		for (fn in type2name$all) {
			mfl = getMatchingFactorLevels(bpRes$config$factors, bpFact, factor = fn, collapse="/")
			design[ , fn ] = substituteValues(design$group, gfact$group, mfl)
		}
	}
	
	# Convert matrix to DF
	reshapeMatrixToDF(caSum, design)
	
}

getCatActiveDataFrame.WP = function(results) {
	
	post = convertPosteriorsToMatrices(results, param = "catActive")
	
	# Keep iteration and pnum, but sum across catActive parameters
	caSum = apply(post$catActive, c(1, 3), sum)
	caSum = t(caSum)
	
	design = data.frame(group = defaultGroupName(), 
											pnum = dimnames(caSum)[[2]], 
											stringsAsFactors = FALSE)
	
	for (fn in getAllFactorNames(results$config$factors)) {
		f = results$config$factors[ , fn ]
		design[ , fn ] = paste( sort(unique(f)), collapse="/")
	}
	
	# Convert matrix to DF
	reshapeMatrixToDF(caSum, design)
	
}

###############################################################################

#' Convert Posterior Distributions to a Single Matrix
#' 
#' This converts raw posteriors into a single matrix. This matrix can then be used with the boa or coda packages for assessing convergence. Some of the convergence diagnostics require you to make separate matrices for separate runs of the Gibbs sampler.
#' 
#' @section BP designs:
#' For BP designs, the parameter names in the resulting matrix are preceded by "grp:", where "grp" is the name of the group that the parameter is from.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param stripConstantParameters Remove all parameters with a constant value. Constant parameters cannot converge.
#' @param stripCatActive Remove all of the cat active parameters. It is difficult to assess convergence for indicator parameters that are either 0 or 1.
#' @param stripCatMu Remove all of the category mean parameters. It is difficult to assess convergence for these parameters because they often have multi-modal posterior distributions.
#' 
#' @return A matrix containing all of the posterior distributions for the selected parameters. Each column is one parameter.
#' 
#' @family generic functions
#' 
#' @export
#' 
convertPosteriorsToMatrix = function(res, stripConstantParameters=TRUE, stripCatActive=TRUE, stripCatMu=TRUE) {
	
	if (resultIsType(res, "WP")) {
		fun = convertPosteriorsToMatrix.WP
	} else if (resultIsType(res, "BP")) {
		fun = convertPosteriorsToMatrix.BP
	} else if (resultIsType(res, "Parallel")) {
	  stop("This function does not support Parallel results.")
	}
	
	fun(res, stripConstantParameters, stripCatActive, stripCatMu)
}

convertPosteriorsToMatrix.WP = function(results, stripConstantParameters=TRUE, stripCatActive=TRUE, stripCatMu=TRUE) {
	
	rawPost = results$posteriors
	
	iterations = 0
	constantPar = NULL
	catActivePar = NULL
	catMuPar = NULL
	participantLLPar = NULL
	
	for (n in names(rawPost)) {
		iterations = max(c(iterations, length(rawPost[[n]])))
		
		if (all(rawPost[[n]] == rawPost[[n]][1])) {
			constantPar = c(constantPar, n)
		}
		
		if (grepl("catActive", n, fixed=TRUE)) {
			catActivePar = c(catActivePar, n)
		}
		
		if (grepl("catMu", n, fixed=TRUE)) {
			catMuPar = c(catMuPar, n)
		}
		
		if (grepl("participantLL", n, fixed=TRUE)) {
			participantLLPar = c(participantLLPar, n)
		}
		
	}
	
	excludedPar = participantLLPar
	if (stripConstantParameters) {
		excludedPar = c(excludedPar, constantPar)
	}
	if (stripCatActive) {
		excludedPar = c(excludedPar, catActivePar)
	}
	if (stripCatMu) {
		excludedPar = c(excludedPar, catMuPar)
	}
	
	usedParam = names(rawPost)
	usedParam = usedParam[ !(usedParam %in% excludedPar) ]
	
	np = length(usedParam)
	
	m = matrix(0, nrow=iterations, ncol=np)
	colnames(m) = usedParam
	
	for (i in 1:np) {
		
		temp = rawPost[[usedParam[i]]]
		
		if (length(temp) == 1) {
			temp = rep(temp, iterations)
		} else if (length(temp) != iterations) {
			warning(paste("Data vector for parameter ", usedParam[i], " of incorrect length.", sep=""))
		}
		
		m[,i] = temp
		
	}
	
	m
}

convertPosteriorsToMatrix.BP = function(bpRes, stripConstantParameters=TRUE, stripCatActive=TRUE, stripCatMu=TRUE) {
	
	mat = NULL
	
	for (grp in names(bpRes$groups)) {
		temp = convertPosteriorsToMatrix.WP(bpRes$groups[[ grp ]], stripConstantParameters, stripCatActive, stripCatMu)
		colnames(temp) = paste(grp, ":", colnames(temp), sep="")
		mat = cbind(mat, temp)
	}
	
	mat
}

###############################################################################