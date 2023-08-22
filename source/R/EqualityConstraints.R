
# Equality constraints are fundamentally a WP concept, so there are no BP functions
# related to equality constraints.



getConstrainedConditionEffects = function(config) {

	cec = NULL
	
	for (parName in names(config$conditionEffects)) {
		
		thisFactors = getFactorsForConditionEffect.WP(config, parName)

		thisCEC = getConstrainedConditionEffectList(config=config, parName=parName, usedFactors=thisFactors)
		cec = rbind(cec, thisCEC)
		
	}
	
	if (length(unique(cec$target)) != nrow(cec)) {
		#TODO: Remove duplicates?
		#TODO: Put this elsewhere?
		warning("Duplicate equality constraints.")
	}
	
	#strip non-condition constraints
	#TODO: This really does not belong here
	rejectRows = NULL
	for (i in 1:nrow(cec)) {
		target = cec$target[i]
		source = cec$source[i]
		
		rejectTarget = target != "FREE_PARAMETER" && !grepl("_cond[", target, fixed = TRUE)
		rejectSource = source != "FREE_PARAMETER" && !grepl("_cond[", source, fixed = TRUE)
		
		if (rejectTarget || rejectSource) {
			rejectRows = c(rejectRows, i)
		}
	}
	cec = cec[ !(1:nrow(cec) %in% rejectRows), ]
	
	#convert to list
	rval = list()
	for (i in 1:nrow(cec)) {
		rval[[ cec$target[i] ]] = cec$source[i]
	}
	
	rval
}


getConstrainedConditionEffectList = function(config, parName, usedFactors) {
	
	factors = config$factors
	
	rval = NULL
	
	#special case: no factors used.
	if (length(usedFactors) == 0) {
		for (i in 1:nrow(factors)) {
			target = paste(parName, "_cond[", factors$cond[i], "]", sep="")
			source = paste(parName, "_cond[", config$cornerstoneConditionName, "]", sep="")
			
			if (factors$cond[i] == config$cornerstoneConditionName) {
				source = "FREE_PARAMETER"
			}
			
			rval = rbind(rval, data.frame(target=target, source=source, stringsAsFactors=FALSE))
		}
		return(rval)
	}
	
	
	cefcond = factors[ , c(usedFactors, "cond") ]
	unicef = unique(factors[ , usedFactors ]) #get unique combinations of free condition effects
	if (length(usedFactors) == 1) {
		unicef = data.frame(a=unique(factors[ , usedFactors ]))
		names(unicef) = usedFactors
	}
	
	for (i in 1:nrow(unicef)) {
		
		rows = NULL
		for (j in 1:nrow(cefcond)) {
			match = all(unicef[i, usedFactors] == cefcond[j, usedFactors])
			if (match) {
				rows = c(rows, j)
			}
		}
		
		theseConds = cefcond[rows,"cond"]
		
		#select which condition will be used as the source parameter for this factor level
		thisLevelSource = theseConds[1]
		if (config$cornerstoneConditionName %in% theseConds) {
			thisLevelSource = config$cornerstoneConditionName
		}
		
		for (cond in theseConds) {
			target = paste(parName, "_cond[", cond, "]", sep="")
			
			if (cond == thisLevelSource) {
				source = "FREE_PARAMETER"
			} else {
				source = paste(parName, "_cond[", thisLevelSource, "]", sep="")
			}
			
			rval = rbind(rval, data.frame(target=target, source=source, stringsAsFactors=FALSE))
		}
		
	}

	rval
	
}


getConditionParameterParts = function(parName) {
	
	if (!grepl("_cond[", parName, fixed=TRUE)) {
		warning("Non-condition parameter provided to getConditionParameterParts().")
		return("")
	}
	
	bits = strsplit(parName, "_cond[", fixed=TRUE)[[1]]
	
	rval = list(parName = bits[1])
	
	rval$cond = strsplit(bits[2], "]", fixed=TRUE)[[1]]
	
	rval
}


getRootSourceConditionParameter = function(results, parName, cond, fullParam = NULL) {
	
	target = paste(parName, "_cond[", cond, "]", sep="")
	if (!is.null(fullParam)) {
		target = fullParam
	}
	
	source = results$equalityConstraints[[ target ]]
	
	if (is.null(source)) {
		warning(paste("No source parameter specified for \"", target, "\". It must be a free parameter.", sep=""))
		source = "FREE_PARAMETER"
	}
	
	rval = ""
	
	if (source == "FREE_PARAMETER") {
		rval = target
	} else {
		rval = getRootSourceConditionParameter(results, "", "", fullParam = source)
	}
	
	rval
	
}

# WP only
getEqualConditionParameters = function(results, parName) {
	
  if (!resultIsType(results, "WP")) {
    stop("getEqualConditionParameters can only be used with WP results objects.")
  }
  
	fcopy = results$config$factors
	fcopy$equalityGroup = 0
	
	for (i in 1:nrow(fcopy)) {
		target = paste(parName, "_cond[", fcopy$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, parName, fcopy$cond[i])
		
		if (target == rootSource) {
			fcopy$equalityGroup[i] = max(fcopy$equalityGroup) + 1
		}
		
	}
	
	for (i in 1:nrow(fcopy)) {
		target = paste(parName, "_cond[", fcopy$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, parName, fcopy$cond[i])
		
		if (target != rootSource) {
			
			sourceCond = getConditionParameterParts(rootSource)$cond
			
			fcopy$equalityGroup[i] = fcopy$equalityGroup[ fcopy$cond == sourceCond ]
		}
		
	}
	
	fcopy[ , c("cond", "equalityGroup") ]
}


# This function accounts for equality constraints
getConditionParameterPrior = function(results, parName, cond, fullParam = NULL) {
	
	target = paste(parName, "_cond[", cond, "]", sep="")
	if (!is.null(fullParam)) {
		target = fullParam
	}
	
	constVal = results$constantValueOverrides[[ target ]]
	if (!is.null(constVal)) {
		rval = list(location = constVal, scale = 0)
		return(rval)
	}
	
	source = results$equalityConstraints[[ target ]]
	#rootSource = getRootSourceParameter(results, parName, cond)
	
	if (is.null(source)) {
		warning(paste("No source parameter specified for \"", target, "\". It must be a free parameter.", sep=""))
		source = "FREE_PARAMETER"
	}
	
	rval = list()
	
	if (source == "FREE_PARAMETER") {
		
		targetParts = getConditionParameterParts(target)
		
		if (targetParts$cond == results$config$cornerstoneConditionName) {
			
			rval$location = 0
			rval$scale = 0
			
		} else {
			
			rval$location = results$priors[[ paste(targetParts$parName, "_cond.loc", sep="") ]]
			rval$scale = results$priors[[ paste(targetParts$parName, "_cond.scale", sep="") ]]
			
		}
		
	} else {
		rval = getConditionParameterPrior(results, "", "", fullParam = source)
	}
	
	rval
}


getParametersWithConditionEffects = function(conditionEffects) {
	rval = character(0)
	for (n in names(conditionEffects)) {
		if (all(conditionEffects[[n]] != "none")) {
			rval = c(rval, n)
		}
	}
	rval
}

