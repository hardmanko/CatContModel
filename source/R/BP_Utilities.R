
#' Combine Multiple Within-Participants Groups into a BP Results Object
#' 
#' @param groups A named list in which each element is the return value of [`runParameterEstimation`].
#' 
#' @return An object that combines the results into a form that is usuable by the between-participants and generic functions in this package.
#' 
#' @family BP functions
#'
#' @export
combineGroupResults.BP = function(groups) {
	
	if (length(names(groups)) != length(groups)) {
		stop('"groups" must be a *named* list.')
	}
	
	checkResults.BP(groups)
	
	config = groups[[1]]$config
	config$factors = makeDefaultFactors.BP(groups)
	config$conditionEffects = makeConditionEffects.BP(groups)
	
	runConfig = groups[[1]]$runConfig
	# TODO: Double check that all groups have same iterations
	
	# Things that must go in general config:
	# modelVariant
	# dataType
	
	# Things that really should be the same for all groups:
	# priors in general
	# maxCategories
	# minSD
	# catMuPriorApproximationPrecision
	# conditionEffects                     <--- TODO????
	# responseRange (and studyRange?)
	
	# Things that may be different between groups
	# iterationsPerStatusUpdate

	bpRes = list(groups = groups, config = config, runConfig = runConfig)
	class(bpRes) = c(class(bpRes), "CCM_BP")

	bpRes
}

checkResults.BP = function(groups, iterations = TRUE, modelVariant = TRUE, dataType = TRUE) {
  
  #TODO: Also check priors.
  # Use compareResults to check that it's the same model config, including priors and constant values.
  #for (i in 2:length(groups)) {
  #  compareResults(groups[[1]], groups[[i]], exclude=c("data", "mhTuning"))
  # data=FALSE, mhTuning=FALSE, config=TRUE, constantValues=TRUE, equalityConstraints = TRUE, priors=TRUE
  #}
  
  checkEq = c(
    ifelse(iterations, "iterations", NA),
    ifelse(modelVariant, "modelVariant", NA),
    ifelse(dataType, "dataType", NA)
  )
  checkEq = checkEq[ !is.na(checkEq) ]
  
  for (ch in checkEq) {
    values = NULL
    for (n in names(groups)) {
      values = c(values, groups[[n]]$config[[ch]])
    }
    if (length(unique(values)) > 1) {
      stop(paste0("config$", ch, " is not the same between different results."))
    }
  }
  
}

# BP condition effects are created by taking the union of the CE of all WP groups.
makeConditionEffects.BP = function(groups) {
	
	bpCE = list()
	for (grp in names(groups)) {
		
		config = groups[[ grp ]]$config
		
		ce = config$conditionEffects
		
		varyingFactors = getAllFactorNames(config$factors, removeConstant = TRUE)
		
		for (n in names(ce)) {
			if (ce[[n]] == "all") {
				ce[[n]] = varyingFactors # Only use varying factors
			} else if (ce[[n]] == "none") {
				ce[[n]] = character(0)
			}
			
			# Only if there are varying factors, add in those factor names
			if (length(varyingFactors) > 0) {
				bpCE[[n]] = union(bpCE[[n]], ce[[n]])
			}
		}
		
	}
	
	bpCE
	
}

backPropogateBPFactorsToWPFactors = function(bpRes) {
	
	bpFact = bpRes$config$factors
	
	bpFactNames = getFactorTypeToName(bpFact)$bp
	
	for (grp in names(bpRes$groups)) {
		thisBP = bpFact[ bpFact$group == grp, ]
		thisBP[ , c("key", "group", bpFactNames) ] = NULL
		bpRes$groups[[ grp ]]$config$factors = thisBP
	}
	
	bpRes
}

#backPropogateFactorsToConditionEffects = function(factors) {
	# This function is not possible because factors is for the whole design,
	# but each parameter can have different condition effects.
#}



#The name of the between-participants factor is "group"
#Note that depending on the design, group may need to be split up into more factors
makeDefaultFactors.BP = function(groups) {
	
	factors = NULL
	for (grp in names(groups)) {
		f = groups[[grp]]$config$factors
		f$group = grp
		
		f$key = paste(f$group, f$cond, sep=":")
		
		if (!is.null(factors)) {
			for (n in names(f)) {
				if (!(n %in% names(factors))) {
					factors[, n] = NA
				}
			}
			for (n in names(factors)) {
				if (!(n %in% names(f))) {
					f[, n] = NA
				}
			}
			f = f[ , names(factors) ]
		}

		factors = rbind(factors, f)
	}
	
	if (any(is.na(factors))) {
		logWarning("Factor names appear to be inconsistent between the groups. Did you set up the factors correctly before running parameter estimation? You should examine and correct the factors data frames for the individual groups and rerunning parameter estimation. This is not something that can be reliably fixed after running parameter estimation.")
	}
	
	factors = normalizeFactors(factors, removeConstant=TRUE)
	
	type2name = getFactorTypeToName(factors)
	if (length(groups) > 1 && length(type2name$bp) == 0) {
		factors$BP_Group = factors$group
		logMsg("No between-participants factors found. A default between-participants factor named BP_Group has been created. To remove this factor from the design, set bpRes$config$factors$BP_Goup to NULL.")
	}
	
	allFn = getAllFactorNames(factors)
	
	factors = factors[ , c("key", "group", "cond", allFn) ]
	
	factors
}




getAllPnums.BP = function(bpRes) {
	allPnum = NULL
	for (gName in names(bpRes$groups)) {
		gPnums = bpRes$groups[[gName]]$pnums
		gPnums = paste0(gName, ":", gPnums)
		allPnum = c(allPnum, gPnums)
	}
	allPnum
}


pnumListToVector = function(pnumList) {
	ap = NULL
	for (n in names(pnumList)) {
		ap = c(ap, paste0(n, ":", pnumList[[n]]))
	}
	ap
}

pnumVectorToList = function(pnumV) {
	rval = list()
	parts = strsplit(pnumV, split=":", fixed=TRUE)
	for (i in 1:length(parts)) {
		g = parts[[i]][1]
		rval[[ g ]] = c(rval[[ g ]], parts[[i]][2])
	}
	for (n in names(rval)) {
		rval[[n]] = unique(rval[[n]])
	}
	rval
}

splitKeyVector = function(keys) {
	partList = strsplit(keys, ":", TRUE)
	rval = data.frame(key = keys, group = keys, cond = keys, stringsAsFactors = FALSE)
	for (i in 1:length(partList)) {
		rval$group[i] = partList[[i]][1]
		rval$cond[i] = partList[[i]][2]
	}
	rval
}



getMatchingKeysForUniqueFL_row = function(factors, uniqueFL_row) {
	theseKeys = NULL
	
	for (i in 1:nrow(factors)) {

		#factors[ i, names(uniqueFL_row), drop=FALSE ]
		thisFO = subset(factors, select = names(uniqueFL_row), subset = i == 1:nrow(factors))
		
		if (all(thisFO == uniqueFL_row)) {
			theseKeys = c(theseKeys, factors$key[i])
		}
	}
	
	theseKeys
}

# factors should be normalized
# uniqueFL contains only the factors of interest
getMatchingKeysForUniqueFL = function(factors, uniqueFL) {
	
	#Force only wanted columns
	factors = factors[ , c(names(uniqueFL), "key") ]
	
	allKeys = list()
	
	for (i in 1:nrow(uniqueFL)) {
		
		uniqueFL_row = subset(uniqueFL, subset = i == 1:nrow(uniqueFL))
		
		allKeys[[i]] = getMatchingKeysForUniqueFL_row(factors, uniqueFL_row)
	}
	
	allKeys
}

# factors should be normalized
# factor is scalar character
getMatchingFactorLevels = function(factors, uniqueFL, factor, collapse=NULL) {
	
	keys = getMatchingKeysForUniqueFL(factors, uniqueFL)
	if (length(keys) == 0) {
		return(list())
	}
	
	matchingFL = list()
	
	for (i in 1:length(keys)) {
		matchingFL[[i]] = factors[ factors$key %in% keys[[ i ]], factor ]
	}
	
	if (!is.null(collapse)) {
		rval = rep("", length(matchingFL))
		for (i in 1:length(matchingFL)) {
			rval[i] = paste(sort(unique(matchingFL[[i]])), collapse=collapse)
		}
		matchingFL = rval
	}
	
	matchingFL
	
}





# unused
matchFactorRows = function(factors, urow) {
	
	factors = subset(factors, select = names(urow))
	
	ind = 1:nrow(factors)
	for (i in ind) {
		if (any(subset(factors, subset = i == ind) != urow)) {
			ind[i] = NA
		}
	}
	
	ind[ !is.na(ind) ]
}


defaultGroupName = function() {
	"default"
}


