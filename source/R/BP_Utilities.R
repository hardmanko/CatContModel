
#' Merge Multiple Within-Participants Groups
#' 
#' @param groups A NAMED list in which each element is the return value of [`runParameterEstimation`].
#' 
#' @return An object that merges the results into a form that is usuable by the between-participants functions in this package.
#' 
#' @family BP functions
#' @md
#' @export
mergeGroupResults.BP = function(groups) {
	
	if (length(names(groups)) != length(groups)) {
		stop("groups must be a *named* list.")
	}
	
	CatContModel:::checkResults.BP(groups)
	
	# This needs a bunch of careful thought still, 
	# but might actually be good enough to work.
	config = groups[[1]]$config
	config$factors = makeDefaultFactors.BP(groups)
	config$conditionEffects = makeConditionEffects.BP(groups)
	
	# Things that must go in general config:
	# modelVariant
	# iterations
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

	bpRes = list(groups = groups, config = config)
	class(bpRes) = c(class(bpRes), "CCM_BP")

	bpRes
}

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
	
	#for (n in names(bpCE)) {
	#	if (length(bpCE[[n]]) == 0) {
	#		bpCE[[n]] = "none"
	#	}
	#}
	
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

#TODO: Also check priors
checkResults.BP = function(groups, iterations = TRUE, modelVariant = TRUE, dataType = TRUE) {
	
	# Consider using the new compareResults function
	#compareResults(groups[[1]], groups[[2]], data=FALSE, mhTuning=FALSE, constantValueOverrides=FALSE, equalityConstraints = FALSE, priors=TRUE?, config=TRUE, configIgnore = c("iterationsPerStatusUpdate"))
	
	
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


#The name of the between-participants factor is "group"
#Note that depending on the design, group may need to be split up into more factors
makeDefaultFactors.BP = function(groups) {
	
	factors = NULL
	for (grp in names(groups)) {
		f = groups[[grp]]$config$factors
		f$group = grp
		f$BP_Factor = grp
		
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
		warning("Factor names appear to be inconsistent between the groups. Did you set up the factors correctly before running parameter estimation? That is the best approach to use. You must manually examine and correct the factors data frame to correctly represent the design, possibly also doing so for each group individually.")
	}
	
	allFn = getAllFactorNames(factors)
	
	factors = factors[ , c("key", "group", "cond", allFn) ]
	
	factors
}




getAllPnums.BP = function(bpRes) {
	allPnum = NULL
	for (n in names(bpRes$groups)) {
		pn = bpRes$groups[[n]]$pnums
		pn = paste0(n, ":", pn)
		allPnum = c(allPnum, pn)
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

getMatchingKeysForUniqueFL_row = function(foFactors, uniqueFL_row) {
	theseKeys = NULL
	
	for (i in 1:nrow(foFactors)) {
		
		thisFO = subset(foFactors, select = names(uniqueFL_row), subset = i == 1:nrow(foFactors))
		
		if (all(thisFO == uniqueFL_row)) {
			theseKeys = c(theseKeys, foFactors$key[i])
		}
	}
	
	theseKeys
}

# foFactors contains only the factors of interest and "key"
# uniqueFL contains only the factors of interest
getMatchingKeysForUniqueFL = function(foFactors, uniqueFL) {
	
	#Force only wanted columns
	foFactors = foFactors[ , c(names(uniqueFL), "key") ]
	
	allKeys = list()
	
	for (i in 1:nrow(uniqueFL)) {
		
		uniqueFL_row = subset(uniqueFL, subset = i == 1:nrow(uniqueFL))
		
		allKeys[[i]] = getMatchingKeysForUniqueFL_row(foFactors, uniqueFL_row)
	}
	
	allKeys
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


