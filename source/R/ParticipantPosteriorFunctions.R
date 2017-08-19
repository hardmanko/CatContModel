

###############################################################################

#' Vectors of Participant-Level Posterior Parameter Chains
#' 
#' This function helps with retrieving posterior chains for parameters.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param param The name of the parameter.
#' @param pnum The participant number. If `NULL`, the population mean of the parameter will be used. For between-participants designs, may be provided in a `"group:pnum"` format, in which case `group` may be omitted.
#' @param cond The condition to use. If `NULL`, no condition effect will be applied.
#' @param group Only for between-participants designs. The name of the group of participants that `pnum` is in.
#' @param manifest If `TRUE`, the default, the chain will be converted from the latent to the manifest space. If `FALSE`, it will be left in the latent space.
#' 
#' @return A vector containing the parameter chain.
#' 
#' @family generic functions
#' @md
#' @export
getParameterPosterior = function(res, param, pnum, cond, group = NULL, manifest = TRUE) {
	
	if (grepl(":", pnum, fixed=TRUE)) {
		grp_pnum = strsplit(pnum, ":", fixed=TRUE)[[1]]
		group = grp_pnum[1]
		pnum = grp_pnum[2]
	}
	
	if (resultIsType(res, "WP")) {
		if (!is.null(group)) {
			warning('The "group" argument is ignored for within-participants designs.')
		}
		rval = getParameterPosterior.WP(res, param, pnum, cond, manifest = manifest)
	} else if (resultIsType(res, "BP")) {
		if (is.null(group)) {
			stop('The "group" argument must be provided for between-participants designs.')
		}
		rval = getParameterPosterior.BP(res, param, pnum, cond, group = group, manifest = manifest)
	}
	
	rval
	
}

getParameterPosterior.WP = function(results, param, pnum, cond, manifest=TRUE) {
	
	if (!is.null(pnum)) {
		param_base = results$posteriors[[ paste(param, "[", pnum, "]", sep="") ]]
	} else {
		param_base = results$posteriors[[ paste(param, ".mu", sep="") ]]
	}
	
	param_cond = rep(0, length(param_base))
	if (!is.null(cond)) {
		param_cond = results$posteriors[[ paste(param, "_cond[", cond, "]", sep="") ]]
	}
	
	param_latent = param_base + param_cond
	
	transformationFunction = getParameterTransformation(results, param)
	if (!manifest) {
		transformationFunction = function(x) { x }
	}
	transformationFunction(param_latent)
	
}

getParameterPosterior.BP = function(bpRes, param, pnum, cond, group, manifest=TRUE) {
	getParameterPosterior.WP(bpRes$groups[[ group ]], param, pnum, cond, manifest = manifest)
}

###############################################################################

#' Participant-Level Posterior Parameter Chains for all Participants, Conditions, and Groups
#' 
#' The results of this function are essentially those of [`getParameterPosterior`], but for all participants, conditions, and groups.
#' 
#' @family generic functions
#' @export
getAllParameterPosteriors = function(res, param, manifest, format = "data.frame") {
	
	if (resultIsType(res, "WP")) {
		rval = getAllParameterPosteriors.WP(res, param, manifest)
	} else if (resultIsType(res, "BP")) {
		rval = getAllParameterPosteriors.BP(res, param, manifest)
	}
	
	if (format == "matrix") {
		# do nothing
	} else if (format == "data.frame") {
		rval = reshapeMatrixToDF(rval$post, rval$design)
	}
	
	rval
}

getAllParameterPosteriors.BP = function(bpRes, param, manifest = FALSE) {
	
	allPost = NULL
	allDesign = NULL
	
	for (group in names(bpRes$groups)) {
		
		ppce = getAllParameterPosteriors.WP(bpRes$groups[[ group ]], param, manifest = manifest)
		
		ppce$design$group = group
		ppce$design$key = paste(ppce$design$group, ppce$design$cond, sep=":")
		
		allDesign = rbind(allDesign, ppce$design)
		
		allPost = cbind(allPost, ppce$post)
		
	}
	
	list(design = allDesign, post = allPost)
}

getAllParameterPosteriors.WP = function(results, param, manifest = FALSE) {
	
	design = unique(subset(results$data, select=c("pnum", "cond")))
	
	postCombined = matrix(NA, nrow=results$config$iterations, ncol=nrow(design))
	colnames(postCombined) = paste(design$pnum, design$cond, sep=":")
	
	for (i in 1:nrow(design)) {
		
		postCombined[,i] = getParameterPosterior.WP(results, param, pnum = design$pnum[i], cond = design$cond[i], manifest = manifest)
		
	}
	
	list(design = design, post = postCombined)
	
}




###############################################################################

#' Transformed Participant-Level Parameters for a Participant, Condition, and Iteration
#' 
#' Note that inactive categories have already been removed, which is why the return value does not include the `catActive` parameters.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param pnum A single participant number. For between-participants designs, may be provided in a `"group:pnum"` format, in which case `group` may be omitted.
#' @param cond A single condition name.
#' @param iteration The index of an iteration of the chain.
#' @param group Only for between-participants designs. The name of the group of participants that `pnum` is in.
#' @param removeInactiveCategories Whether `catMu` parameters for inactive categories should be removed. If `TRUE`, `catMu` parameters for inactive categories are removed and `catActive` is not provided because it would be 1 for all included categories. If `FALSE`, both `catMu` and `catActive` parameters are provided.
#' 
#' @return A list containing transformed parameters.
#' 
#' @md
#' @family generic functions
#' @export
getSingleIterationParameters = function(res, pnum, cond, iteration, group = NULL, removeInactiveCategories = TRUE) {
	
	if (grepl(":", pnum, fixed=TRUE)) {
		grp_pnum = strsplit(pnum, ":", fixed=TRUE)[[1]]
		group = grp_pnum[1]
		pnum = grp_pnum[2]
	}
	
	if (resultIsType(res, "WP")) {
		if (!is.null(group)) {
			warning('The "group" argument is ignored for within-participants designs.')
		}
		rval = getSingleIterationParameters.WP(res, pnum, cond, iteration = iteration, removeInactiveCategories = removeInactiveCategories)
	} else if (resultIsType(res, "BP")) {
		if (is.null(group)) {
			stop('The "group" argument must be provided for between-participants designs.')
		}
		rval = getSingleIterationParameters.BP(res, pnum, cond, group = group, iteration = iteration, removeInactiveCategories = removeInactiveCategories)
	}
	
	rval
}

getSingleIterationParameters.WP = function(results, pnum, cond, iteration, removeInactiveCategories = TRUE) {
	
	allParam = getAllParams(results, filter=TRUE)
	
	combinedParam = list()
	for (pp in allParam) {
		condParam = results$posteriors[[ paste(pp, "_cond[", cond, "]", sep="") ]][iteration]
		
		partParam = results$posteriors[[ paste(pp, "[", pnum, "]", sep="") ]][iteration]
		
		trans = getParameterTransformation(results, pp)
		
		combinedParam[[pp]] = trans(partParam + condParam)
	}
	
	#Do catMu
	ca = cm = rep(NA, results$config$maxCategories)
	for (i in 1:results$config$maxCategories) {
		istr = paste("[", pnum, ",", i - 1, "]", sep="")
		cm[i] = results$posteriors[[ paste("catMu", istr, sep="") ]][iteration]
		ca[i] = results$posteriors[[ paste("catActive", istr, sep="") ]][iteration]
	}
	
	if (removeInactiveCategories) {
		combinedParam$catMu = cm[ ca == 1 ]
	} else {
		combinedParam$catMu = cm
		combinedParam$catActive = ca
	}
	
	if (results$config$dataType == "circular") {
		combinedParam$catMu = combinedParam$catMu %% 360 #limit to the interval [0, 360)
	}
	
	combinedParam
}

getSingleIterationParameters.BP = function(bpRes, pnum, cond, group, iteration, removeInactiveCategories = TRUE) {
	getSingleIterationParameters.WP(bpRes$groups[[ group ]], pnum, cond, group, iteration, removeInactiveCategories = removeInactiveCategories)
}

###############################################################################

#' Convert Participant Parameters to Matrices or Arrays
#' 
#' For standard parameters (for which each participant has only 1 parameter), the result is a matrix where each row of the matrix is an iteration and each column of the matrix is a participant.
#' 
#' For the category parameters (`catMu` and `catActive`), the result is an array where the 3 dimensions are participant, category, and iteration, in that order.
#' 
#' The column index for a pnum can be gotten with `which(results$pnums == pnum)`.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param param A list of parameter names. If `NULL`, the default, all parameters are done.
#' 
#' @return A named list of matrices and/or arrays.
#' 
#' @md
#' @family generic functions
#' @export
convertPosteriorsToMatrices = function(res, param = NULL) {
	
	if (resultIsType(res, "WP")) {
		rval = convertPosteriorsToMatrices.WP(res, param)
	} else if (resultIsType(res, "BP")) {
		rval = convertPosteriorsToMatrices.BP(res, param)
	}
	
	rval
}

convertPosteriorsToMatrices.WP = function(results, param = NULL) {
	
	matrixParams = c("pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess", 
									 "contSD", "catSelectivity", "catSD")
	arrayParams = c("catActive", "catMu")
	
	if (!is.null(param)) {
		matrixParams = matrixParams[ matrixParams %in% param ]
		arrayParams = arrayParams[ arrayParams %in% param ]
	}
	
	post = list()
	for (mp in matrixParams) {
		post[[mp]] = matrix(0, nrow=results$config$iterations, ncol=length(results$pnums))
		colnames(post[[mp]]) = results$pnums
	}
	for (ap in arrayParams) {
		post[[ap]] = array(0, dim=c(length(results$pnums), results$config$maxCategories, results$config$iterations))
		dimnames(post[[ap]]) = list(results$pnums)
	}
	
	for (pInd in 1:length(results$pnums)) {
		
		pnum = results$pnums[pInd]
		
		istr = paste("[", pnum, "]", sep="")
		
		for (mp in matrixParams) {
			post[[mp]][,pInd] = results$posteriors[[paste(mp, istr, sep="")]]
		}
		
		if (results$config$maxCategories > 0) {
			for (cat in 1:results$config$maxCategories) {
				
				istr = paste("[", pnum, ",", cat - 1, "]", sep="")
				
				for (ap in arrayParams) {
					
					pname = paste(ap, istr, sep="")
					
					x = results$posteriors[[pname]]
					
					post[[ap]][pInd,cat,] = x
					
				}
				
			}
		}
	}
	
	post
}

convertPosteriorsToMatrices.BP = function(bpRes, param=NULL) {
	
	matrixParams = c("pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess", 
									 "contSD", "catSelectivity", "catSD")
	arrayParams = c("catActive", "catMu")
	
	if (is.null(param)) {
		param = c(matrixParams, arrayParams)
	}
	
	post = list()
	
	for (n in names(bpRes$groups)) {
		
		results = bpRes$groups[[n]]
		tp = convertPosteriorsToMatrices.WP(results, param)
		
		for (mp in matrixParams[ matrixParams %in% param ]) {
			colnames(tp[[mp]]) = paste0(n, ":", colnames(tp[[mp]]))
			post[[mp]] = cbind(post[[mp]], tp[[mp]])
		}
		for (ap in arrayParams[ arrayParams %in% param ]) {
			dn = dimnames(tp[[ap]])
			dn[[1]] = paste0(n, ":", dn[[1]])
			dimnames(tp[[ap]]) = dn
			
			a = abind::abind(post[[ap]], tp[[ap]], along = 1)
			
			post[[ap]] = a
		}
	}
	
	post
}

