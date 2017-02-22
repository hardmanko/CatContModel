
#The functions in this file are eventually going into their own package, which is why they
#do not share interface conventions with the rest of CatContModel.
#
#One important constraint is that "factors" in this file is not allowed to have any columns
#other than the factors (NO "cond" COLUMN!).


# Sort of internal function.
# only appropriate if x has length 2.
devianceFunction_absDif = function(x) {
	abs(x[1] - x[2])
}


# Sort of internal function. Could be visible to users. If so, should be renamed.
# prior and posterior are each vectors
# it is assumed that the prior may be diffuse, such as a Cauchy prior.
# In that case, it is possible to have prior values that are very far from the
# largest posterior value (like 10^9 farther). This hurts the ability of the
# density estimation to estimate the density in the same way for the prior and
# posterior, which could result in a bias.
testFunction_savageDickey_diffusePrior = function(prior, posterior) {
	
	# If the prior extends to very large values relative to the posterior,
	# cut it off. This helps with density estimation at smaller values.
	priorMax = min( 100 * max(posterior), max(prior) * 1.1 )
	keptPriorValues = (prior <= priorMax)
	pKept = mean(keptPriorValues)
	prior = prior[ keptPriorValues ]
	
	priorLS = polspline::logspline(prior, lbound = 0, ubound = priorMax)
	priorDens = polspline::dlogspline(0, priorLS)
	
	# Account for the fact that there is some density in the upper area 
	# that isn't accounted for by the logspline. Scale down the density
	# in the lower area because the upper area was not included and the 
	# lower area thus had too much data in it.
	priorDens = priorDens * pKept
	
	#Do the posterior in a tryCatch because failures can happen.
	success = tryCatch({
		postLS = polspline::logspline(posterior, lbound = 0)
		TRUE
	}, error = function(e) {
		print(e)
		return(FALSE)
	})
	
	if (success) {
		postDens = polspline::dlogspline(0, postLS)
		bf10 = priorDens / postDens
	} else {
		bf10 = NA
	}
	
	list(bf01 = 1 / bf10, bf10 = bf10, success = success)
}


splitCellName = function(str) {
	parts = strsplit(str, split = ":", fixed=TRUE)[[1]]
	
	fl = list()
	
	for (i in 1:length(parts)) {
		
		parts2 = strsplit(parts[i], split=".", fixed=TRUE)[[1]]
		
		fact = parts2[1]
		lev = parts2[2]
		
		fl[[ fact ]] = lev
		
	}
	fl
}


#cnl can be a list or data.frame.
makeCellName = function(cnl) {
	n = ""
	fNames = names(cnl)
	for (i in 1:length(fNames)) {
		fact = fNames[i]
		lev = cnl[[fact]]
		n = paste(n, fact, ".", lev, sep="")
		if (i < length(fNames)) {
			n = paste(n, ":", sep="")
		}
	}
	n
}


# Maybe internal function?
# Renaming the columns doesn't really make sense, because the terms are the terms
# as named as used by the contrasts function. As such, there is no natural renaming,
# other than for contr.treatment (and contr.SAS).
makeDesignMatrix = function(factors, dmFactors, contrastType, renameCols=FALSE) {
	
	form = paste("~", paste(dmFactors, collapse = " * "))
	form = stats::formula(form)
	
	tt = stats::terms(form)
	
	contrasts = list()
	for (fn in dmFactors) {
		contrasts[[fn]] = contrastType
	}
	
	m = stats::model.matrix(form, factors, contrasts)
	
	if (renameCols && (contrastType == "contr.treatment" || contrastType == "contr.SAS")) {
		colnames(m) = renameDesignMatrixColumns( colnames(m), factors, dmFactors, contrastType)
	}
	
	list(mat=m, terms=tt)
}


# Internal function.
# If the design matrix is not full rank, this function removes columns 
# that contribute to it being not full rank.
stripExcessTermsFromDM = function(dm) {
	
	mat = dm$mat
	
	q = qr(mat)
	
	kept = q$pivot[ seq(q$rank) ]
	dropped = q$pivot[ !(seq(ncol(q$qr)) %in% seq(q$rank)) ]

	newMat = mat[ , kept ]
	
	attr(newMat, "assign") = attr(mat, "assign")[ kept ]
	attr(newMat, "contrasts") = attr(mat, "contrasts")
	
	list(mat=newMat, kept = colnames(mat)[kept], dropped = colnames(mat)[dropped])
	
}


# Internal function.
# Gets the columns of the design matrix associated with particular effects.
# Those columns are the same as the elements of the effect parameters.
getEffectAssignmentColumns = function(dm, fNames) {
	
	labels = attr(dm$terms, "term.labels")
	
	#Select first order (main effect) terms.
	o1labels = labels[ attr(dm$terms, "order") == 1 ]
	
	if (!all(fNames %in% o1labels)) {
		stop("fNames and design matrix term labels do not match.")
	}
	
	fNames = o1labels[ o1labels %in% fNames ] # sort fNames by label names
	
	effect = paste0(fNames, collapse=":")
	
	assignments = attr(dm$mat, "assign")
	
	effectAssignment = which(labels == effect)
	mCols = which(assignments == effectAssignment)
	
	mCols
	
}


# Maybe make external?
isDesignFullyCrossed = function(factors, warnOnDuplicate = TRUE) {
	
	totalCells = 1
	for (n in names(factors)) {
		totalCells = totalCells * length(unique(factors[,n]))
	}
	
	uf = unique(factors)
	if (warnOnDuplicate && nrow(uf) != nrow(factors)) {
		warning("There are duplicate rows in factors. This may be ok.")
	}
	
	factors = uf
	if (nrow(factors) == totalCells) {
		return(TRUE)
	}
	
	FALSE
}


#API function
#
# 1. Calculate the design matrix, X.
# 1b. If design is not fully crossed, strip excess terms from X.
# 2. Calculate S = (X'X)^-1 X'
# 3. Calculate beta = S %*% mu
# 4. Select out only the needed beta, beta_s
# 5. Select the corresponding columns (and rows) in X, X_s
# 6. Complete the betas by X_s %*% beta_s
#
# mu is matrix where each column is one cell mean and each row is an iteration, each column corresponding to a row in factors
# testedFactors: Names of the factors to perform a test of. If length > 1, an interaction is tested.
getEffectParameters = function(cellMeans, factors, testedFactors, dmFactors = testedFactors, 
															 contrastType = NULL, warnOnDrop=FALSE) {
	
	if (is.null(contrastType)) {
		if (isDesignFullyCrossed(factors, warnOnDuplicate=FALSE)) {
			contrastType = "contr.sum"
		} else {
			contrastType = "contr.treatment"
		}
	}
	
	# 1. Calculate the design matrix, X.
	dm = makeDesignMatrix(factors, dmFactors, contrastType, renameCols=FALSE)

	# 1b. If design is not fully crossed, strip excess terms from X.
	strippedInfo = stripExcessTermsFromDM(dm)
	
	if (warnOnDrop && length(strippedInfo$dropped) > 0) {
		warning( paste0("The design is not fully crossed (unbalanced). As a result, some terms were dropped from the design matrix: ", paste(strippedInfo$dropped, collapse=", "), ". The following terms were kept: ", paste(strippedInfo$kept, collapse=", "), ". The naming of these terms depends on the contrastType that was used.") )
	}
	
	X = dm$mat = strippedInfo$mat

	
	# 2. Calculate S = (X'X)^-1 X'
	# If mu = X beta
	# then beta = (X'X)^-1 X'mu = S mu
	S = solve(t(X) %*% X) %*% t(X)
	
	# 3. Calculate beta = S %*% mu
	# Some extra transposition because mu and beta are both transposed
	beta = t(S %*% t(cellMeans))
	
	# 4. Select out only the needed beta, beta_s
	mCols = getEffectAssignmentColumns(dm, testedFactors)
	
	beta_s = subset(beta, select = mCols)
	
	# 5. Select the corresponding columns (and rows) in X, X_s
	X_s = subset(X, select = mCols)
	
	#select only some rows of X_s
	uniqueFL = unique(subset(factors, select = testedFactors))
	X_s_rows = rep(NA, nrow(uniqueFL))
	
	for (i in 1:nrow(uniqueFL)) {
		
		rows = NULL
		for (j in 1:nrow(factors)) {
			if (all(factors[j, testedFactors] == uniqueFL[i,])) {
				rows = c(rows, j)
			}
		}
		X_s_rows[i] = rows[1] #for duplicate rows, the DM elements are the same
	}
	
	X_s = X_s[ X_s_rows, ]
	
	# 6. Complete the betas by X_s %*% beta_s
	fullEffects = beta_s %*% t(X_s)
	
	#Name the columns with cell names
	colnames(fullEffects) = makeCellName(uniqueFL)
	
	fullEffects
}


# Perform Hypothesis Test from Effect Parameters
# 
# Rather than using \code{\link{getEffectParameters}} to calculate effect parameters from cell means,
# you may have directly estimated the effect parameters in your model. In that case, you can use
# this function to test hypotheses related to those effect parameters.
# 
# If using the default Savage-Dickey test function, you should use the same number of samples from the prior and posterior. That is to say, \code{priorEffects} and \code{postEffects} should have the same number of rows. A warning wil be emitted if the number of rows is not equal.
# 
# @param priorEffects Numeric matrix of effect parameters sampled from the priors. Each column is one parameter and each row is one iteration. It must have at least two columns.
# @param postEffects Numeric matrix of effect parameters sampled from the posteriors. Must have the same number of columns as priorEffects. Mathematically, it works even with different numbers of columns, but such a test would be bizarre and meaningless (some parameters didn't have priors?). If the number of columns is not the same, an error will be emitted.
# @param devianceFunction A function used for calculating the deviation of the effect parameters. It takes a vector of effect parameters and calculates some measure of how dispersed they are. One example of such a function is the built-in R function \code{var}.
# @param testFunction A function of two numeric vector arguments. Each element in a vector is a measure of the dispersion, as calculated by \code{devianceFunction}, of the effect parameters. The first argument is from the prior and the second is from the posterior. It performs a hypothesis test of whether the prior and posterior deviations differ. It can return anything as its result is returned directly from this function. See the source code of \code{testFunction_savageDickey_diffusePrior} for an example.
# 
# @return The return value depends on the choice of \code{testFunction}. See \code{\link{testFunction_savageDickey_diffusePrior}} for an example.
testHypothesis_effect = function(priorEffects, postEffects, devianceFunction = NULL, testFunction = NULL) 
{
	
	if (ncol(priorEffects) != ncol(postEffects)) {
		stop("The number of columns (parameters) in priorEffects does not equal the number of columns in postEffects. See the documentation for this function. This is probably the result of a serious conceptual error.")
	}
	
	if (ncol(priorEffects) < 2 || ncol(postEffects) < 2) {
		stop("priorEffects or postEffects has less than two columns. You must have at least two effect parameters to perform a hypothesis test of main effects or interactions.")
	}
	
	if (is.null(devianceFunction)) {
		nc = ncol(priorEffects)
		if (nc == 2) {
			devianceFunction = devianceFunction_absDif
		} else {
			devianceFunction = stats::var #yeah, var
		}
	}
	
	if (is.null(testFunction)) {
		testFunction = testFunction_savageDickey_diffusePrior
	}
	
	priorDeviance = apply(priorEffects, 1, devianceFunction)
	postDeviance = apply(postEffects, 1, devianceFunction)
	
	testFunction(priorDeviance, postDeviance)
}


# Perform Hypothesis Test from Cell Means
# 
# @param priorCMs Numeric matrix. Cell means sampled from the priors. The columns must correspond to the rows of factors but do not need to be named.
# @param postCMs Numeric matrix. Cell means sampled from the posterior distribution. The columns must correspond to the rows of factors but do not need to be named. 
# @param factors A \code{data.frame} containing information about the experimental design. Each column is a factor of the design. Each row contains the levels of the factors that define a cell of the design. No additional columns may be included in factors.
# @param testedFactors Character vector. The factors for which to perform the hypothesis test as a vector of factor names. A single factor name results in the test of the main effect of the factor. Multiple factor names result in the test of the interaction of all of those factors.
# @param dmFactors Character vector. The factors to use to construct the design matrix. For a fully-crossed (balanced) design, this can always be equal to \code{testFactors} (the default). For non-fully-crossed designs, you may sometimes want to create a design matrix using some factors, but perform a hypothesis test with only some of those factors (\code{testedFactors} must be a subset of \code{dmFactors}).
# @param contrastType Character (or function). The contrast to use to create the design matrix. Can be any of the function names on the documentation page for \code{contr.sum}. For a non-fully-crossed (unbalanced) design, you should use either "contr.treatment" or "contr.SAS". For a balanced design, you can use anything, but psychologists are most used to "contr.sum", which uses sums-to-zero constraints.
# @param devianceFunction See the documentation for \code{testHypothesis_effect}. You can probably leave this at the default value.
# @param testFunction See the documentation for \code{testHypothesis_effect}. You can probably leave this at the default value.
testHypothesis = function(priorCMs, postCMs, factors, testedFactors, dmFactors = testedFactors,
													contrastType = NULL, devianceFunction = NULL, testFunction = NULL) {
	
	if (!all(testedFactors %in% dmFactors)) {
		stop("testedFactors must be a subset of dmFactors.")
	}
	
	priorBetas = getEffectParameters(priorCMs, factors, 
																	 testedFactors=testedFactors, dmFactors=dmFactors, contrastType=contrastType)
	postBetas = getEffectParameters(postCMs, factors, 
																	testedFactors=testedFactors, dmFactors=dmFactors, contrastType=contrastType)
	
	ht = testHypothesis_effect(priorBetas, postBetas)
	
	ht
}


# Summary Statistics of Effect Parameters
# 
# Calculates summary statistics, like mean, median, and credible interval, of effect parameters.
# 
# @param effects A matrix, such as from \code{\link{getEffectParameters}}, where each column is one parameter and each row is an iteration. If the columns are named with the name of the effect, those names will be used.
# @param fList A named list of functions of one vector argument that will be applied to the effects individually.
# @param quantiles Quantiles to be calculated. By default, the median and 95% credible interval is calculated.
effectSummary = function(effects, fList = list(mean=mean), quantiles = c(0.025, 0.5, 0.975)) {
	
	if (is.null(colnames(effects))) {
		colnames(effects) = paste("effect", 1:ncol(effects), sep="_")
	}
	
	res = NULL
	for (i in 1:ncol(effects)) {
		
		ef = effects[,i]
		
		temp = data.frame(effect = colnames(effects)[i])
		
		for (n in names(fList)) {
			temp[ , n ] = fList[[n]](ef)
		}
		
		qs = stats::quantile(ef, quantiles )
		for (j in 1:length(qs)) {
			temp[,names(qs)[j]] = qs[j]
		}
		
		res = rbind(res, temp)
		
	}
	
	res
	
}

# Unwanted function
#This function may not be needed (probably shouldn't be needed).
#fNames should be enough fNames to account for everything in oldNames
renameDesignMatrixColumns = function(oldNames, factors, fNames, contrastType) {
	
	newNames = rep(NA, length(oldNames))
	
	for (i in 1:length(oldNames)) {
		
		if (oldNames[i] == "(Intercept)") {
			newNames[i] = oldNames[i]
			next
		}
		
		parts = strsplit(oldNames[i], split = ":", fixed=TRUE)[[1]]
		
		factorsLevels = list()
		
		for (j in 1:length(parts)) {
			dlm = ":"
			
			#Find the factor that starts the name of this part of this term
			for (fn in fNames) {
				if (substr(parts[j], 1, nchar(fn)) == fn) {
					break
				}
			}
			
			level_something = gsub(fn, "", parts[j])
			
			if (contrastType == "contr.sum" || contrastType == "contr.helmert") {
				levInd = as.integer(level_something)
				
				lev = unique(factors[, fn])[levInd]
			} else if (contrastType == "contr.treatment" || contrastType == "contr.SAS") {
				lev = level_something
			} else {
				stop("Invalid contrast type.")
			}
			
			factorsLevels[[fn]] = lev
		}
		
		newNames[i] = makeCellName(factorsLevels)
		
	}
	
	newNames
}

# Unwanted function
# assumes sums to 0 and maybe fully crossed design
# unfilled values should be NA
# can be used with betas or weights
fillOutFullValues_sumsToZero = function(fv, uniqueFL) {
	for (col in colnames(fv)) {
		if (all(!is.na(fv[, col]))) {
			next #if this one is known, do nothing
		}
		
		fl = splitCellName(col)
		for (variedF in names(fl)) {
			variedLevels = unique(uniqueFL[,variedF])
			variedLevels = variedLevels[ variedLevels != fl[[variedF]] ]
			
			othersKnown = rep(FALSE, length(variedLevels))
			otherCols = rep("", length(variedLevels))
			
			for (vli in 1:length(variedLevels)) {
				flc = fl
				flc[[variedF]] = variedLevels[vli]
				variedCol = makeCellName(flc)
				
				othersKnown[vli] = all(!is.na(fv[, variedCol]))
				otherCols[vli] = variedCol
				
			}
			
			if (all(othersKnown)) {
				#we're in business
				
				otherColVals = fv[ , otherCols ]
				if (length(otherCols) == 1) {
					fv[ , col ] = -otherColVals
				} else {
					fv[ , col ] = apply(otherColVals, 1, function(x) { -sum(x) })
				}
				
			}
			
		}
		
	}
	
	colIsNA = apply(fv, 2, function(x) {any(is.na(x))})
	if (any(colIsNA)) {
		fv = fillOutFullValues_sumsToZero(fv, uniqueFL)
	}
	
	fv
}



# Unwanted function
#assumes that partialWeights has proper names
fillInWeights = function(partialWeights, factors, fNames, contrastType, uniqueFL = NULL) {
	
	if (is.null(uniqueFL)) {
		uniqueFL = unique( subset(factors, select = fNames) )
	}
	
	# This does not assume fully crossed, but if the design is fully, crossed, it should work.
	fullWeights = matrix(NA, nrow=nrow(partialWeights), ncol=nrow(uniqueFL))
	
	# Set column names for fullWeights
	# This isn't really right, because the names of partialWeights may not be cell names
	allColNames = rep("", ncol(fullWeights))
	for (i in 1:nrow(uniqueFL)) {
		temp = subset(uniqueFL, subset = (i == 1:nrow(uniqueFL)) )
		allColNames[i] = makeCellName(as.list(temp))
	}
	colnames(fullWeights) = allColNames
	
	for (colname in colnames(partialWeights)) {
		fullWeights[ , colname ] = partialWeights[ , colname ]
	}
	
	if (contrastType == "contr.sum") {
		
		if (!isDesignFullyCrossed(factors, warnOnDuplicate=FALSE)) {
			stop("Design is not fully crossed: You can't use sums-to-zero contrasts for this purpose. Use contr.treatment instead.")
		}
		
		fullWeights = fillOutFullValues_sumsToZero(fullWeights, uniqueFL)
	} else if (contrastType == "contr.treatment") {
		fullWeights[ is.na(fullWeights) ] = 0
	} else {
		stop("Unsupported contrastType.")
	}
	
	fullWeights
	
}

# Unwanted function
getPartialFilledS = function(factors, testedFactors, dmFactors = testedFactors, contrastType = NULL, warnOnDrop = FALSE) {
	
	################################################
	# This section is a C/P from getEffectParameters
	
	if (is.null(contrastType)) {
		if (isDesignFullyCrossed(factors, warnOnDuplicate=FALSE)) {
			contrastType = "contr.sum"
		} else {
			contrastType = "contr.treatment"
		}
	}
	
	# 1. Calculate the design matrix, X.
	dm = makeDesignMatrix(factors, dmFactors, contrastType, renameCols=FALSE)
	
	# 1b. If design is not fully crossed, strip excess terms from X.
	strippedInfo = stripExcessTermsFromDM(dm)
	
	if (warnOnDrop && length(strippedInfo$dropped) > 0) {
		warning( paste0("The design is not fully crossed (unbalanced). As a result, some terms were dropped from the design matrix: ", paste(strippedInfo$dropped, collapse=", "), ". The following terms were kept: ", paste(strippedInfo$kept, collapse=", "), ". The naming of these terms depends on the contrastType that was used.") )
	}
	
	X = dm$mat = strippedInfo$mat
	
	
	# 2. Calculate S = (X'X)^-1 X'
	# If mu = X beta
	# then beta = (X'X)^-1 X'mu = S mu
	S = solve(t(X) %*% X) %*% t(X)
	
	# End C/P
	###############################
	
	
	S = t(S)
	
	mCols = getEffectAssignmentColumns(dm, testedFactors)
	
	S_s = subset(S, select=mCols)
	
	colnames(S_s) = renameDesignMatrixColumns(colnames(S_s), factors, testedFactors, contrastType = contrastType)
	
	full_s = fillInWeights(S_s, factors, testedFactors, contrastType = contrastType)
	
	full_s
}






