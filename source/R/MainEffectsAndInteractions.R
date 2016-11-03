





guessFactorNames = function(factors) {
	n = names(factors)
	n[ n != "cond" ]
}


#appropriate only when v sums to 0
devianceFunction_vectorLength = function(v) {
	sqrt(sum(v^2))
}


#only appropriate if x has length 2
devianceFunction_absDif = function(x) {
	abs(x[1] - x[2])
}

#always appropriate
devianceFunction_absSum = function(x) {
	sum(abs(x))
}

devianceFunction_sd = function(x) {
	stats::sd(x)
}

devianceFunction_var = function(x) {
	#var(x)
	#or
	mean(x^2) # as long as sum(x) == 0
}


testFunction_savageDickey = function(prior, posterior) {
	
	priorMax = min( 100 * max(posterior), max(prior) * 1.1 )
	pKept = mean(prior < priorMax)
	prior = prior[ prior < priorMax ]
	
	
	priorLS = polspline::logspline(prior, lbound = 0, ubound = priorMax)
	priorDens = polspline::dlogspline(0, priorLS)
	
	#account for the fact that there is some density in the upper area that isn't accounted for by the logspline
	priorDens = priorDens * pKept
	
	
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

# The interval test works by finding the 1st percentile of
# both the prior and posterior distributions. An interval from 0 to that
# quantile defines the null hypothesis interval. 
#Tests that the parameter is within an interval from 0 to
#a quantile defined by p.
#
# This doesn't work when the prior is diffuse because
# the entire posterior can be below the 1st percentile of the
# prior. Thus, the posterior is nearer to 0 than the prior, but
# the posterior is far from 0.
testFunction_interval = function(prior, posterior, p = 0.01) {
	
	qprior = stats::quantile(prior, p)
	qpost = stats::quantile(posterior, p)
	
	larger = max(c(qprior, qpost))
	
	pbPrior = mean(prior < larger)
	pbPost = mean(posterior < larger)
	
	bfBelow = ((1 - pbPrior) / pbPrior) * (pbPost / (1 - pbPost))
	
	bfAbove = 1 / bfBelow
	
	list(bf10 = bfAbove, bf01 = bfBelow, p = p, q = larger, success=TRUE)
}



sampleFromConditionEffectPriors = function(results, factors, param, priorSamples) {
	
	priorDs = matrix(NA, nrow=priorSamples, ncol=nrow(factors))
	
	#first pass: Sample from free parameters and the cornerstone condition
	for (i in 1:nrow(factors)) {
		
		target = paste(param, "_cond[", factors$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, param, factors$cond[i])
		
		#If this is its own source (i.e. it is a free parameter
		if (target == rootSource) {
			prior = getConditionParameterPrior(results, param, factors$cond[i])
			if (prior$scale == 0) {
				priorDs[,i] = prior$location
			} else {
				priorDs[,i] = stats::rcauchy(priorSamples, prior$location, prior$scale)
			}
		}
		
	}
	
	#second pass: Copy equality constrained parameters
	for (i in 1:nrow(factors)) {
		
		target = paste(param, "_cond[", factors$cond[i], "]", sep="")
		rootSource = getRootSourceConditionParameter(results, param, factors$cond[i])
		
		if (target != rootSource) {
			
			parts = getConditionParameterParts(rootSource)
			
			sourceCol = which(factors$cond == parts$cond)
			priorDs[,i] = priorDs[,sourceCol]
			
		}
	}
	
	colnames(priorDs) = factors$cond
	
	priorDs
	
}


#only works for fully-crossed designs
#fNames is a character vector of factor names
#lNames is a character vector of the same length as fNames of names of levels in the factors
getFactorByLevelWeights = function(factors, fNames, lNames) {

	#Do a little setup
	combinations = list(GRAND_MEAN = character(0))
	
	combinationLayers = data.frame(combination="GRAND_MEAN", layer=0, stringsAsFactors=FALSE)
	
	#get all combinations at layers less than (or equal to) the current test
	for (layer in 1:length(fNames)) {
		r = utils::combn(fNames, layer)
		r = t(r)
		for (j in 1:nrow(r)) {
			x = sort(r[j,])
			y = paste(x, collapse=":")
			combinations[[y]] = x
			
			temp = data.frame(combination=y, layer=length(x), stringsAsFactors=FALSE)
			
			combinationLayers = rbind(combinationLayers, temp)
		}
	}
	
	
	getWeights_internal = function(fcopy, fname, lname) {
		
		fcopy$select = TRUE
		fcopy$weights = 0
		
		if (length(fname) == 0) {
			
			#grand mean
			fcopy$weights = 1 / nrow(fcopy)
			
		} else {
			
			for (i in 1:length(fname)) {
				fcopy$select[ fcopy[, fname[i] ] != lname[i] ] = FALSE
			}
			
			fcopy$weights[ fcopy$select ] = 1 / sum(fcopy$select)
			
			#get lower layer weights and subtract them
			
			lowerLayers = combinationLayers[ combinationLayers$layer < length(fname), ]
			for (i in 1:nrow(lowerLayers)) {
				
				fnameNext = combinations[[ lowerLayers[i,"combination"] ]]
				lnameNext = lname[fname %in% fnameNext] #handles fnameNext == character(0) appropriately.
				
				nextWeights = getWeights_internal(factors, fnameNext, lnameNext)
				
				fcopy$weights = fcopy$weights - nextWeights$weights
			}
			
		}
		
		fcopy
		
	}
	
	getWeights_internal(fcopy=factors, fname=fNames, lname=lNames)
}

#' Posterior Distributions of Main Effect and Interaction Parameters
#' 
#' @param results A results object.
#' @param param The name of a parameter with condition effect.
#' @param fNames A character vector giving the names of factors to use. The interaction of the factors will be the effect that is used. If there is only one factor, the main effect will be used.
#' 
#' @return A matrix with row being iterations and columns being effect parameters. The columns are named with the following scheme: "F1.L1:F2.L2" where "Fn" is the name of a factor and "Ln" is the level of that factor.
#' 
#' @export
getEffectParameterPosteriors = function(results, param, fNames) {
	
	factors = results$config$factors
	
	postDs = matrix(NA, nrow=results$config$iterations, ncol=nrow(factors))
	for (i in 1:nrow(factors)) {
		postDs[,i] = results$posteriors[[ paste(param, "_cond[", factors$cond[i], "]", sep="") ]]
	}
	
	getEffectParameters_general(cellMeans = postDs, factors = factors, fNames = fNames)
}

getEffectWeightsMatrix = function(factors, fNames, uniqueFL) {
	
	colNames = rep("", nrow(uniqueFL))
	weights = NULL
	
	for (i in 1:nrow(uniqueFL)) {
		
		thisLevels = as.character(uniqueFL[i,]) #coerce to char vector
		
		w = getFactorByLevelWeights(factors, fNames, thisLevels)
		
		weights = cbind(weights, w$weights)
		
		colNames[i] = paste0(paste0(fNames, ".", thisLevels), collapse=":")
	}
	
	colnames(weights) = colNames
	
	weights
	
}

getEffectParameters_general = function(cellMeans, factors, fNames, uniqueFL = NULL, stripRedundant = FALSE) {
	
	if (is.null(uniqueFL)) {
		uniqueFL = unique(subset(factors, select = fNames))
	}
	
	#uniqueFL = unique(factors[ , fNames ])
	#if (length(fNames) == 1) {
	#	temp = list()
	#	temp[[fNames]] = uniqueFL
	#	uniqueFL = as.data.frame(temp, stringsAsFactors=FALSE)
	#}
	
	weights = getEffectWeightsMatrix(factors, fNames, uniqueFL)

	effects = cellMeans %*% weights
	
	if (stripRedundant) {
		
		keep = rep(TRUE, nrow(uniqueFL))
		for (f in fNames) {
			uniqueLevels = unique(uniqueFL[,f])
			keepLevels = uniqueLevels[1:(length(uniqueLevels) - 1)]
			
			keep = keep & (uniqueFL[,f] %in% keepLevels)
		}
		
		#keepUFL = unique(subset(uniqueFL, subset = keep, select = fNames))
		keepUFL = unique(uniqueFL[ keep, fNames ])
		if (length(fNames) == 1) {
			temp = list()
			temp[[fNames]] = keepUFL
			keepUFL = as.data.frame(temp, stringsAsFactors=FALSE)
		}

		keepColNames = rep("", nrow(keepUFL))
		for (i in 1:nrow(keepUFL)) {
			fact = names(keepUFL)
			levels = as.character(keepUFL[i,]) #coerce to char vector
			keepColNames[i] = paste0(paste0(fact, ".", levels), collapse=":")
		}
		
		effects = effects[, keepColNames]
		
	}
	
	effects
}

testEffect_general = function(priorDs, postDs, factors, fNames, uniqueFL = NULL, 
											devianceFunction = NULL, testFunction = NULL) 
{

	priorEffects = getEffectParameters_general(priorDs, factors, fNames, uniqueFL=uniqueFL, stripRedundant=FALSE)
	postEffects = getEffectParameters_general(postDs, factors, fNames, uniqueFL=uniqueFL, stripRedundant=FALSE)
	
	if (is.null(devianceFunction)) {
		nc = ncol(priorEffects)
		if (nc == 2) {
			devianceFunction = devianceFunction_absDif
		} else if (nc >= 3 && nc <= 4) { #this is not based on much testing...
			devianceFunction = devianceFunction_absSum
		} else {
			devianceFunction = stats::var #yeah, var
		}
	}
	
	if (is.null(testFunction)) {
		testFunction = testFunction_savageDickey
	}
	
	priorDistance = apply(priorEffects, 1, devianceFunction)
	postDistance = apply(postEffects, 1, devianceFunction)
	
	testFunction(priorDistance, postDistance)
}


testEffect_specialized = function(results, factors, param, fNames, uniqueFL = NULL, 
																	priorSamples = 1e5, devianceFunction = NULL, testFunction = NULL) 
{
	
	#get priors and posteriors
	priorDs = sampleFromConditionEffectPriors(results, factors, param, priorSamples)
	
	postDs = matrix(NA, nrow=results$config$iterations, ncol=nrow(factors))
	for (i in 1:nrow(factors)) {
		postDs[,i] = results$posteriors[[ paste(param, "_cond[", factors$cond[i], "]", sep="") ]]
	}
	
	testEffect_general(priorDs, postDs, factors=factors, fNames=fNames, uniqueFL=uniqueFL,
										 devianceFunction=devianceFunction, testFunction=testFunction)
	
}

testMEI_single = function(results, param, priorSamples = 1e5, doPairwise = FALSE, devianceFunction = NULL, testFunction = "Savage-Dickey") {
	
	if (is.character(testFunction)) {
		if (testFunction == "Savage-Dickey") {
			testFunction = testFunction_savageDickey
		} else if (testFunction == "Interval") {
			testFunction = function(prior, posterior) {
				testFunction_interval(prior, posterior, p=0.01)
			}
		} else {
			stop("Invalid test function name.")
		}
	}
	if (!is.function(testFunction)) {
		stop("testFunction is not a function.")
	}
	
	factorsToTest = getFactorsForConditionEffect(results$config, param)
	if (length(factorsToTest) == 0) {
		return(NULL)
	}

	#This is basically C/P from another function. Maybe wrap it up?
	combinationLayers = NULL
	
	for (layer in 1:length(factorsToTest)) {
		r = utils::combn(factorsToTest, layer)
		r = t(r)
		for (j in 1:nrow(r)) {
			x = sort(r[j,])
			y = paste(x, collapse=":")
			
			combinationLayers = rbind(combinationLayers, data.frame(combination=y, layer=length(x), stringsAsFactors=FALSE))
		}
	}
	
	allTests = NULL
	for (i in 1:nrow(combinationLayers)) {
		
		theseFactors = strsplit(combinationLayers$combination[i], ":", fixed=TRUE)[[1]]
		
		factorName = paste0(theseFactors, collapse=":")
		
		thisRes = testEffect_specialized(results, results$config$factors, param, 
																		 fNames = theseFactors, priorSamples = priorSamples, 
																		 devianceFunction = devianceFunction,
																		 testFunction = testFunction)
		
		omnibus = data.frame(param=param, factor=factorName, levels="Omnibus", 
												 bf10=thisRes$bf10, bf01=thisRes$bf01, success=thisRes$success)
		allTests = rbind(allTests, omnibus)
		
		
		if (doPairwise && combinationLayers$layer[i] == 1) {
			#do pairwise comps
			
			thisFactor = theseFactors #just to be clear that there is only 1 factor
			
			allLevels = unique(results$config$factors[ , thisFactor ])
			comb = t(utils::combn(allLevels, 2))
			for (j in 1:nrow(comb)) {
				
				uniqueFL = list()
				uniqueFL[[thisFactor]] = comb[j,]
				uniqueFL = as.data.frame(uniqueFL, stringsAsFactors = FALSE)
				
				thisRes = testEffect_specialized(results, results$config$factors, param, 
																				 fNames = thisFactor, uniqueFL = uniqueFL, 
																				 priorSamples = priorSamples, 
																				 devianceFunction = devianceFunction_absDif, 
																				 testFunction = testFunction)
				
				levelNames = paste0(comb[j,], collapse=", ")
				
				pairwise = data.frame(param=param, factor=factorName, levels=levelNames, 
															bf10=thisRes$bf10, bf01=thisRes$bf01, success=thisRes$success)
				allTests = rbind(allTests, pairwise)
				
			} #j
		} #if doPairwise
	} #i
	
	allTests = allTests[ order(allTests$levels, allTests$factor), ]
	
	allTests
	
}



#' Test Main Effects and Interactions of Factors
#' 
#' This only supports one-factor designs or fully-crossed multi-factor designs. If your design is not fully crossed, you can use \code{\link{testConditionEffects}} to examine pairwise comparisons.
#' 
#' You must provide a \code{data.frame} containing the mapping from conditions to factor levels. This should be provided as in \code{results$config$factors}. See \code{\link{runParameterEstimation}} for more information. If you are using a one-factor design, this will have been created for you and you don't need to do anything.
#' 
#' This function uses kernel density estimation to estimate the densities of some relevant quantities. This procedure is somewhat noisy. As such, I recommend that you perform the procedure many times, the number of which can be configured with the \code{subsamples} argument. Then, aggregate results from the many repetitions of the procedure can be analyzed, which is done by default but can be changed by setting \code{summarize} to \code{FALSE}.
#' 
#' @param results A results object.
#' @param param Optional. Character vector of names of parameters to perform tests for. If NULL (default), is set to all parameters with condition effects.
#' @param summarize If TRUE (default), the results across subsamples will be summarized. If FALSE, the results from each of the subsamples will be returned. Those results can be later summarized with \code{\link{summarizeSubsampleResults}}.
#' @param subsamples
#' @param subsampleProportion
#' @param doPairwise Do pairwise tests of differences between levels of main effects (these are often called "post-hoc" tests).
#' @param devianceFunction A function used for calculating the deviation of the effect parameters.
#' 
#' @export
testMainEffectsAndInteractions = function(results, param=NULL, summarize=TRUE,
																					subsamples = 50, subsampleProportion = 1, 
																					doPairwise = FALSE, devianceFunction = NULL) 
{
	
	if (is.null(results$config$factors)) {
		stop('You must provide "results$config$factors".')
	}
	
	if (is.null(param)) {
		param = getParametersWithConditionEffects(results$config$conditionEffects)
	}
	
	subsampleIterationsToRemove = getSubsampleIterationsToRemove(results$config$iterations, subsamples, subsampleProportion)
	
	allBFs = NULL

	pb = utils::txtProgressBar(0, 1, 0, style=3)
	currentStep = 1
	lastStep = length(subsampleIterationsToRemove) * length(param)
	
	for (sub in 1:length(subsampleIterationsToRemove)) {
		
		if (length(subsampleIterationsToRemove[[sub]]) > 0) {
			resultSubsample = removeBurnIn(results, subsampleIterationsToRemove[[sub]])
		} else {
			resultSubsample = results
		}
		
		#Very important: The kernel density estimation procedure has a problem.
		#For the way in which it is used, the density depends on the number of
		#samples. More samples results in less density. Thus, the prior and posterior
		#sample counts must match.
		priorSamples = resultSubsample$config$iterations

		
		for (pInd in 1:length(param)) {
			
			result = testMEI_single(results=resultSubsample, param=param[pInd], 
															priorSamples=priorSamples, doPairwise=doPairwise, testFunction = "Savage-Dickey")
			allBFs = rbind(allBFs, result)
			
			utils::setTxtProgressBar(pb, value = currentStep / lastStep)
			currentStep = currentStep + 1
			
		}
	}
	
	close(pb)
	
	if (summarize) {
		rval = summarizeSubsampleResults(allBFs, aggregateBy = c("param", "factor", "levels"))
	
		rval$notOmnibus = rval$levels != "Omnibus"
		rval = rval[ order(rval$notOmnibus, rval$param, rval$factor, rval$levels), ]
		
		rval$notOmnibus = NULL
	} else {
		rval = allBFs
		attr(rval, "aggregateColumns") = c("param", "factor", "levels")
	}
	
	rval
	
}

#' Summarize Results from Multiple Subsamples
#' 
#' This function should be used with the values returned by \code{\link{testMainEffectsAndInteractions}} and \code{\link{testConditionEffects}} when the \code{summarize} argument is \code{FALSE}.
#' 
#' @param allBFs A data.frame containing the individual Bayes factors. It should have a format like the result of \code{\link{testMainEffectsAndInteractions}} or \code{\link{testConditionEffects}}. 
#' @param proportioniles Percentiles devided by 100 to calculate.
#' @param geometricZs A numberic vector of Z-values. BF quantiles will be calculated based on the geometric mean and standard deviation based on these Z-values.
#' @param consistencyCutoff A numeric vector of cutoffs. The proportion of Bayes factors above each cutoff is calculated. Defaults to 1. If you want a more stringent criterion in which consistent results have to be not just in the right direction, but also reasonably large, you could set it to 10, for example.
#' @param logBF Summarize log Bayes factors? If \code{FALSE}, no logs will be taken. If \code{TRUE}, log base 10 BFs will be used. If a numeric value, that value will be used as the base for the logarithm.
#' @param aggregateBy Columns of allBFs to aggregate by. Typically not required as it is read from an attribute of allBfs called aggregateColumns.
#' 
#' @return A data frame containing summarized test results. It has the following columns:
#' \tabular{ll}{
#' 	\code{...} \tab The columns that were aggregated by. See the \code{aggregateBy} argument. \cr
#' 	\code{bfType} \tab The type of Bayes factor in this row. "10" means that the Bayes factor is in favor of the alternative hypothesis that there is an effect (e.g. two conditions differed; there is a main effect). "01" means that the Bayes factor is in favor of the null hypothesis that there is no effect.\cr
#' 	\code{bf} \tab The arithmetic mean Bayes factor. \cr
#' 	\code{sd} \tab Standard deviation of the Bayes factors. \cr
#' 	\code{geo.mean} \tab The geometric mean Bayes factor. Bayes factors estimated with the approach used in this package tend to vary exponentially, which makes the geometric mean a possibly better measure than the arithmetic mean. Note that the geometric mean and the median tend to be in closer agreement than the arithmetic mean and the median. \cr
#' 	\code{geo.sd} \tab Geometric standard deviation of the Bayes factors. Multiply the geo.mean by the geo.sd to go up one standard deviation. In general, \code{geo.mean * geo.sd^z} will give you the geometric value corresponding to the given z score. \cr
#' 	\code{Min, Median, Max} \tab The minimum, median, and maximum of the Bayes factors. \cr
#' 	\code{n\%} \tab Other percentiles, as given in the \code{proportioniles} argument.
#' }
#' 
#' @export
summarizeSubsampleResults = function(allBFs, proportioniles = c(0, 0.025, 0.5, 0.975, 1), 
															geometricZs = NULL, consistencyCutoff = 1, logBF = FALSE, aggregateBy = NULL) {
	

	if (logBF == TRUE) {
		logBF = 10
		if (consistencyCutoff == 1) {
			consistencyCutoff = 0
		}
	}
	
	if (is.null(aggregateBy)) {
		aggregateBy = attr(allBFs, "aggregateColumns")
		if (is.null(aggregateBy)) {
			stop("No columns to aggregate by.")
		}
	}
	
	uniqueBF = unique(allBFs[ , aggregateBy ])
	
	rval = NULL
	for (i in 1:nrow(uniqueBF)) {

		for (bf in c("bf01", "bf10")) {
			
			theseAgg = list()
			useRows = rep(TRUE, nrow(allBFs))
			for (agg in aggregateBy) {
				theseAgg[[agg]] = uniqueBF[ i, agg ]
				
				useRows = useRows & (allBFs[ , agg ] == uniqueBF[ i, agg ])
			}

			dfl = theseAgg
			dfl[["bfType"]] = substr(bf, 3, 4)
			
			x = allBFs[ useRows, bf ]
			if (logBF != FALSE) {
				x = log(x, base=logBF)
			}
			
			if (length(x) > 1) {

				dfl[["bf"]] = mean(x)
				dfl[["sd"]] = stats::sd(x)
				

				dfl[["geo.mean"]] = geoMean(x)
				dfl[["geo.sd"]] = geoSD(x)
				if (logBF != FALSE) {
					dfl[["geo.mean"]] = NA
					dfl[["geo.sd"]] = NA
				}
				if (length(geometricZs) > 0) {
					for (j in 1:length(geometricZs)) {
						gq = geoQ(geometricZs[j], dfl[["geo.mean"]], dfl[["geo.sd"]])
						n = paste0("geo.mean ", ifelse( geometricZs[j] > 0, "+ ", " "), geometricZs[j], " sd")
						dfl[[ n ]] = gq
					}
				}
				
				for (j in 1:length(consistencyCutoff)) {
					name = paste0("p > ", consistencyCutoff[j])
					dfl[[name]] = mean(x > consistencyCutoff[j])
				}
				

				#percentiles
				if (length(proportioniles) > 0) {
					qs = stats::quantile(x, proportioniles)
					if (proportioniles[1] == 0) {
						names(qs) = c("Min", names(qs)[-1])
					}
					if (any(proportioniles == 0.5)) {
						ind = which(proportioniles == 0.5)
						n = names(qs)
						n[ind] = "Median"
						names(qs) = n
					}
					if (proportioniles[length(proportioniles)] == 1) {
						names(qs) = c(names(qs)[-length(proportioniles)], "Max")
					}
					for (j in 1:length(qs)) {
						dfl[[ names(qs)[j] ]] = as.numeric(qs[j])
					}
				}
				
			} else {
				
				dfl[["bf"]] = x

			}
			
			if (length(names(rval)) > length(names(dfl))) {
				for (n in names(rval)) {
					if (!(n %in% names(dfl))) {
						dfl[[n]] = NA
					}
				}
			}
			
			dfl = as.data.frame(dfl, stringsAsFactors=FALSE, check.names=FALSE)
			
			rval = rbind(rval, dfl)
		}
		
	}
	
	rval
}


prettyPrintBFResults = function(bfRes, aggregateBy, aggregateByLevels, bfType = "10", 
																quantiles = c(0, 0.025, 0.5, 0.975, 1), geometricZs = c(-2, 2)) 
{

	
	keep = rep(TRUE, nrow(bfRes))
	for (i in 1:length(aggregateBy)) {
		keep = keep & bfRes[ , aggregateBy[i] ] == aggregateByLevels[i]
	}
	bfRes = bfRes[ keep, ]
	
	bfs = bfRes[ , paste0("bf", bfType) ]
	logbfs = log(bfs, base=10)

	
	qs = stats::quantile(bfs, quantiles)
	logqs = stats::quantile(logbfs, quantiles)
	
	qs = as.matrix(qs)
	logqs = as.matrix(logqs)
	
	qs = cbind(qs, logqs)
	colnames(qs) = c("Linear", "Log")
	
	if (quantiles[1] == 0) {
		rownames(qs) = c("Min", rownames(qs)[-1])
	}
	if (quantiles[length(quantiles)] == 1) {
		rownames(qs) = c(rownames(qs)[-length(quantiles)], "Max")
	}
	if (any(quantiles == 0.5)) {
		ind = which(quantiles == 0.5)
		n = rownames(qs)
		n[ind] = "Median"
		rownames(qs) = n
	}
	

	cat("Bayes factors in favor of the hypothesis that there ")
	if (bfType == "01") {
		cat("*is not* ")
	} else {
		cat("*is* ")
	}
	cat("an effect for ")
	cat( paste(paste(aggregateBy, aggregateByLevels, sep=": "), collapse=", ") )
	cat(".\n")
	
	cat("\nQuantiles:\n")
	print(qs)
	cat("\n")
	
	
	ms = matrix(0, nrow=3, ncol=2)
	rownames(ms) = c("Linear", "Log Linear", "Geometric")
	colnames(ms) = c("Mean", "SD")
	
	ms["Linear", "Mean"] = mean(bfs)
	ms["Linear", "SD"] = stats::sd(bfs)
	
	ms["Log Linear", "Mean"] = mean(logbfs)
	ms["Log Linear", "SD"] = stats::sd(logbfs)
	
	ms["Geometric", "Mean"] = geoMean(bfs)
	ms["Geometric", "SD"] = geoSD(bfs)

	print(ms)
	
	pInFavor = mean(bfs > 1)
	cat("\nProportion of BFs in favor of hypothesis: ")
	cat(pInFavor)
	cat(".\n")
	
	
	title = paste0("Log BF for ", paste(paste(aggregateBy, aggregateByLevels, sep=": "), collapse=", "))
	graphics::hist(logbfs, xlab="Log BF", main=title)
	graphics::abline(v = 0, lty=2)
}



#all(x > 0)
geoMean = function(x) {
	#prod(x)^(1 / length(x)) 
	#which translates to
	exp( 1/length(x) * sum(log(x)) )
}


geoSD = function(x) {
	m = geoMean(x)
	
	exp(sqrt(sum(log(x / m)^2) / length(x)))
}

geoZ = function(x, mu, sigma) {
	log(x / mu, base=sigma)
}

geoQ = function(z, mu, sigma) {
	mu * sigma ^ z
}



