getMultilevelFactorCrossing = function(factorsToTest) {
	combinationLayers = NULL
	
	for (layer in 1:length(factorsToTest)) {
		r = utils::combn(factorsToTest, layer)
		for (j in 1:ncol(r)) {
			x = sort(r[,j])
			y = paste(x, collapse=":")
			
			combinationLayers = rbind(combinationLayers, data.frame(combination=y, layer=length(x), stringsAsFactors=FALSE))
		}
	}
	
	combinationLayers
}

eval.string = function(x, envir=parent.frame()) {
	eval(parse(text=x), envir = envir)
}

cleanAndSummarizeMEIResults = function(BFs, summarize, aggregateBy) {
	if (any(BFs$success == FALSE)) {
		
		form = stats::formula(paste0("success ~ ", paste(aggregateBy, collapse=" * ")))
		ff = stats::aggregate(form, BFs, function(x) { sum(!x) })
		ff$failures = ff$success
		ff$success = NULL
		
		wmsg = paste0("Bayes factor estimation failed for ", sum(BFs$success == FALSE), " subsamples. The cases with failures have been printed to the console below. Failures have been stripped from the results.")
		warning( wmsg, immediate. = TRUE )
		cat("\n\nFailures:\n")
		print( ff[ ff$failures > 0, ] )
		cat("\n\n")
		
		BFs = BFs[ BFs$success, ]
	}
	
	if (summarize) {
		rval = summarizeSubsampleResults(BFs, aggregateBy = aggregateBy)
		
	} else {
		rval = BFs
		attr(rval, "aggregateColumns") = aggregateBy
	}
	
	ocode = paste0("order(", paste(paste0("rval$", aggregateBy), collapse=", "), ")")
	eval.string(ocode)
	
	rval = rval[ eval.string(ocode), ]
	
	rval
}


#' Summarize Results from Multiple Subsamples
#' 
#' This function should be used with the values returned by \code{\link{testMainEffectsAndInteractions}} and \code{\link{testConditionEffects}} when the \code{summarize} argument is \code{FALSE}. It summarizes Bayes factors across many repeated estimates of those Bayes factors.
#' 
#' @param BFs A data.frame containing the individual Bayes factors. It should have a format like the result of \code{\link{testMainEffectsAndInteractions}} or \code{\link{testConditionEffects}}. 
#' @param proportioniles Percentiles divided by 100 to calculate.
#' @param geometricZs A numeric vector of Z-values. Geometric BF quantiles will be calculated based on the geometric mean and standard deviation for each of these provided Z-values.
#' @param consistencyCutoff A numeric vector of cutoffs. The proportion of Bayes factors above each cutoff is calculated. Defaults to \code{c(1, 3, 10)}.
#' @param logBF Summarize log Bayes factors? If \code{FALSE}, no logs will be taken. If \code{TRUE}, log base 10 BFs will be used. If a numeric value, that value will be used as the base for the logarithm.
#' @param aggregateBy Columns of \code{BFs} to aggregate by. Typically not required as it is read from an attribute of \code{BFs} called \code{aggregateColumns}.
#' 
#' @return A data frame containing summarized test results. It has the following columns:
#' \tabular{ll}{
#' 	\code{...} \tab The columns that were aggregated by. See the \code{aggregateBy} argument. \cr
#' 	\code{bfType} \tab The type of Bayes factor in this row. "10" means that the Bayes factor is in favor of the alternative hypothesis that there is an effect (e.g. two conditions differed; there is a main effect). "01" means that the Bayes factor is in favor of the null hypothesis that there is no effect.\cr
#' 	\code{bf} \tab The arithmetic mean Bayes factor. \cr
#' 	\code{sd} \tab Standard deviation of the Bayes factors. \cr
#' 	\code{geo.mean} \tab The geometric mean Bayes factor. Bayes factors estimated with the approach used in this package tend to vary exponentially, which makes the geometric mean a possibly better measure than the arithmetic mean. Note that the geometric mean and the median tend to be in closer agreement than the arithmetic mean and the median. \cr
#' 	\code{geo.sd} \tab Geometric standard deviation of the Bayes factors. Multiply the geo.mean by the geo.sd to go up one standard deviation. In general, \code{geo.mean * geo.sd^z} will give you the geometric value corresponding to the given z score. \cr
#' 	\code{geo.mean + n SD} \tab The geometric mean "plus" \code{n} standard deviations, where the \code{n} values are given by the \code{geometricZs} argument. \cr
#' 	`p(BF > n)` \tab The proportion of Bayes factors greater than `n`. See `consistencyCutoff` to set `n`. \cr
#' 	\code{Min, Median, Max} \tab The minimum, median, and maximum of the Bayes factors. \cr
#' 	\code{n\%} \tab Other percentiles, as given in the \code{proportioniles} argument.
#' }
#' 
#' @md
#' @export
summarizeSubsampleResults = function(BFs, proportioniles = c(0, 0.025, 0.5, 0.975, 1), 
																		 geometricZs = NULL, consistencyCutoff = c(1, 3, 10), logBF = FALSE, 
																		 aggregateBy = NULL) {
	
	
	if (logBF == TRUE) {
		logBF = 10
		if (consistencyCutoff == 1) {
			consistencyCutoff = 0
		}
	}
	
	if (is.null(aggregateBy)) {
		aggregateBy = attr(BFs, "aggregateColumns")
		if (is.null(aggregateBy)) {
			stop("No columns to aggregate by.")
		}
	}
	
	uniqueBF = unique(BFs[ , aggregateBy ])
	
	rval = NULL
	for (i in 1:nrow(uniqueBF)) {
		
		for (bf in c("bf01", "bf10")) {
			
			theseAgg = list()
			useRows = rep(TRUE, nrow(BFs))
			for (agg in aggregateBy) {
				theseAgg[[agg]] = uniqueBF[ i, agg ]
				
				useRows = useRows & (BFs[ , agg ] == uniqueBF[ i, agg ])
			}
			
			dfl = theseAgg
			dfl[["bfType"]] = substr(bf, 3, 4)
			
			x = BFs[ useRows, bf ]
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
					name = paste0("p(BF > ", consistencyCutoff[j], ")")
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


#all(x > 0)
geoMean = function(x) {
	if (any(x <= 0)) {
		warning("All x must be > 0 to calculate the geometric mean.")
		return(NA)
	}
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

#######################################################################