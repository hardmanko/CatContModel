

#' Posterior Means and Credible Intervals for Participant Parameters
#' 
#' The credible intervals for catActive should be interpreted in the context of them only taking on integer values.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param params A vector of parameter names. If \code{NULL}, all valid parameters are used.
#' @param doCatActive If \code{TRUE}, the posterior mean and credible intervals of the \code{catActive} parameters will be calculated.
#' @param credLevel The credibility level of the credible intervals. Defaults to 0.95.
#' @param fun A user provided function that will be passed a vector of a single participant by condition posterior distribution and that should return a single value.
#' 
#' @return A data frame with several columns:
#' \tabular{ll}{
#'	\code{pnum} \tab The participant number.\cr
#' 	\code{cond} \tab The zero-indexed condition index.\cr
#' 	\code{param} \tab The parameter name.\cr
#' 	\code{mean} \tab The posterior mean.\cr
#' 	\code{ciLower,ciUpper} \tab The lower and upper bounds of the credible interval.\cr
#' 	\code{fun} \tab If \code{fun} was provided, the results of that function call.
#' }
#' 
#' @family WP functions
#' @export
participantPosteriorSummary = function(results, params=NULL, doCatActive=TRUE, credLevel = 0.95, fun=NULL) {
	
	resDf = NULL
	
	if (is.null(params)) {
		params = getAllParams(results, filter=TRUE)
	}
	
	cips = c( (1 - credLevel) / 2, (1 + credLevel) / 2 )
	
	for (param in params) {
		
		for (pnum in results$pnums) {
			
			for (cond in results$config$factors$cond) {
				
				partParam = getParameterPosterior(results, param, pnum=pnum, cond=cond, manifest=TRUE)
				
				qs = as.numeric( stats::quantile(partParam, cips) )
				
				temp = data.frame(pnum=pnum, cond=cond, param=param, 
													mean = mean(partParam), ciLower = qs[1], ciUpper = qs[2], stringsAsFactors = FALSE)
				
				if (!is.null(fun)) {
					temp$fun = fun(partParam)
				}
				
				resDf = rbind(resDf, temp)
			}
		}
	}
	
	if (doCatActive) {
		post = convertPosteriorsToMatrices(results, "catActive")
		
		for (pnum in results$pnums) {
			
			ca = post$catActive[ which(results$pnums == pnum), ,  ]
			
			itCount = apply(ca, 2, sum)
			
			qs = as.numeric( stats::quantile(itCount, cips) )
			
			temp = data.frame(pnum=pnum, cond="ALL_CONDS", param="catActive", 
												mean = mean(itCount), ciLower = qs[1], ciUpper = qs[2], stringsAsFactors = FALSE)
			
			if (!is.null(fun)) {
				temp$fun = fun(itCount)
			}
			
			resDf = rbind(resDf, temp)
			
		}
	}
	
	# Remove participant X condition cells with no data
	uniPC = unique(results$data[ , c("pnum", "cond") ])
	keep = rep(FALSE, nrow(resDf))
	for (pnum in unique(uniPC$pnum)) {
		conds = as.character(uniPC$cond[ uniPC$pnum == pnum ])
		conds = c(conds, "ALL_CONDS")
		keep = keep | (resDf$pnum == pnum & resDf$cond %in% conds)
	}
	resDf = resDf[ keep, ]
	
	resDf = resDf[ order(resDf$pnum, resDf$param, resDf$cond), ]
	
	resDf
}



# STOP! You cannot aggregate across participants to get means! You must do mu + condition effect
singleParamPMCI = function(res, param, aggregateBy = c("group", "cond", "pnum"), credLevel = 0.95, manifest = TRUE) {
	
	credIntFun = function(x) {
		qs = stats::quantile(x, c((1 - credLevel) / 2, (1 + credLevel) / 2))
		rval = c(base::mean(x), qs)
		names(rval) = c("mean", "lower", "upper")
		rval
	}
	
	if (resultIsType(res, "WP")) {
		aggregateBy = aggregateBy[ aggregateBy != "group" ]
	}
	
	df = getAllParameterPosteriors(res, param, manifest = manifest, format = "data.frame")
	
	
	lostAgg = aggregateBy[ !(aggregateBy %in% names(df)) ]
	if (length(lostAgg) > 0) {
		warning(paste0("Could not aggregate by the following variables because they were not available in the data: ", paste(lostAgg, collapse = ", ")))
		aggregateBy = aggregateBy[ aggregateBy %in% names(df) ]
	}

	
	form = formula( paste0("x ~ ", paste(aggregateBy, collapse = " * ")) )
	agg = aggregate(form, df, credIntFun)
	agg$param = param
	agg = cbind(agg[ , c("param", aggregateBy) ], agg$x)
	
	agg
}

multiParamPMCI = function(res, param = NULL, aggregateBy = c("group", "cond", "pnum"), credLevel = 0.95, manifest = TRUE) {
	
	if (is.null(param)) {
		param = getAllParams(res, filter=TRUE)
	}
	
	allCI = NULL
	for (pp in param) {
		thisCI = singleParamPMCI(res, pp, aggregateBy = aggregateBy, credLevel = credLevel, manifest = manifest)
		
		allCI = rbind(allCI, thisCI)
	}

	allCI
}

#partPMCI_1 = participantPosteriorSummary(res, c("pMem", "contSD"), doCatActive = FALSE)
#partPMCI_1 = partPMCI_1[ partPMCI_1$param != "catActive", ]

#partPMCI_2 = multiParamPMCI(res, c("pMem", "contSD"))

#pmci1 = posteriorMeansAndCredibleIntervals(res, c("pMem", "contSD"))


#participantPosteriorSummary = function(res, params=NULL, doCatActive=TRUE, credLevel = 0.95, fun=NULL) {
	
#}


