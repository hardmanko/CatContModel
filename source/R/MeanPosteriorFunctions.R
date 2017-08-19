

#' Population/Condition Posterior Means and Credible Intervals
#' 
#' Calculates posterior means and credible intervals for the population means in each condition for the
#' given parameters. For each condition, condition effects are added to population means, the result is 
#' transformed to the manifest space, and the mean and credible interval for the manifest value is calculated.
#' Note that this is different from adding condition effects to participant-level parameters, tranforming 
#' the result, calculating on each iteration the mean of the transformed participant parameters, and 
#' calculating the posterior mean and credible interval of the iteration means. 
#' Using iteration means rather than population means will generally result in less than the true
#' amount of variability, which is why population means are used. Note, however, that this is a little strange,
#' because in the model, condition effects are not added to population means, but participant means.
#' 
#' In the return value of this function, the lower and upper columns give the endpoints of the credible interval.
#' 
#' Note that parameters without condition effects do not have values specific to conditions, so the `cond` will be `NA` for those parameters.
#'
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param params A vector of parameter names. If `NULL`, the default, all valid parameters are used.
#' @param credLevel The credibility level of the credible intervals. Defaults to 0.95.
#' 
#' @return A data.frame containing the results.
#' 
#' @md
#' @export
posteriorMeansAndCredibleIntervals = function(results, params=NULL, credLevel=0.95) {
	
	if (is.null(params)) {
		params = getAllParams(results, filter=TRUE)
	}
	
	rval = NULL
	
	halfa = (1 - credLevel) / 2
	
	for (param in params) {
		
		trans = getParameterTransformation(results, param)
		
		mu = results$posteriors[[ paste(param, ".mu", sep="") ]]
		
		for (cond in results$config$factors$cond) {
			
			condEff = results$posteriors[[ paste(param, "_cond[", cond, "]", sep="") ]]
			
			combined = trans(mu + condEff)
			
			quants = as.vector(stats::quantile(combined, c(halfa, 1 - halfa)))
			
			temp = data.frame(param=param, cond=cond, mean=mean(combined), lower=quants[1], upper=quants[2], stringsAsFactors=FALSE)
			
			paramWithConditionEffects = getParametersWithConditionEffects(results$config$conditionEffects)
			if (!(param %in% paramWithConditionEffects)) {
				temp$cond = "ALL_CONDS"
			}
			
			#For parameters without condition effects, only include one condition 
			#(the cornerstone condition, but it doesn't matter)
			if (param %in% paramWithConditionEffects || cond == results$config$cornerstoneConditionName) {
				rval = rbind(rval, temp)
			}
			
		}
	}
	
	rval = rval[ order(rval$param, rval$cond), ]
	
	rval
}



#PMCI collapsing across a set of usedConds
getMultiConditionPMCI = function(results, param, usedConds) {
	
	condMat = matrix(NA, nrow=results$config$iterations, ncol=length(usedConds))
	paramMu = results$posteriors[[ paste0(param, ".mu") ]]
	for (j in 1:length(usedConds)) {
		condEff = results$posteriors[[ paste0(param, "_cond[", usedConds[j], "]") ]]
		trans = getParameterTransformation(results, param)
		condMat[,j] = trans(paramMu + condEff)
	}
	
	#Take the mean of the condition means as the average across those conditions.
	meanOfConds = apply(condMat, 1, mean)
	
	qs = as.numeric( stats::quantile(meanOfConds, c(0.025, 0.975)) )
	
	pmci = data.frame(param = param, cond = paste(usedConds, collapse=","),
										mean = mean(meanOfConds),
										lower = qs[1], upper = qs[2])
	pmci
}

getMultiConditionPosterior_Matrix = function(postCondEff, usedKeys) {
	subCondEff = subset(postCondEff, select=usedKeys)
	
	#Take the mean of the condition means as the average across those conditions.
	#Note that you get a different result if you average then transform vs transform then average
	apply(subCondEff, 1, mean)
	
}

# unused
getMultiConditionPMCI_Matrix = function(postCondEff, usedKeys) {
	
	meanOfConds = getMultiConditionPosterior_Matrix(postCondEff, usedKeys)
	
	qs = as.numeric( stats::quantile(meanOfConds, c(0.025, 0.975)) )
	
	pmci = data.frame(key = paste(usedKeys, collapse=","),
										mean = mean(meanOfConds),
										lower = qs[1], upper = qs[2])
	pmci
}



