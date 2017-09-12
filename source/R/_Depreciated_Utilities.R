
######################
# The functions addFactorsToDataFrame and getFactorLevelsFromConds 
# only really make sense for WP designs


#' Add Factor Levels to Data Frame
#' 
#' @param df A \code{data.frame} with a column named \code{cond} with condition names like those in \code{factors$cond}. \code{df} should not have any columns with the names of the factors.
#' @param factors A \code{data.frame} formatted like \code{results$config$factors}.
#' 
#' @return A \code{data.frame} the same as \code{df} except with additional columns with factor levels.
#' 
#' @export
addFactorsToDataFrame = function(df, factors) {
	
	tf = getFactorLevelsFromConds(factors, df$cond)
	tf$cond = NULL
	
	cbind(df, tf)
}

#' Factor Levels from Condition Names
#' 
#' @param factors A \code{data.frame} formatted like \code{results$config$factors}.
#' @param conds A character vector of condition names like those in \code{factors$cond}.
#' 
#' @return A \code{data.frame} with one row for each element in \code{conds}. The return value will have columns for each factor plus a column named \code{cond}.
#' 
#' @export
getFactorLevelsFromConds = function(factors, conds) {
	res = NULL
	for (i in 1:length(conds)) {
		thisRows = which(factors$cond == conds[i])
		
		if (length(thisRows) > 0) {
			
			res = rbind(res, factors[ thisRows, ])
		} else {
			#make NA factor levels
			temp = factors[ 1, ]
			temp[ , names(temp) ] = NA
			temp$cond = conds[i]
			res = rbind(res, temp)
		}
		
	}
	
	res
}

############################################################

# I guess this function is really slow for some reason
matrixAggregation = function(formula, factors, mat, rowfun, colfun = NULL, colfirst = FALSE) {
	
	f = stats::formula(formula)
	t = stats::terms(f)
	
	fot = attr(t,"term.labels")[ attr(t,"order") == 1 ]
	ufact = subset(factors, select = fot)
	ufact = unique(ufact)
	
	aggMat = NULL
	for (i in 1:nrow(ufact)) {
		
		urow = subset(ufact, subset = i == 1:nrow(ufact))
		fr = matchFactorRows(factors, urow)
		
		temp = subset(mat, select = fr)
		
		if (colfirst) {
			
			if (!is.null(colfun)) {
				temp = apply(temp, 2, colfun)
			}
			
			if (!is.null(rowfun)) {
				if (is.vector(temp)) {
					temp = matrix(temp, nrow=1, ncol=length(temp))
				}
				temp = apply(temp, 1, rowfun)
			}
			
		} else {
			
			if (!is.null(rowfun)) {
				temp = apply(temp, 1, rowfun)
			}
			
			if (!is.null(colfun)) {
				if (is.vector(temp)) {
					temp = matrix(temp, nrow=length(temp), ncol=1)
				}
				temp = apply(temp, 2, colfun)
			}
		}
		
		aggMat = cbind(aggMat, temp)
	}
	
	colnames(aggMat) = NULL #or something based on ufact
	rownames(ufact) = NULL
	
	list(m=aggMat, f=ufact)
	
}


############################################################

participantPosteriorSummary.Old = function(results, params=NULL, doCatActive=TRUE, credLevel = 0.95, fun=NULL) {
	
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