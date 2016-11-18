

getSSCellMatrix = function(nrow, ncol) {
	screenCells = matrix(0, nrow=nrow * ncol, ncol=4)
	for (i in 1:nrow) {
		for (j in 1:ncol) {
			
			xj = j
			yi = nrow - i + 1
			
			x = c((xj - 1)/ncol, xj/ncol)
			y = c((yi - 1)/nrow, yi/nrow)
			
			screenCells[ (i - 1)*ncol + j, ] = c(x[1], x[2], y[1], y[2])
			
		}
	}
	screenCells
}

mergeCells = function(cells, i1, i2) {
	dim1 = cells[i1,]
	dim2 = cells[i2,]
	
	mins = pmin(dim1, dim2)
	maxs = pmax(dim1, dim2)
	
	cells[ i1, ] = c(mins[1], maxs[2], mins[3], maxs[4])
	
	cells[ -i2, ]
}

#' Plot Histogram with Mean Indicated
#' 
#' The mean is indicated with a strong vertical line.
#' 
#' @param x Vector of values to plot.
#' @param xlab Label for the x-axis.
#' @param xlim Limits for the x-axis.
#' @param breaks Passed on to the \code{breaks} argument of \code{hist}.
#' 
#' @export
plotHistWithMean = function(x, xlab=NULL, xlim=NULL, breaks=10) {
	
	graphics::par(xpd=FALSE)
	
	if (is.null(xlim)) {
		d = max(x) - min(x)
		xlim = range(x) + 0.2 * c(-d, d)
	}
	
	graphics::hist(x, breaks=breaks, main="", xlab=xlab, xlim=xlim, 
								 ylab="Number of participants", col=grDevices::rgb(0.85, 0.85, 0.85))
	graphics::box()
	graphics::axis(4, labels=FALSE)
	graphics::abline(v=mean(x), lwd=2, lty=2)
}



#' Plot Posterior Means and Credible Intervals for a Single-Factor Design
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param The parameter for which to plot.
#' @param xlab The label to put on the x-axis.
#' @param ylab The label to put on the y-axis.
#' 
#' @export
plotPosteriorLineChart = function(results, param, xlab=NULL, ylab=param) {
	
	usedFactor = getFactorsForConditionEffect(results$config, param)
	if (length(usedFactor) > 1) {
		stop("The number of used factors is too high to use plotPosteriorLineChart. See plotFactorialPosteriorLineChart()")
	} else if (length(usedFactor) < 1) {
		#if no factors used
		stop("The number of used factors is too low to use plotPosteriorLineChart.")
	}
	
	if (is.null(xlab)) {
		xlab = usedFactor
	}
	
	usedLevels = unique(results$config$factors[ , usedFactor ])
	usedConds = rep(NA, length(usedLevels))
	usedCondLabels = rep(NA, length(usedLevels))
	for (i in 1:length(usedLevels)) {
		whichRows = which(results$config$factors[ , usedFactor ] == usedLevels[i])
		possibleCondsForThisLevel = results$config$factors[ whichRows, "cond" ]
		usedConds[i] = possibleCondsForThisLevel[1]
		
		if (length(results$config$factorNames) == 1) {
			usedCondLabels[i] = paste(possibleCondsForThisLevel, collapse=",")
		} else {
			usedCondLabels[i] = usedLevels[i]
		}
	}
	
	ci = posteriorMeansAndCredibleIntervals(results, param, credLevel = 0.95)
	
	ci = ci[ ci$cond %in% usedConds, ]
	
	maxH = max(ci$upper)
	minH = min(ci$lower)
	
	xPos = 1:nrow(ci)
	
	graphics::plot(xPos, ci$mean, pch=16, ylim=c(min(ci$lower), max(ci$upper)), axes=FALSE, xlab=xlab, ylab=ylab )
	graphics::box()
	graphics::axis(1, at=xPos, labels = usedCondLabels)
	graphics::axis(2, las=1)
	graphics::axis(4, labels=FALSE)
	
	graphics::lines(xPos, ci$mean)
	
	simpleErrorBars = function(x, yStart, yEnd, headWidth) {
		graphics::lines( rep(x, 2), c(yStart, yEnd))
		graphics::lines(c( x - headWidth / 2, x + headWidth / 2 ), rep(yEnd, 2))
	}
	
	for (i in 1:nrow(ci)) {
		simpleErrorBars(xPos[i], ci[i, "mean"], ci[i,"lower"], 0.1)
		simpleErrorBars(xPos[i], ci[i, "mean"], ci[i,"upper"], 0.1)
	}
	
}

#' Plot Posterior Means and Credible Intervals for a Factorial Design
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param param The parameter for which to plot.
#' @param factorOrder The order in which the factors are plotted. The first factor is put on the x-axis. The order of the others doesn't really matter.
#' @param xlab The label to put on the x-axis.
#' @param ylab The label to put on the y-axis.
#' @param legendPosition Where to put the legend in the plot. \code{NULL} means no legend.
#' @param plotSettings A data.frame of plotting settings such as made by LineChart::buildGroupSettings(). See that package for more information.
#' 
#' @note This function requres the LineChart package that is an optional part of the installation of this package. 
#'
#'  @export
plotFactorialPosteriorLineChart = function(results, param, factorOrder = NULL, xlab = NULL, ylab = param, legendPosition = "CHOOSE_BEST", plotSettings = NULL) {
	
	#require(LineChart)
	
	usedFactor = getFactorsForConditionEffect(results$config, param)
	if (length(usedFactor) < 1) {
		stop("No factors have condition effects.")
	}
	
	if (is.null(xlab)) {
		xlab = factorOrder[1]
	}
	
	factors = results$config$factors
	
	ci = posteriorMeansAndCredibleIntervals(results, param, credLevel = 0.95)
	ci$keep = TRUE
	for (i in 1:nrow(ci)) {
		thisFullName = paste(param, "_cond[", ci$cond[i], "]", sep="")
		source = getRootSourceConditionParameter(results, param, ci$cond[i])
		ci$keep[i] = (thisFullName == source)
	}
	
	if (is.null(factorOrder)) {
		factorOrder = getFactorsForConditionEffect(results$config, param)
	}
	
	#add factors to ci
	for (i in 1:nrow(factors)) {
		row = which(ci$cond == factors$cond[i])
		for (f in factorOrder) {
			ci[ row, f ] = factors[ i, f ]
		}
	}
	
	ci = ci[ ci$keep, ]
	
	form = stats::as.formula( paste0("mean ~ ", paste0(factorOrder, collapse = " * ")) )

	plotDf = LineChart::createPlottingDf(form, ci, settings=plotSettings)
	
	ci$group = "0"
	if (length(factorOrder) >= 2) {
		for (i in 1:nrow(ci)) {
			ci$group[i] = paste0(ci[ i, factorOrder[2:length(factorOrder)] ], collapse=":")
		}
	}
	
	for (i in 1:nrow(plotDf)) {
		
		ciRow = which(plotDf[i,"xLabels"] == ci[ , factorOrder[1] ] & plotDf[i,"group"] == ci$group)
		
		plotDf$errBarLower[i] = ci$lower[ciRow] - ci$mean[ciRow]
		plotDf$errBar[i] = ci$upper[ciRow] - ci$mean[ciRow]
	}
	LineChart::lineChartDf(plotDf, xlab = xlab, ylab = ylab)
	
	if (length(factorOrder) >= 2 && !is.null(legendPosition)) {
		LineChart::legendFromPlottingDf(legendPosition, plotDf)
	}
	
	invisible(plotDf)
}



#' Plot Posterior Densities of Category Means
#' 
#' Plot a histogram of the posterior densities of the category mean parameters. Due to the nature of category mean paramters being about to trade locations with one another, there is no distinction between the different parameters. Thus, this collapses across all of the individual parameters.
#' 
#' Only the active categories are plotted.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param precision The width of the bars used to calculate densities.
#' @param pnums A vector of the participant numbers for which to make the plot.
#' 
#' @export
plotCatMu = function(results, precision, pnums = NULL) {
	
	if (is.null(pnums)) {
		pnums = results$pnums
	}
	
	post = convertPosteriorsToMatrices(results, param=c("catMu", "catActive"))
	
	select = which(dimnames(post$catMu)[[1]] %in% pnums)
	
	post$catMu = post$catMu[ select, , ]
	post$catActive = post$catActive[ select, , ]
	
	actCm = as.vector(post$catMu[ post$catActive == 1 ])
	if (results$config$dataType == "circular") {
		actCm = actCm %% 360
	}
	
	if (results$config$dataType == "circular") {
		angles = seq(0, 360, precision)
	} else {
		angles = seq(results$config$responseRange[1], results$config$responseRange[2], precision)
	}
	
	heights = rep(0, length(angles) - 1)
	for (i in 1:(length(angles) - 1)) {
		heights[i] = mean(actCm >= angles[i] & actCm < (angles[i] + 1))
	}
	
	xStart = c(0, 360)
	if (results$config$dataType == "linear") {
		xStart = results$config$responseRange
	}
	
	graphics::plot(xStart, c(0, max(heights) * 1.05), type='n', 
								 xlab=bquote("Category Location ("*mu*")"), ylab="Posterior Density", axes=FALSE)
	graphics::box()
	if (results$config$dataType == "circular") {
		graphics::axis(1, at=seq(0, 360, 60))
	} else {
		graphics::axis(1)
	}
	graphics::axis(2)
	
	colorGeneratingFunction = function(a) { grDevices::rgb(0.5, 0.5, 0.5) }
	if (!is.null(results$colorGeneratingFunction)) {
		colorGeneratingFunction = results$colorGeneratingFunction
	}
	for (i in 1:length(heights)) {
		a = angles[i]
		
		col = colorGeneratingFunction(a + precision/2) #halfway between ends
		
		graphics::polygon( c(a, a, a+precision, a+precision), c(0, heights[i], heights[i], 0), col=col, border=FALSE)
	}
	
	angles = angles + precision/2
 	temp = cbind(angles[-length(angles)], heights)
 	colnames(temp) = c("angle", "density")
	invisible(as.data.frame(temp))
}


#' Plot Parameter Summaries
#' 
#' For each parameter with condition effects, a line chart is plotted with error bars giving a 95%
#' credible interval of the mean. For parameters without condition effects, a histogram of participant
#' parameters is plotted, with the mean indicated by a vertical dashed line.
#' 
#' In addition, a histogram of the number of active categories per participant is plotted, as is the
#' collapsed posterior distribution of the category locations. This is done by selecting all active 
#' category locations for all participants and plotting them. No means are taken, it's just the raw
#' posterior distribution of all categories collapsed across all participants.
#' 
#' If you want the category location part of the plot to be colored, you need to add a 
#' colorGeneratingFunction function to the results object (i.e. results$colorGeneratingFunction).
#' This is a function that takes an angle in the interval [0, 360) and produces a color value corresponding
#' to that angle.
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param paramSymbols A named list like that returned by \code{\link{getParameterSymbols}} that gives plotting symbols for the parameters of the model.
#' @param catMuPrec The width of each bin used to plot the catMu parameters, in degrees.
#'
#' @export
#'
plotParameterSummary = function(results, paramSymbols=NULL, catMuPrec=2) {

	
	if (is.null(paramSymbols)) {
		paramSymbols = getParameterSymbols(results$config$modelVariant)
	}

	
	graphics::close.screen(all.screens=TRUE) #in case it is already open


	parameterConfiguration = list()
	parameterConfiguration$pMem = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
																		 label=bquote("Prob. in memory ("*.(paramSymbols$pMem)*")"))
	
	parameterConfiguration$pBetween = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
																				 label=bquote("Prob. Between-Item ("*.(paramSymbols$pBetween)*")"))
	parameterConfiguration$pContBetween = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
																						 label=bquote("Prob. continuous WM ("*.(paramSymbols$pContBetween)*")"))
	parameterConfiguration$pContWithin = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
																						label=bquote("Prob. continuous WM ("*.(paramSymbols$pContWithin)*")"))
	
	parameterConfiguration$pCatGuess = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
																					label=bquote("Prob. of categorical guess ("*.(paramSymbols$pCatGuess)*")"))
	
	parameterConfiguration$catSD = list(breaks=10, label=bquote("Categorical Imprecision ("*.(paramSymbols$catSD)*")"))
	parameterConfiguration$catSelectivity = list(breaks=10, 
																							 label=bquote("Categorical Selectivity ("*.(paramSymbols$catSelectivity)*")"))
	parameterConfiguration$contSD = list(breaks=10, label=bquote("Continuous imprecision ("*.(paramSymbols$contSD)*")"))
	


	#do setup for betweenAndWithin model
	parameterOrder = c("pMem", "contSD", "pCatGuess", 
										 "pBetween", "pContBetween", "pContWithin", 
										 "catMu",
										 "catActive", "catSD", "catSelectivity")
	
	
	#setup the split screen
	screenCells = getSSCellMatrix(4, 3)
	screenCells = mergeCells(screenCells, 7, 8)
	screenCells = mergeCells(screenCells, 7, 8) #yes, twice. This makes a 3 panel wide one for catMu
	
	if (results$config$modelVariant == "betweenItem") {
		parameterOrder = c("pMem", "contSD", "pContBetween", 
											 "catActive", "catMu",
											 "catSD", "catSelectivity", "pCatGuess")
		
		screenCells = getSSCellMatrix(3, 3)
		screenCells = mergeCells(screenCells, 5, 6)
		
	} else if (results$config$modelVariant == "withinItem") {
		parameterOrder = c("pMem", "contSD", "pContWithin", 
											 "catActive", "catMu",
											 "catSD", "catSelectivity", "pCatGuess")
		
		screenCells = getSSCellMatrix(3, 3)
		screenCells = mergeCells(screenCells, 5, 6)
		
	} else if (results$config$modelVariant == "ZL") {
		parameterOrder = c("pMem", "contSD")
		
		screenCells = getSSCellMatrix(1, 2)
	}
	
	
	graphics::split.screen(screenCells)
	#par(mar=c(4,5,2,1))
	
	post = convertPosteriorsToMatrices(results)
	
	for (paramInd in 1:length(parameterOrder)) {

		graphics::screen(paramInd)
		graphics::par(mar=c(4,5,2,1))
		
		param = parameterOrder[paramInd]
		
		pc = parameterConfiguration[[param]]

		parameterFactors = getFactorsForConditionEffect(results$config, param)
		
		if (length(parameterFactors) == 0) {
			#no condition effects
			
			if (param == "catMu") {
				
				plotCatMu(results, precision = catMuPrec)
				
			} else if (param == "catActive") {
				
				meanCatActive = apply(post$catActive, 1, mean) * results$config$maxCategories
				plotHistWithMean(meanCatActive, xlab="Number of Categories", breaks=10)
				
			} else {
				#do histogram
				trans = getParameterTransformation(param, results)
				meanForParts = apply( trans(post[[param]]), 2, mean)
				plotHistWithMean(meanForParts, xlab=pc$label, breaks=pc$breaks, xlim=pc$range)
			}
			
		} else if (length(parameterFactors) == 1) {
			
			plotPosteriorLineChart(results, param, xlab=parameterFactors, ylab=pc$label)
			
		} else if (length(parameterFactors) >= 2) {
			
			plotFactorialPosteriorLineChart(results, param, factorOrder = parameterFactors, ylab=pc$label)
			
		}

		
		graphics::mtext(paste(LETTERS[paramInd], ".", sep=""), side=3, line=0.2, adj=0, cex=1.2 * graphics::par()$cex)
	}
	
	#This doesn't destroy anything, it just lets new plots go in a new, non-split surface
	graphics::close.screen(all.screens=TRUE)

}



# depreciated
if (FALSE) {
plotTwoFactorPosteriorLineChart = function(results, param, factorOrder = NULL, ylab = param, legendPosition = "CHOOSE_BEST", plotSettings = NULL) {
	require(LineChart)
	
	
	factors = results$config$factors
	
	ci = posteriorMeansAndCredibleIntervals(results, param, credLevel = 0.95)
	
	if (is.null(factorOrder)) {
		factorOrder = getFactorsForConditionEffect(results$config, param)
	}
	if (length(factorOrder) > 2) {
		factorOrder = factorOrder[1:2]
		warning("More than two factors provided. Only the first two will be used.")
	}
	
	#add factors to ci
	for (i in 1:nrow(factors)) {
		row = which(ci$cond == factors$cond[i])
		for (f in factorOrder) {
			ci[ row, f ] = factors[ i, f ]
		}
	}
	
	ci[, "F1"] = factors[ , factorOrder[1] ]
	ci[, "F2"] = factors[ , factorOrder[2] ]
	
	groupLevels = unique(factors[ , factorOrder[2]])
	if (is.null(plotSettings)) {
		plotSettings = LineChart::buildGroupSettings(groupLevels, suppressWarnings=TRUE)
	}
	
	plotDf = LineChart::createPlottingDf(mean ~ F1 * F2, ci, settings=plotSettings)
	for (i in 1:nrow(plotDf)) {
		ciRow = which(plotDf[i,"xLabels"] == ci$F1 & plotDf[i,"group"] == ci$F2)
		plotDf$errBarLower[i] = ci$lower[ciRow] - ci$mean[ciRow]
		plotDf$errBar[i] = ci$upper[ciRow] - ci$mean[ciRow]
	}
	LineChart::lineChartDf(plotDf, xlab = factorOrder[1], ylab = ylab)
	if (!is.null(legendPosition)) {
		LineChart::legendFromPlottingDf(legendPosition, plotDf)
	}
	
	invisible(plotDf)
}
}