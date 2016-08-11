

histWithMean = function(x, xlab=NULL, xlim=NULL, breaks=10) {
	
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

plotPosteriorLineChart = function(results, param, ylab=param) {
	
	ci = posteriorMeansAndCredibleIntervals(results, param, credLevel = 0.95)
	
	maxH = max(ci$upper)
	minH = min(ci$lower)
	
	xPos = 1:nrow(ci)
	
	graphics::plot(xPos, ci$mean, pch=16, ylim=c(min(ci$lower), max(ci$upper)), axes=FALSE, ylab=ylab, xlab=results$conditions$type )
	graphics::box()
	graphics::axis(1, at=xPos, labels = ci$cond)
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

plotCatMu = function(results, catMuPrec) {
	
	post = convertPosteriorsToMatrices(results, param=c("catMu", "catActive"))
	
	actCm = as.vector(post$catMu[ post$catActive == 1 ])
	if (results$config$dataType == "circular") {
		actCm = actCm %% 360
	}
	
	if (results$config$dataType == "circular") {
		angles = seq(0, 360, catMuPrec)
	} else {
		angles = seq(results$config$responseRange[1], results$config$responseRange[2], catMuPrec)
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
		
		col = colorGeneratingFunction(a + catMuPrec/2) #halfway between ends
		
		graphics::polygon( c(a, a, a+catMuPrec, a+catMuPrec), c(0, heights[i], heights[i], 0), col=col, border=FALSE)
	}
	
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
#' @param results A results object.
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

		
		if (param %in% results$config$parametersWithConditionEffects) {

			plotPosteriorLineChart(results, param, pc$label)
			
		} else {
			
			if (param == "catMu") {
				
				plotCatMu(results, catMuPrec)
				
			} else if (param == "catActive") {
				
				meanCatActive = apply(post$catActive, 1, mean) * results$config$maxCategories
				histWithMean(meanCatActive, xlab="Number of Categories", breaks=10)
				
			} else {
				#do histogram
				trans = getParameterTransformation(param, results)
				meanForParts = apply( trans(post[[param]]), 2, mean)
				histWithMean(meanForParts, xlab=pc$label, breaks=pc$breaks, xlim=pc$range)
			}
			
		}
		
		graphics::mtext(paste(LETTERS[paramInd], ".", sep=""), side=3, line=0.2, adj=0, cex=1.2 * graphics::par()$cex)
	}
	
	#This doesn't destroy anything, it just lets new plots go in a new, non-split surface
	graphics::close.screen(all.screens=TRUE)

}


