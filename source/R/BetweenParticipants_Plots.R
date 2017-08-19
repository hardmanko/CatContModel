


# ... is names of parameters, with a list for each giving breaks, range, and/or label for that parameter.
#
# createParameterSummaryPlotConfiguration(paramSymbols, pMem = list(label = "pMem"), catSD = list(breaks = seq(0, 50, 1)))
createParameterSummaryPlotConfiguration = function(paramSymbols, ...) {
	
	baseConfig = list()
	baseConfig$pMem = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
												 label=bquote("Prob. in memory ("*.(paramSymbols$pMem)*")"))
	
	baseConfig$pBetween = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
														 label=bquote("Prob. Between-Item ("*.(paramSymbols$pBetween)*")"))
	baseConfig$pContBetween = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
																 label=bquote("Prob. continuous WM ("*.(paramSymbols$pContBetween)*")"))
	baseConfig$pContWithin = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
																label=bquote("Proportion cont. WM ("*.(paramSymbols$pContWithin)*")"))
	
	baseConfig$pCatGuess = list(breaks=seq(0, 1, 0.1), range=c(0,1), 
															label=bquote("Prob. of categorical guess ("*.(paramSymbols$pCatGuess)*")"))
	
	baseConfig$catSD = list(breaks=10, label=bquote("Categorical Imprecision ("*.(paramSymbols$catSD)*")"))
	baseConfig$catSelectivity = list(breaks=10, 
																	 label=bquote("Categorical Selectivity ("*.(paramSymbols$catSelectivity)*")"))
	baseConfig$contSD = list(breaks=10, label=bquote("Continuous imprecision ("*.(paramSymbols$contSD)*")"))
	
	
	changes = list(...)
	for (param in names(changes)) {
		for (n in names(changes[[param]])) {
			baseConfig[[param]][[n]] = changes[[param]][[n]]
		}
	}
	
	baseConfig
}


getParameterSummaryPlotLayout = function(modelVariant) {
	
	if (modelVariant == "betweenAndWithin") {
		parameterOrder = c("pMem", "contSD", "pCatGuess", 
											 "pBetween", "pContBetween", "pContWithin", 
											 "catMu",
											 "catActive", "catSD", "catSelectivity")
		
		screenCells = getSSCellMatrix(4, 3)
		screenCells = mergeCells(screenCells, 7, 8)
		screenCells = mergeCells(screenCells, 7, 8) #yes, twice. This makes a 3 panel wide one for catMu
		
	} else if (modelVariant == "betweenItem") {
		parameterOrder = c("pMem", "contSD", "pContBetween", 
											 "catActive", "catMu",
											 "catSD", "catSelectivity", "pCatGuess")
		
		screenCells = getSSCellMatrix(3, 3)
		screenCells = mergeCells(screenCells, 5, 6)
		
	} else if (modelVariant == "withinItem") {
		parameterOrder = c("pMem", "contSD", "pContWithin", 
											 "catActive", "catMu",
											 "catSD", "catSelectivity", "pCatGuess")
		
		screenCells = getSSCellMatrix(3, 3)
		screenCells = mergeCells(screenCells, 5, 6)
		
	} else if (modelVariant == "ZL") {
		parameterOrder = c("pMem", "contSD")
		
		screenCells = getSSCellMatrix(1, 2)
	}
	
	list(parameterOrder = parameterOrder, screenCells = screenCells)
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




#' Plot Parameter Summaries
#' 
#' Parameters that vary with factors, a line chart is plotted with error bars giving a 95%
#' credible interval of the mean. For parameters that do not vary with factors, a histogram of participant
#' parameters is plotted, with the mean indicated by a vertical dashed line.
#' 
#' Plotting to the default plotting surface in R often does not work well due to the large amount of plots that are made. To plot to a pdf instead, which tends to work better, you can either 1) use the `pdf` function before calling `plotParameterSummary` or 2) use the `asPdf` argument to plot to a pdf rather than the default plotting surface.
#' 
#' @section `catMu` and `catActive`:
#' In addition, a histogram of the number of active categories per participant is plotted, as is the
#' collapsed posterior distribution of the category locations. This is done by selecting all active 
#' category locations for all participants and plotting them. No means are taken, it's just the raw
#' posterior distribution of all categories collapsed across all participants.
#' 
#' If you want the category location part of the plot to be colored, you need to add a 
#' `colorGeneratingFunction` function to the results object (i.e. `results$colorGeneratingFunction`).
#' This is a function that takes an angle in the interval [0, 360) and produces a color value corresponding
#' to that angle.
#' 
#' @section Between-Participants Designs:
#' The between-participants version of this function does some things slightly differently than the standard version of the function. The main differences are related to the `catMu` and `catActive` parameters which cannot vary as a function of within-participants factors. When multiple groups, however, are used, those parameters are able to vary across groups. Thus, parameter summary plots should distinguish between the groups for `catMu` and `catActive`.
#' 
#' For `catMu`, each group is plotted as its own line. For `catActive`, the posterior mean and credible interval for the mean number of active categories is plotted. To be precise, the number of active categories per participant is calculated for each iteration. Then, the mean number of active categories across all participants is calculated for each iteration, producing the posterior distribution of the mean number of active categories, which is used to calculate the posterior mean and credible intervals.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param paramSymbols A named list like that returned by [`getParameterSymbols`] that gives plotting symbols for the parameters of the model.
#' @param catMuPrec The width of each bin used to plot the `catMu` parameters.
#' @param factorOrder The order in which the factors are plotted. Only useful if there is more than one factor. The first factor is put on the x-axis of factorial line charts. The order of the others factors doesn't really matter.
#' @param asPdf If `TRUE`, rather than plotting to the current plotting device, the plot is made in a pdf file that is opened in the default pdf viewer for your system.
#' @param pdfScale Multiplicative scale factor for the size of the pdf. A larger value makes the plotting area larger, which makes all of the contents appear to be smaller.
#' @param pdfDir The directory in which to make the pdf file. Defaults to a temporary directory.
#'
#' 
#' @family generic functions
#' @family plotting functions
#' @md
#' @export
plotParameterSummary = function(bpRes, paramSymbols = NULL, catMuPrec = 2, factorOrder = NULL, parameterConfiguration = NULL, asPdf = FALSE, pdfScale = 1.5, pdfDir = tempdir()) {
	
	resultWasWP = FALSE
	wpResults = NULL
	
	# Note that this function can take single group results (so it could just be renamed to plotParameterSummary)
	if (resultIsType(bpRes, "WP")) {
		resultWasWP = TRUE
		wpResults = bpRes
		
		bpRes = mergeGroupResults.BP(list(default = bpRes))
		bpRes$config$factors$BP_Factor = NULL # remove default extra factor
	}
	
	
	if (asPdf) {
		pdfFile = paste0(pdfDir, "/parameterSummary.pdf")
		
		if (bpRes$config$modelVariant == "ZL") {
			pdfSize = c(6, 3) * pdfScale
		} else if (bpRes$config$modelVariant == "betweenItem") {
			pdfSize = c(9, 9) * pdfScale
		} else if (bpRes$config$modelVariant == "withinItem") {
			pdfSize = c(9, 9) * pdfScale
		} else if (bpRes$config$modelVariant == "betweenAndWithin") {
			pdfSize = c(9, 12) * pdfScale
		}
		
		pdf(pdfFile, width = pdfSize[1], height = pdfSize[2])
	}
	

	availableFactorNames = getAllFactorNames(bpRes$config$factors, removeConstant = TRUE)
	
	if (is.null(paramSymbols)) {
		paramSymbols = getParameterSymbols(bpRes$config$modelVariant)
	}
	
	if (is.null(factorOrder)) {
		factorOrder = availableFactorNames
	}
	
	if (is.null(parameterConfiguration)) {
		parameterConfiguration = createParameterSummaryPlotConfiguration(paramSymbols)
	}
	
	layout = CatContModel:::getParameterSummaryPlotLayout(bpRes$config$modelVariant)
	
	graphics::close.screen(all.screens=TRUE) #in case it is already open
	graphics::split.screen(layout$screenCells)
	
	for (paramInd in 1:length(layout$parameterOrder)) {
		
		graphics::screen(paramInd)
		graphics::par(mar=c(4,5,2,1))
		
		param = layout$parameterOrder[paramInd]
		
		pc = parameterConfiguration[[param]]
		
		#Get the factors for this parameter that 1) are estimated and 2) are in the available factors in factors.
		parameterFactors = getFactorsForConditionEffect(bpRes, param)
		parameterFactors = parameterFactors[ parameterFactors %in% availableFactorNames ]
		
		if (param == "catMu") {
			
			if (resultWasWP) {
				plotCatMu.WP(wpResults, precision = catMuPrec)
			} else {
				plotCatMu.BP(bpRes, precision = catMuPrec)
			}
			
		} else if (param == "catActive") {
			
			if (length(parameterFactors) == 0) {
				
				post = convertPosteriorsToMatrices(bpRes, param = "catActive")
				
				meanCatActive = apply(post$catActive, 1, mean) * bpRes$config$maxCategories
				plotHistWithMean(meanCatActive, xlab="Number of Categories", breaks=10)
				
			} else {
				factorialCatActivePlot.BP(bpRes, factorOrder = parameterFactors)
			}
			
		} else {
			if (length(parameterFactors) == 0) {
				
				plotHistogram(bpRes, param, xlab = pc$label, breaks = pc$breaks, xlim = pc$range)
				
			} else {
				if (any(!(parameterFactors %in% factorOrder))) {
					warning("factorOrder does not contain some or all factors.")
				} else {
					#reorder only if factorOrder contains all factors.
					parameterFactors = factorOrder[ factorOrder %in% parameterFactors ]
				}
				
				plotFactorialLineChart(bpRes, param, factorOrder = parameterFactors, ylab=pc$label)
			}
		}
		
		# Whatever was plotted, add panel letters
		graphics::mtext(paste(LETTERS[paramInd], ".", sep=""), side=3, line=0.2, adj=0, cex=1.2 * graphics::par()$cex)
	}
	
	#This doesn't destroy anything, it just lets any following plots go in a new, non-split surface
	graphics::close.screen(all.screens=TRUE)
	
	if (asPdf) {
		dev.off()
		
		os = Sys.info()[['sysname']]
		if (os == "Windows") {
			system(paste0('cmd /c "', pdfFile, '"'), wait=FALSE)
		} else if (os == "Darwin" || os == "Linux") {
			system(paste0('open "', pdfFile, '"'), wait=FALSE)
		}

	}
}



#########################
# Histogram (0 factors) #
#########################

#' Plot Histogram of Participant Mean Parameter Values
#' 
#' Works for both WP and BP designs.
#' 
#' @param res A results object produced by either [`runParameterEstimation`] or [`mergeGroupResults.BP`].
#' @param param The name of the parameter to use.
#' @param xlab Label to place on the x-axis. Defaults to the value obtained from [`getParameterSymbols`].
#' @param breaks Passed as `breaks` argument of `hist`. Defaults to 10 for standard deviation parameters or `seq(0, 1, 0.1)` for probability parameters.
#' @param xlim Passed as `xlim` argument of `hist`.
#' 
#' @md
#' @family generic functions
#' @family plotting functions
#' @export
plotHistogram = function(res, param, xlab = NULL, breaks = NULL, xlim = NULL) {
	
	if (is.null(xlab)) {
		paramSymbols = getParameterSymbols(res$config$modelVariant)
		xlab = paramSymbols[[param]]
	}
	
	if (is.null(breaks)) {
		if (param %in% getProbParams(res)) {
			breaks = seq(0, 1, 0.1)
		} else {
			breaks = 10
		}
	}
	
	ppce = getAllParameterPosteriors(res, param, manifest = TRUE, format = "data.frame")
	
	dfagg = stats::aggregate(x ~ pnum, ppce, mean) # mean of all conditions and iterations is equivalent
		# to taking the within iteration mean of conditions and then the mean of all iterations
	
	plotHistWithMean(dfagg$x, xlab = xlab, breaks = breaks, xlim = xlim)
	
	invisible(dfagg)
	
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
		xlim = range(x) + 0.1 * c(-d, d)
	}
	
	graphics::hist(x, breaks=breaks, main="", xlab=xlab, xlim=xlim, 
								 ylab="Number of participants", col=grDevices::rgb(0.85, 0.85, 0.85))
	graphics::box()
	graphics::axis(4, labels=FALSE)
	graphics::abline(v=mean(x), lwd=2, lty=2)
}



#########################
# Factorial line charts #
#########################

# This function is a helper specific to the LineChart package
credIntErrBarFun_Base = function(x, alpha) {
	ps = c(alpha / 2, 1 - alpha / 2)
	qs = quantile(x, ps)
	list(eb = qs, includesCenter = TRUE)
}


#' Plot Factor-Varying Parameter with a Line Chart
#' 
#' Use this function to plot line charts of posterior means and credible intervals 
#' for a single model parameter that varies with at least one factor. This function
#' can handle any number of factors, but the plots become messy.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param param The parameter to plot.
#' @param factorOrder The order in which the factors are plotted. Only useful if there is more than one factor. The first factor is put on the x-axis of factorial line charts. The order of the others factors doesn't really matter.
#' @param xlab The label to put on the x-axis. Defaults to `factorOrder[1]`.
#' @param ylab The label to put on the y-axis. Defaults to the value from [`getParameterSymbols`] for `param`.
#' @param legendPosition Where to put the legend in the plot. `NULL` means no legend and `"CHOOSE_BEST"` (default) means to try to place the legend so as to not overlap with the points or error bar ends.
#' @param plotSettings A `data.frame` of plotting settings such as made by `LineChart::buildGroupSettings()`. See that package for more information.
#' 
#' @return Invisibly, the plotting data frame used to create the plot.
#' 
#' @md
#' @family generic functions
#' @family plotting functions
#' @export
plotFactorialLineChart = function(res, param, factorOrder = NULL,
																	xlab = NULL, ylab = NULL, 
																	legendPosition = "CHOOSE_BEST", plotSettings = NULL)
{
	
	if (is.null(factorOrder)) {
		factorOrder = getFactorsForConditionEffect(res, param)
	}
	
	if (is.null(xlab)) {
		xlab = factorOrder[1]
	}
	
	if (is.null(ylab)) {
		paramSymbols = getParameterSymbols(res$config$modelVariant)
		ylab = paramSymbols[[param]]
	}
	
	factors = updateFactorsForConditionEffects(res, param)
	factors = normalizeFactors(factors)
	
	condEff = getConditionEffects(res, param, addMu = TRUE, manifest = TRUE)
	# Collapsing condition effects is done to deal with missing factors, averaging across factor levels
	cce = collapseConditionEffects(condEff, factors, usedFactors = factorOrder, aggFun = mean)
	condEff = cce$condEff
	
	allPost = CatContModel:::reshapeMatrixToDF(cce$condEff$post, cce$uniqueFL)

	ebf = function(x) {
		CatContModel:::credIntErrBarFun_Base(x, alpha=0.05)
	}
	
	form = stats::as.formula( paste0("x ~ ", paste0(factorOrder, collapse = " * ")) )
	
	plotDf = LineChart::lineChart(form, allPost, errBarType = ebf, settings = plotSettings, legendPosition = legendPosition, xlab = xlab, ylab = ylab)
	axis(4, labels = FALSE)
	
	invisible(plotDf)
	
	

	
	
	#rval = CatContModel:::plotFactorialLineChart_Matrix(condEff$post, factors, factorOrder = factorOrder,
	#																		 xlab = xlab, ylab = ylab,
	#																		 legendPosition = legendPosition, plotSettings = plotSettings)
	
	#invisible(rval)
	
	
	#if (resultIsType(res, "WP")) {
	#	plotFun = plotFactorialLineChart.WP
	#} else if (resultIsType(res, "BP")) {
	#	plotFun = plotFactorialLineChart.BP
	#}
	
	#rval = plotFun(res, param = param, factors = factors, factorOrder = factorOrder,
	#				xlab = xlab, ylab = ylab,
	#				legendPosition = legendPosition, plotSettings = plotSettings)
	
	#invisible(rval)
}



plotFactorialLineChart.WP = function(results, param, factors = results$config$factors, 
																		 factorOrder = NULL, xlab = NULL, ylab = "Parameter", 
																		 legendPosition = "CHOOSE_BEST", plotSettings = NULL)
{

	factors = normalizeFactors(factors)
	
	condEff = getConditionEffects.WP(results, param, addMu = TRUE, manifest=TRUE)
	
	rval = plotFactorialLineChart_Matrix(condEff$post, factors, factorOrder = factorOrder,
																xlab = xlab, ylab = ylab,
																legendPosition = legendPosition, plotSettings = plotSettings)
	
	invisible(rval)
}

plotFactorialLineChart.BP = function(bpRes, param, factors = bpRes$config$factors, 
																		 factorOrder = NULL, xlab = NULL, ylab = "Parameter", 
																		 legendPosition = "CHOOSE_BEST", plotSettings = NULL)
{
	
	condEff = getConditionEffects.BP(bpRes, param, addMu = TRUE, manifest=TRUE)
	
	rval = plotFactorialLineChart_Matrix(condEff$post, factors, factorOrder = factorOrder,
																xlab = xlab, ylab = ylab,
																legendPosition = legendPosition, plotSettings = plotSettings)
	
	invisible(rval)
}

#depreciated
plotFactorialLineChart_Matrix = function(postCondEff, factors, factorOrder = NULL, 
																				 xlab = NULL, ylab = "Parameter", 
																				 legendPosition = "CHOOSE_BEST", plotSettings = NULL) 
{
	
	if (is.null(factorOrder)) {
		factorOrder = getAllFactorNames(factors)
	}
	
	if (is.null(xlab)) {
		xlab = factorOrder[1]
	}
	
	foFactors = subset(factors, select=c(factorOrder, "key"))
	uniqueFL = unique(subset(factors, select=factorOrder))
	
	uflKeys = CatContModel:::getMatchingKeysForUniqueFL(foFactors, uniqueFL)
	
	allPost = NULL
	for (i in 1:nrow(uniqueFL)) {
		
		pp = CatContModel:::getMultiConditionPosterior_Matrix(postCondEff, uflKeys[[i]])
		
		temp = data.frame(x = pp)
		for (n in names(uniqueFL)) {
			temp[,n] = uniqueFL[i,n]
		}
		allPost = rbind(allPost, temp)
		
	}
	
	ebf = function(x) {
		credIntErrBarFun_Base(x, alpha=0.05)
	}
	
	form = stats::as.formula( paste0("x ~ ", paste0(factorOrder, collapse = " * ")) )
	
	plotDf = LineChart::lineChart(form, allPost, errBarType = ebf, settings = plotSettings, legendPosition = legendPosition, xlab = xlab, ylab = ylab)
	axis(4, labels = FALSE)
	
	invisible(plotDf)
}


#############
# catActive #
#############


#' @md
#' @family BP functions
#' @family plotting functions
#' @export
factorialCatActivePlot.BP = function(bpRes, factorOrder = NULL) {
	
	fNames = getFactorTypeToName(bpRes$config$factors)
	factorOrder = valueIfNull(factorOrder, fNames$bp)
	factorOrder = factorOrder[ factorOrder %in% fNames$bp ]
	
	df = getCatActivePPCE_DF.BP(bpRes, factorOrder = factorOrder)
	
	# Average across participants
	df$iteration = 1:bpRes$config$iterations
	form = paste0("x ~ iteration * ", paste(factorOrder, collapse=" * "))
	agg = aggregate(formula(form), df, mean)
	
	ebf = function(x) {
		credIntErrBarFun_Base(x, alpha=0.05)
	}
	
	formula = formula(paste0("x ~ ", paste(factorOrder, collapse=" * ")))
	plotDf = LineChart::lineChart(formula, agg, errBarType = ebf, ylab = "Number of Categories")
	axis(4, labels = FALSE)

	rval = list(raw = df, plotted = agg)
	invisible(rval)
}

getCatActivePPCE_DF.BP = function(bpRes, factorOrder = NULL) {
	if (is.null(factorOrder)) {
		fNames = getFactorTypeToName(bpRes$config$factors)
		factorOrder = fNames$bp
	}
	
	gfact = bpRes$config$factors
	gfact = gfact[ , c("group", factorOrder) ]
	gfact = unique(gfact)
	
	post = convertPosteriorsToMatrices(bpRes, param = "catActive")
	ca = post$catActive
	
	meanCatActive = apply(ca, c(1, 3), sum)
	meanCatActive = t(meanCatActive)
	
	group_part = dimnames(meanCatActive)[[2]]
	parts = strsplit(group_part, split = ":", fixed=TRUE)
	design = NULL
	for (i in 1:length(parts)) {
		design = rbind(design, data.frame(group = parts[[i]][1], pnum = parts[[i]][2]))
	}
	for (n in factorOrder) {
		design[ , n ] = CatContModel:::substituteValues(design$group, gfact$group, gfact[, n])
	}
	
	reshapeMatrixToDF(meanCatActive, design)
}


##################
# CatMu plotting #
##################

catMu_plotLines = function(x, y, col = "black", type = "line", lwd = 1, lty = 1) {
	if (length(col) == 1) {
		col = rep(col, length(x))
	}
	if (length(lwd) == 1) {
		lwd = rep(lwd, length(x))
	}
	if (length(lty) == 1) {
		lty = rep(lty, length(x))
	}
	
	for (i in 1:(length(x) - 1)) {
		xt = x[c(i, i+1)]
		yt = y[c(i, i+1)]
		
		if (type == "line") {
			graphics::lines(xt, yt, col = col[i], lwd = lwd[i], lty = lty[i])
		} else if (type == "polygon") {
			graphics::polygon( rep(xt, each = 2), c(0, yt, 0), col=col[i], border=FALSE)
		}
	}
}


catMu_plotSetup = function(plotData, dataType, ylim_1 = 0) {
	
	graphics::plot(plotData$x, plotData$y, type='n', ylim = c(ylim_1, max(plotData$y) * 1.05),
								 xlab=bquote("Category Location ("*mu*")"), ylab="Posterior Density", axes=FALSE)

	graphics::box()
	if (dataType == "circular") {
		graphics::axis(1, at=seq(0, 360, 45))
	} else {
		graphics::axis(1)
	}
	graphics::axis(2)
	graphics::axis(4, labels = FALSE)

}

catMu_getPlotData = function(catMu, catActive, precision, dataType, responseRange = NULL, colorGeneratingFunction=NULL)
{
	
	if (dataType == "linear" && is.null(responseRange)) {
		stop("If dataType == linear, responseRange must be provided.")
	}
	
	if (dataType == "circular") {
		catMu = catMu %% 360
	}
	
	if (is.null(colorGeneratingFunction)) {
		cgf = function(a) { grDevices::rgb(0.5, 0.5, 0.5) }
	} else {
		cgf = colorGeneratingFunction
	}
	
	if (dataType == "circular") {
		angles = seq(0, 360, precision)
	} else {
		angles = seq(responseRange[1], responseRange[2], precision)
	}
	
	activeCatMu = catMu[ catActive == 1 ]
	
	heights = rep(0, length(angles))
	colors = rep("", length(angles))
	for (i in 1:(length(angles) - 1)) {
		heights[i] = mean(activeCatMu >= angles[i] & activeCatMu < (angles[i] + 1))
		colors[i] = cgf(angles[i] + precision/2)
	}
	heights[length(heights)] = heights[1]
	colors[length(colors)] = colors[1]
	
	# Center the x-value on the middle of the bin
	#angles = angles + precision / 2
	
	data.frame(x = angles, y = heights, color = colors, stringsAsFactors = FALSE)
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
#' @param groupColor Ignored for WP designs.
#' @param legendPosition Ignored for WP designs.
#' 
#' @md
#' @family plotting functions
#' @family generic functions
#' @export
plotCatMu = function(res, precision = 2, pnums = NULL, type = NULL, 
															 lwd = 1, lty = 1, groupColor = NULL, legendPosition = "topright") 
{
	
	if (resultIsType(res, "WP")) {
		
		rval = plotCatMu.WP(res, precision = precision, pnums = valueIfNull(pnums, res$pnums), 
															type = valueIfNull(type, "polygon"), 
															lwd = lwd, lty = lty)
		
	} else if (resultIsType(res, "BP")) {

		rval = plotCatMu.BP(res, precision = precision, pnums = pnums, 
															type = valueIfNull(type, "line"), 
															groupColor = groupColor, lwd = lwd, lty = lty, 
															legendPosition = legendPosition)
	}
	
	invisible(rval)
}


plotCatMu.WP = function(results, precision = 2, pnums = results$pnums, type = "polygon", lwd = 1, lty = 1) {
	
	post = convertPosteriorsToMatrices(results, param = c("catMu", "catActive"))
	
	catMu = post$catMu[ pnums, , ]
	catActive = post$catActive[ pnums, , ]
	
	plotData = CatContModel:::catMu_getPlotData(catMu, catActive, precision = precision, dataType = results$config$dataType, responseRange = results$config$responseRange, colorGeneratingFunction = results$colorGeneratingFunction)
	
	CatContModel:::catMu_plotSetup(plotData, results$config$dataType)
	CatContModel:::catMu_plotLines(plotData$x, plotData$y, col = plotData$color, type = type, lwd = lwd, lty = lty)

	invisible(plotData)
}


plotCatMu.BP = function(bpRes, precision = 2, pnums = NULL, type = "line", groupColor = NULL, lwd = 1, lty = 1, legendPosition = "topright") 
{
	
	ng = length(names(bpRes$groups))
	
	if (is.null(groupColor)) {
		groupColor = rainbow(ng)
	}
	if (length(groupColor) == 1) {
		groupColor = rep(groupColor, ng)
	}
	if (length(lwd) == 1) {
		lwd = rep(lwd, ng)
	}
	if (length(lty) == 1) {
		lty = rep(lty, ng)
	}
	
	
	if (is.null(pnums)) {
		pnums = getAllPnums.BP(bpRes)
	}
	pnl = pnumVectorToList(pnums)
	
	allPlotData = NULL
	groupPlotData = list()
	
	for (group in names(bpRes$groups)) {
		if (is.null(pnl[[group]])) {
			next
		}
		
		post = convertPosteriorsToMatrices(bpRes$groups[[group]], param = c("catMu", "catActive"))
		
		catMu = post$catMu[ pnl[[group]], , ]
		catActive = post$catActive[ pnl[[group]], , ]
		
		plotData = catMu_getPlotData(catMu, catActive, precision = precision, dataType = bpRes$config$dataType, responseRange = bpRes$config$responseRange, colorGeneratingFunction = bpRes$colorGeneratingFunction)
		
		groupPlotData[[group]] = plotData
		plotData$group = group
		allPlotData = rbind(allPlotData, plotData)
	}
	
	plotColorBar = !is.null(bpRes$colorGeneratingFunction)
	if (plotColorBar) {
		colorBarYlim = c(-0.03, -0.11) * (max(allPlotData$y) - min(allPlotData$y))
	} else {
		colorBarYlim = c(0, 0)
	}
	
	catMu_plotSetup(allPlotData, bpRes$config$dataType, ylim_1 = colorBarYlim[2])
	
	for (i in 1:length(names(bpRes$groups))) {
		
		group = names(bpRes$groups)[i]
		
		if (is.null(pnl[[group]])) {
			next
		}
		
		plotData = groupPlotData[[ group ]]
		catMu_plotLines(plotData$x, plotData$y, col = groupColor[i], type = type, lwd = lwd[i], lty = lty[i])
	}
	
	if (plotColorBar) {
		xpos = groupPlotData[[1]]$x
		plotColorWheelBar(xpos, xpos, colorBarYlim, colorGeneratingFunction = bpRes$colorGeneratingFunction, horiz=TRUE)
	}
	
	if (!is.null(legendPosition)) {
		legend(legendPosition, legend = names(bpRes$groups), col = groupColor, lty = lty, lwd = lwd)
	}
	
	invisible(allPlotData)
}



