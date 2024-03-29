
# ... is names of parameters, with a list for each giving breaks, range, and/or label for that parameter.
#
# createParameterSummaryPlotConfiguration(paramSymbols, pMem = list(label = "pMem"), catSD = list(breaks = seq(0, 50, 1)))
createParameterSummaryPlotConfiguration = function(paramSymbols, ...) {
	
	baseConfig = list()
	baseConfig$pMem = list(breaks=seq(0, 1, 0.1), range=c(0,1), label="Prob. in memory")
	
	baseConfig$pBetween = list(breaks=seq(0, 1, 0.1), range=c(0,1), label="Prob. Between-Item")
	baseConfig$pContBetween = list(breaks=seq(0, 1, 0.1), range=c(0,1), label="Prob. continuous WM")
	baseConfig$pContWithin = list(breaks=seq(0, 1, 0.1), range=c(0,1), label="Proportion cont. WM")
	
	baseConfig$pCatGuess = list(breaks=seq(0, 1, 0.1), range=c(0,1), label="Prob. categorical guess")
	
	baseConfig$catSD = list(breaks=10, label="Categorical Imprecision")
	baseConfig$catSelectivity = list(breaks=10, label="Categorical Selectivity")
	baseConfig$contSD = list(breaks=10, label="Continuous imprecision")
	
	baseConfig$catMu = list(label = "Category Location")
	baseConfig$catActive = list(breaks=10, label = "Number of Categories")
	
	
	for (n in names(baseConfig)) {
		
		lab = baseConfig[[ n ]]$label
		sym = paramSymbols[[ n ]]
		
		baseConfig[[n]]$label = bquote(.(lab)*" ("*.(sym)*")")
	}
	
	
	changes = list(...)
	for (parName in names(changes)) {
		for (n in names(changes[[parName]])) {
			baseConfig[[parName]][[n]] = changes[[parName]][[n]]
		}
	}
	
	baseConfig
}

getSingleParameterPlotConfig = function(res, parName) {
	symbols = getParameterSymbols(res$config$modelVariant)
	config = createParameterSummaryPlotConfiguration(symbols)
	config[[ parName ]]
}


getParameterSummaryPlotLayout = function(modelVariant) {
	
	if (modelVariant == "betweenAndWithin") {
		parameterOrder = c("pMem", "contSD", "pCatGuess", 
											 "pBetween", "pContBetween", "pContWithin", 
											 "catMu",
											 "catActive", "catSD", "catSelectivity")
		
		screenCells = getSSCellMatrix(4, 3)
		screenCells = mergeCells(screenCells, 7, 8)
		screenCells = mergeCells(screenCells, 7, 8) #yes, merge twice. This makes a 3 panel wide one for catMu
		
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
#' Each of the panels of the overall parameter summary plot can be made with different underlying plotting functions. 
#' These are [`plotParameterHistogram`], [`plotParameterLineChart`], and [`plotCatMu`].
#' 
#' @details
#' 
#' Plotting to the default plotting surface in R often does not work well due to the large number of plots that are made. 
#' To plot to a pdf instead, which tends to work better, you can use the `pdfFile` argument.
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

#' @param catMuPrec The width of each bin used to plot the `catMu` parameters.
#' @param factorOrder The order in which the factors are plotted. Only useful if there is more than one factor. The first factor is put on the x-axis of factorial line charts. The order of the others factors doesn't really matter.
#' @param cip The proportion of the posterior in the credible interval.
#' @param paramSymbols A named list like that returned by [`getParameterSymbols`] that gives plotting symbols for the parameters of the model.
#' @param pdfFile Name of a pdf file to plot to. If `"TEMP"`, a temporary file is used and opened if `openPdf == TRUE`.
#' @param pdfScale Multiplicative scale factor for the size of the pdf. A larger value makes the plotting area larger, which makes all of the contents appear to be smaller.
#' @param openPdf If `TRUE`, an attempt will be made to open the pdf file in a viewer.
#' 
#' @family generic functions
#' @family plotting functions
#'
#' @export
plotParameterSummary = function(res, catMuPrec = 2, factorOrder = NULL, cip = 0.95, paramSymbols = NULL, pdfFile = NULL, pdfScale = 1.5, openPdf=TRUE) 
{
	
	# This isn't used anywhere
	resultWasWP = FALSE
	
	if (resultIsType(res, "WP")) {
		
		# TODO: Why are you doing this conversion to BP???
		resultWasWP = TRUE
		
		groups = list()
		groups[[ defaultGroupName() ]] = res
		
		bpRes = combineGroupResults.BP(groups)
		
	} else if (resultIsType(res, "BP")) {
		bpRes = res
	} else if (resultIsType(res, "Parallel")) {
	  stop("plotParameterSummary does not support Parallel results objects.")
	}
	
	# Set up the pdf
	if (!is.null(pdfFile)) {
		
		if (pdfFile == "TEMP") {
			pdfFile = paste0(tempdir(), "/Parameter Summary.pdf")
		}
		
	  if (bpRes$config$modelVariant == "ZL") {
	    pdfSize = c(6, 3) * pdfScale
	  } else if (bpRes$config$modelVariant == "betweenItem") {
	    pdfSize = c(9, 9) * pdfScale
	  } else if (bpRes$config$modelVariant == "withinItem") {
	    pdfSize = c(9, 9) * pdfScale
	  } else if (bpRes$config$modelVariant == "betweenAndWithin") {
	    pdfSize = c(9, 12) * pdfScale
	  }
	  
	  grDevices::pdf(pdfFile, width = pdfSize[1], height = pdfSize[2])
	}

	
	availableFactorNames = getAllFactorNames(bpRes$config$factors, removeConstant = TRUE)
	
	if (is.null(paramSymbols)) {
		paramSymbols = getParameterSymbols(bpRes$config$modelVariant)
	}
	
	parameterConfiguration = createParameterSummaryPlotConfiguration(paramSymbols)
	
	factorOrder = valueIfNull(factorOrder, availableFactorNames)
	
	layout = getParameterSummaryPlotLayout(bpRes$config$modelVariant)
	
	graphics::close.screen(all.screens=TRUE) #in case it is already open
	graphics::split.screen(layout$screenCells)
	
	for (paramInd in 1:length(layout$parameterOrder)) {
		
		graphics::screen(paramInd)
		graphics::par(mar=c(4,5,2,1))
		
		parName = layout$parameterOrder[paramInd]
		
		pc = parameterConfiguration[[parName]]
		
		#Get the factors for this parameter that 1) are estimated and 2) are in the available factors in factors.
		parameterFactors = getFactorsForConditionEffect(bpRes, parName)
		parameterFactors = parameterFactors[ parameterFactors %in% availableFactorNames ]
		
		if (parName == "catMu") {
			
			plotCatMu(res, precision = catMuPrec)
			
		} else {
			if (length(parameterFactors) == 0) {
				
				plotParameterHistogram(bpRes, parName, xlab = pc$label, breaks = pc$breaks, xlim = pc$range)
				
			} else {
				if (any(!(parameterFactors %in% factorOrder))) {
					logWarning("factorOrder does not contain some or all factors.")
				} else {
					#reorder only if factorOrder contains all factors.
					parameterFactors = factorOrder[ factorOrder %in% parameterFactors ]
				}
				
				plotParameterLineChart(bpRes, parName, factorOrder = parameterFactors, ylab = pc$label, cip = cip)
			}
		}
		
		# Whatever was plotted, add panel letters
		graphics::mtext(paste(LETTERS[paramInd], ".", sep=""), side=3, line=0.2, adj=0, cex=1.2 * graphics::par()$cex)
	}
	
	#This doesn't destroy anything, it just lets any following plots go in a new, non-split surface
	graphics::close.screen(all.screens=TRUE)
	
	if (!is.null(pdfFile)) {
		grDevices::dev.off()
		
	  if (openPdf) {
  	  # Try to open the file
  		os = Sys.info()[['sysname']]
  		if (os == "Windows") {
  			system(paste0('cmd /c "', pdfFile, '"'), wait=FALSE)
  		} else if (os == "Darwin" || os == "Linux") {
  			system(paste0('open "', pdfFile, '"'), wait=FALSE)
  		}
    }
	}
}

#################
# Plot anything #
#################

#' Basic Plot of any Parameter
#' 
#' This function plots a single parameter using default behaviors. If you
#' need more control over how the parameters are plotted, see the underlying 
#' plotting functions which provide more control over plotting behavior: [`plotParameterHistogram`], [`plotParameterLineChart`], and [`plotCatMu`].
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param parName The name of a single parameter (e.g. `"pMem"`).
#' 
#' @family generic functions
#' @family plotting functions
#' 
#' @export
plotParameter = function(res, parName) {
	
	if (parName == "catMu") {
		plotCatMu(res)
	} else {

		parameterFactors = getFactorsForConditionEffect(res, parName)
			
		if (length(parameterFactors) == 0) {
			plotParameterHistogram(res, parName)
		} else {
			plotParameterLineChart(res, parName, factorOrder = parameterFactors)
		}
	}
	
}

#########################
# Histogram (0 factors) #
#########################

#' Plot Histogram of Participant Mean Parameter Values
#' 
#' Plots a histogram of the means of the participants' parameter value. The mean parameter value is indicated with a vertical line.
#' 
#' @param res A generic result object (see [`Glossary`]).
#' @param parName The name of a single parameter (e.g. `"pMem"`).
#' @param xlab Label to place on the x-axis. If `NULL`, defaults to a reasonable value.
#' @param breaks Passed as `breaks` argument of `hist`. Defaults to 10 for standard deviation parameters or `seq(0, 1, 0.1)` for probability parameters.
#' @param xlim Passed as `xlim` argument of `hist`.
#' 
#' @family generic functions
#' @family plotting functions
#' 
#' @export
plotParameterHistogram = function(res, parName, xlab = NULL, breaks = NULL, xlim = NULL) {
	
	if (parName == "catActive") {
		rval = plotParameterHistogram_catActive(res, xlab = xlab, breaks = breaks, xlim = xlim)
		return(invisible(rval))
	}
	
	pc = getSingleParameterPlotConfig(res, parName)
	
	xlab = valueIfNull(xlab, pc$label)
	breaks = valueIfNull(breaks, pc$breaks)
	xlim = valueIfNull(xlim, pc$range)

	ppce = getAllParameterPosteriors(res, parName, manifest = TRUE, format = "data.frame")
	
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
#' @param breaks Passed on to the `breaks` argument of `graphics::hist()`.
#'
#' @family plotting functions
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

# This helper function is used by the LineChart package
credIntErrBarFun_Base = function(x, alpha) {
	ps = c(alpha / 2, 1 - alpha / 2)
	qs = stats::quantile(x, ps)
	list(eb = qs, includesCenter = TRUE)
}


#' Plot Factor-Varying Parameter with a Line Chart
#' 
#' Use this function to plot line charts of posterior means and credible intervals 
#' for a single model parameter that varies with at least one factor. This function
#' can handle any number of factors, but the plots become messy.
#' 
#' @section catActive:
#' If `parName == "catActive"` and it is a BP design, this function plots the mean and credible interval for the number of active categories used by participants in different cells of the design.
#' 
#' @param res A generic results object (see [`Glossary`]).
#' @param parName The parameter to plot. Can be `catActive` for between-participants designs.
#' @param factorOrder The order in which the factors are plotted. Only useful if there is more than one factor. The first factor is put on the x-axis of factorial line charts. The order of the others factors doesn't really matter.
#' @param cip The proportion of the posterior in the credible interval.
#' @param xlab The label to put on the x-axis. Defaults to `factorOrder[1]`.
#' @param ylab The label to put on the y-axis. If `NULL`, defaults to a reasonable value.
#' @param legendPosition Where to put the legend in the plot. `NULL` means no legend and `"CHOOSE_BEST"` (default) means to try to place the legend so as to not overlap with the points or error bar ends.
#' @param plotSettings A `data.frame` of plotting settings such as made by `LineChart::buildGroupSettings()`. See that package for more information.
#' 
#' @return Invisibly, the plotting data frame used to create the plot.
#' 
#' @family generic functions
#' @family plotting functions
#' 
#' @export
plotParameterLineChart = function(res, parName, factorOrder = NULL, cip = 0.95,
																	xlab = NULL, ylab = NULL, 
																	legendPosition = "CHOOSE_BEST", plotSettings = NULL)
{
	
	# Is it possible to use participantPosteriorSummary instead of whatever this function is doing?
	
	if (parName == "catActive") {
		if (resultIsType(res, "BP")) {
			rval = plotParameterLineChart_catActive.BP(res, factorOrder = factorOrder, cip = cip)
			return(invisible(rval))
		} else {
			stop("For WP designs, catActive may only be plotted with a histogram as it does not vary with WP factors. Use plotParameterHistogram() instead.")
		}
	}
	
	if (is.null(factorOrder)) {
		factorOrder = getFactorsForConditionEffect(res, parName)
	}
	
	if (is.null(xlab)) {
		xlab = factorOrder[1]
	}
	
	if (is.null(ylab)) {
		ylab = getSingleParameterPlotConfig(res, parName)$label
	}
	
	factors = updateFactorsForConditionEffects(res, parName)
	factors = normalizeFactors(factors)
	
	condEff = getConditionEffects(res, parName, addMu = TRUE, manifest = TRUE)
	# Collapsing condition effects is done to deal with missing factors, averaging across factor levels
	cce = collapseConditionEffects(condEff, factors, usedFactors = factorOrder, aggFun = mean)
	condEff = cce$condEff
	
	allPost = reshapeMatrixToDF(cce$condEff$post, cce$uniqueFL)

	ebf = function(x) {
		credIntErrBarFun_Base(x, alpha = 1 - cip)
	}
	
	form = stats::as.formula( paste0("x ~ ", paste0(factorOrder, collapse = " * ")) )
	
	plotDf = LineChart::lineChart(form, allPost, errBarType = ebf, settings = plotSettings, legendPosition = legendPosition, xlab = xlab, ylab = ylab)
	graphics::axis(4, labels = FALSE)
	
	invisible(plotDf)
}


#########
# CatMu #
#########


#' Plot Posterior Densities of Category Means
#' 
#' Plot a histogram of the posterior densities of the category mean parameters. Due to the nature of category mean paramters being about to trade locations with one another, there is no distinction between the different parameters. Thus, this collapses across all of the individual parameters.
#' 
#' Only the active categories are plotted as the distribution of inactive categories is irrelevant.
#' 
#' @param res A generic results object (see [`Glossary`]). If the `colorGeneratingFunction` field has been added to the results object, it will be used to add color information to the plot.
#' @param precision The width of the bars used to calculate densities.
#' @param pnums A vector of the participant numbers for which to make the plot.
#' @param type The type of plot to make, either filled polygons (`"polygon"`) or lines (`"line"`). Defaults to `"polygon"` for WP designs and `"line"` for BP designs.
#' @param lwd Line width. Only used if `type == "line"`.
#' @param lty Line type. Only used if `type == "line"`.
#' @param groupColor The color of the lines for each group in BP designs. Should be in the same order as the names of the groups (i.e. the same order as `names(bpRes$groups)`). Ignored for WP designs.
#' @param legendPosition The position of the legend showing the line colors used for the different groups in BP designs. Ignored for WP designs.
#' 
#' @return Invisibly, the data used to make the plot.
#' 
#' @family plotting functions
#' @family generic functions
#' 
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
	
	post = convertPosteriorsToMatrices(results, parNames = c("catMu", "catActive"))
	
	catMu = post$catMu[ pnums, , ]
	catActive = post$catActive[ pnums, , ]
	
	plotData = catMu_getPlotData(catMu, catActive, precision = precision, 
	                             dataType = results$config$dataType, 
	                             responseRange = results$config$responseRange, 
	                             colorGeneratingFunction = results$colorGeneratingFunction)
	
	catMu_plotSetup(plotData, results$config$dataType)
	catMu_plotData(plotData$x, plotData$y, col = plotData$color, type = type, lwd = lwd, lty = lty)

	invisible(plotData)
}


plotCatMu.BP = function(bpRes, precision = 2, pnums = NULL, type = "line", groupColor = NULL, lwd = 1, lty = 1, legendPosition = "topright") 
{
	
	ng = length(names(bpRes$groups))
	
	groupColor = valueIfNull(groupColor, grDevices::rainbow(ng))
	
	if (length(groupColor) == 1) {
		groupColor = rep(groupColor, ng)
	}
	if (length(lwd) == 1) {
		lwd = rep(lwd, ng)
	}
	if (length(lty) == 1) {
		lty = rep(lty, ng)
	}
	
	pnums = valueIfNull(pnums, getAllPnums.BP(bpRes))
	pnl = pnumVectorToList(pnums)
	
	allPlotData = NULL
	groupPlotData = list()
	
	for (group in names(bpRes$groups)) {
		if (is.null(pnl[[group]])) {
			next
		}
		
		post = convertPosteriorsToMatrices(bpRes$groups[[group]], parNames = c("catMu", "catActive"))
		
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
		catMu_plotData(plotData$x, plotData$y, col = groupColor[i], type = type, lwd = lwd[i], lty = lty[i])
	}
	
	if (plotColorBar) {
		xpos = groupPlotData[[1]]$x
		plotColorWheelBar(xpos, xpos, colorBarYlim, colorGeneratingFunction = bpRes$colorGeneratingFunction, horiz=TRUE)
	}
	
	if (!is.null(legendPosition)) {
		graphics::legend(legendPosition, legend = names(bpRes$groups), col = groupColor, lty = lty, lwd = lwd)
	}
	
	invisible(allPlotData)
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
	
	if (length(activeCatMu) > 0) {
  	for (i in 1:(length(angles) - 1)) {
  		heights[i] = mean(activeCatMu >= angles[i] & activeCatMu < angles[i + 1])
  		colors[i] = cgf(angles[i] + precision/2)
  	}
	}
	
	# The last point is the first point. Copy the first point to the last point.
	heights[length(heights)] = heights[1]
	colors[length(colors)] = colors[1]
	
	data.frame(x = angles, y = heights, color = colors, stringsAsFactors = FALSE)
}

catMu_plotData = function(x, y, col = "black", type = "line", lwd = 1, lty = 1) {
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

#############
# catActive #
#############

plotParameterLineChart_catActive.BP = function(bpRes, factorOrder = NULL, cip = 0.95, ylab = NULL) {
	
	ica = getIterationCatActive(bpRes)
	
	name2type = getFactorTypeToName(bpRes$config$factors)
	factorOrder = valueIfNull(factorOrder, name2type$bp)
	factorOrder = factorOrder[ factorOrder %in% name2type$bp ]
	
	ebf = function(x) {
		credIntErrBarFun_Base(x, alpha = 1 - cip)
	}
	
	pc = getSingleParameterPlotConfig(bpRes, "catActive")
	
	ylab = valueIfNull(ylab, pc$label)
	
	formula = formula(paste0("x ~ ", paste(factorOrder, collapse=" * ")))
	plotDf = LineChart::lineChart(formula, ica, errBarType = ebf, ylab = ylab)
	graphics::axis(4, labels = FALSE)
	
	rval = list(raw = ica, plotted = plotDf)
	invisible(rval)
	
}

plotParameterHistogram_catActive = function(res, xlab = NULL, breaks = NULL, xlim = NULL) {
	
	pc = getSingleParameterPlotConfig(res, "catActive")
	
	xlab = valueIfNull(xlab, pc$label)
	breaks = valueIfNull(breaks, pc$breaks)
	xlim = valueIfNull(xlim, pc$xlim)
	
	df = getCatActiveDataFrame(res)
	
	agg = stats::aggregate(x ~ pnum, df, mean)
	
	plotHistWithMean(agg$x, xlab = xlab, xlim = xlim, breaks = breaks)
	
	invisible(agg)
}

