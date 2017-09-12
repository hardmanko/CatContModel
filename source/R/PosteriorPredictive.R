

#' Make a Scatterplot with Color Bars in the Margins
#' 
#' @param data A data frame with \code{study} and \code{response} columns, both in degrees.
#' @param colorGeneratingFunction A function that takes an angle in degrees and returns a color corresponding to that angle. If \code{NULL}, no color bars will be plotted.
#' @param alpha Transparency for plotted points.
#' @param overlap When drawing the color bars, how much adjacent bars should overlap. This overlap helps with pdf rendering issues.
#' @param colorStepSize The distance from the start of one color rectangle in a color bar to the next color rectangle.
#' @param xlim A 2-length vector of the xlim for plotting.
#' @param ylim A 2-length vector of the ylim for plotting.
#' @param xat A vector of the x-values for which axis tick marks are provided.
#' @param yat A vector of the y-values for which axis tick marks are provided.
#' 
#' @family plotting functions
#' @export
scatterplotWithColorBars = function(data, colorGeneratingFunction = NULL, alpha=0.1, overlap=0.1, colorStepSize=2, xlim=c(0,360), ylim=c(0,360), xat = NULL, yat = NULL) {

	graphics::plot(data$study, data$response, xlim=xlim, ylim=ylim, 
			 pch=16, col=grDevices::rgb(0,0,0,alpha), xlab="Study angle", ylab="Response angle", axes=FALSE)
	graphics::box()
	
	if (is.null(xat) && all(xlim == c(0, 360))) {
		xat = seq(0, 360, 60)
	}
	if (is.null(yat) && all(ylim == c(0, 360))) {
		yat = seq(0, 360, 60)
	}
	graphics::axis(1, at=xat)
	graphics::axis(2, at=yat)
	
	opar = graphics::par()
	graphics::par(xpd=TRUE)
	
	if (!is.null(colorGeneratingFunction)) {

		xAngles = seq(xlim[1], xlim[2], colorStepSize)
		ypos = ylim[2] + c(15, 40) / 360 * (ylim[2] - ylim[1])
		plotColorWheelBar(xAngles, xAngles, ypos, colorGeneratingFunction = colorGeneratingFunction, 
											horiz=TRUE, overlap = overlap)
		
		
		yAngles = seq(ylim[1], ylim[2], colorStepSize)
		xpos = xlim[2] + c(15, 40) / 360 * (xlim[2] - xlim[1])
		
		plotColorWheelBar(yAngles, yAngles, xpos, colorGeneratingFunction = colorGeneratingFunction, 
											horiz=FALSE, overlap = overlap)
		
	}
	
	graphics::axis(3, at=xat, labels=FALSE)
	graphics::axis(4, at=yat, labels=FALSE)
	
	graphics::par(xpd=opar$xpd)
}


# @param angles A vector of angles which will be passed to `colorGeneratingFunction` one at a time. 
# @param xpos A vector of x positions to place the colors. Of length equal to the length of `angles`.
# @param ypos A length-2 vector of the bottom and top y positions for the color bar.
# @param colorGeneratingFunction
# @param horiz
# @param overlap
plotColorWheelBar = function(angles, xpos, ypos, colorGeneratingFunction, horiz = TRUE, overlap = 0.1) {
	
	for (i in 1:(length(xpos) - 1)) {
		col = colorGeneratingFunction((angles[i] + angles[i+1])/2)
		
		x = rep(c(xpos[i], xpos[i+1]), each=2)
		y = c(ypos, rev(ypos))
		
		if (i < length(xpos) - 1) {
			x[3:4] = x[3:4] + overlap
		}
		
		if (!horiz) {
			temp = x
			x = y
			y = temp
		}
		
		graphics::polygon(x, y, col=col, border = FALSE)
	}
}

#' Posterior Predictive Distribution Plots
#' 
#' This function samples from the posterior predictive distribution for a participant in all conditions and plots the sampled data alongside the participant's actual data.
#' 
#' @param results The results from the [`runParameterEstimation`] function. The `colorGeneratingFunction` element will be used if available.
#' @param pnums The participant number(s) of the participant for whom you want to predict data. By default, all participants are used.
#' @param conditions If not `NULL`, a vector of the conditions of the experiment to plot. If `NULL` (the default), all conditions are plotted.
#' @param rowLabels A vector of labels of length equal to the number of conditions. The labels are put on each row of plots. If `NULL`, `rowLabels` are made from the `conditions` list in `results`.
#' @param xlim A 2-length vector of the xlim for plotting.
#' @param ylim A 2-length vector of the ylim for plotting.
#' @param xat A vector of the x-values for which axis tick marks are provided.
#' @param yat A vector of the y-values for which axis tick marks are provided.
#' @param alpha Transparency for plotted points.
#' @param plotPnum Boolean. Include the participant number in the plot headers. Only used if there is only 1 pnum.
#'
#' @return Invisibly, the sampled data in a `data.frame`.
#'
#' @md
#' @family WP functions
#' @family plotting functions
#' @export
#'
posteriorPredictivePlot = function(results, pnums = NULL, conditions=NULL, rowLabels=NULL, xlim=NULL, ylim=NULL, xat=NULL, yat=NULL, alpha=0.5, plotPnum=FALSE) {
	
	if (is.null(pnums)) {
		pnums = results$pnums
	}
	
	plotPnum = plotPnum && (length(pnums) == 1) #don't plot pnum if there is more than 1

	pnumData = results$data[ results$data$pnum %in% pnums, ]
	
	if (is.null(xlim)) {
		if (results$config$dataType == "circular") {
			xlim = c(0, 360)
		} else {
			xlim = range(pnumData$study)
		}
	}
	if (is.null(ylim)) {
		if (results$config$dataType == "circular") {
			ylim = c(0, 360)
		} else {
			ylim = results$config$responseRange
		}
	}
	
	
	if (is.null(conditions)) {
		conditions = as.character(sort(unique(pnumData$cond)))
	}
	
	
	allSampled = NULL
	
	for (pnum in pnums) {
		for (cond in conditions) {

			thisCondData = pnumData[ pnumData$cond == cond & pnumData$pnum == pnum, ]

			sampled = NULL
			
			for (i in 1:nrow(thisCondData)) {
				
				s = sample(1:results$config$iterations, 1)
				
				param = getSingleIterationParameters(results, pnum=pnum, cond=cond, iteration=s)
				
				thisSample = sampleDataFromModel(thisCondData$study[i], param, 
																				 modelVariant = results$config$modelVariant, 
																				 dataType = results$config$dataType, 
																				 responseRange = results$config$responseRange)
				
				sampled = rbind(sampled, thisSample)
			}

			sampled$pnum = pnum
			sampled$cond = cond
			allSampled = rbind(allSampled, sampled)
		}
	}
	
	
	graphics::par(mfrow=c(length(conditions), 2), mar=c(5,4,3,1))
	
	for (condIndex in 1:length(conditions)) {
		
		cond = conditions[condIndex]
		
		thisCondData = pnumData[ pnumData$cond == cond, ]
		
		if (!is.null(rowLabels)) {
			rowLabel = rowLabels[condIndex]
		} else {
			rowLabel = ""
			thisFactorLevels = results$config$factors[ results$config$factors$cond == cond, ]
			
			factorNames = getAllFactorNames(results$config$factors)
			
			for (i in 1:length(factorNames)) {
				n = factorNames[i]
				level = thisFactorLevels[1,n]
				
				rowLabel = paste0(rowLabel, n, " ", level)
				
				if (i < length(factorNames)) {
					rowLabel = paste0(rowLabel, ", ")
				}
			}
		}
		labelEnd = paste(ifelse(plotPnum, paste("Part. ", pnums, ", ", sep=""), ""), rowLabel, sep="")

		
		scatterplotWithColorBars(thisCondData, 
														 colorGeneratingFunction = results$colorGeneratingFunction, 
														 alpha=alpha, xlim=xlim, ylim=ylim, xat=xat, yat=yat)
		graphics::mtext(paste("Data - ", labelEnd, sep=""), side=3, line=1.5, cex=graphics::par()$cex * 1.3, adj=0)
		
		scatterplotWithColorBars(allSampled[ allSampled$cond == cond, ], 
														 colorGeneratingFunction = results$colorGeneratingFunction, 
														 alpha=alpha, xlim=xlim, ylim=ylim, xat=xat, yat=yat)
		graphics::mtext(paste("Model - ", labelEnd, sep=""), side=3, line=1.5, cex=graphics::par()$cex * 1.3, adj=0)
		
	}
	
	graphics::par(mfrow=c(1,1))
	
	invisible(allSampled)
}





