

#' Make a Scatterplot with Color Bars in the Margins
#' 
#' @param data A data frame with \code{study} and \code{response} columns, both in degrees.
#' @param colorGeneratingFunction A function that takes an angle in degrees and returns a color corresponding to that angle. If \code{NULL}, no color bars will be plotted.
#' @param alpha Transparency for plotted points.
#' @param xylim A 2-length vector of the xlim and ylim for plotting.
#' @param overlap When drawing the color bars, how much adjacent bars should overlap. This overlap helps with pdf rendering issues.
#' 
#' @export
scatterplotWithColorBars = function(data, colorGeneratingFunction = NULL, alpha=0.1, overlap=0.1, xylim=c(0, 360)) {

	graphics::plot(data$study, data$response, xlim=xylim, ylim=xylim, 
			 pch=16, col=grDevices::rgb(0,0,0,alpha), xlab="Study angle", ylab="Response angle", axes=FALSE)
	graphics::box()
	
	at = NULL
	if (all(xylim == c(0, 360))) {
		at = seq(0, 360, 60)
	}
	graphics::axis(1, at=at)
	graphics::axis(2, at=at)
	
	opar = graphics::par()
	graphics::par(xpd=TRUE)
	
	if (!is.null(colorGeneratingFunction)) {
		step = 2
		angles = seq(xylim[1], xylim[2], step)
		
		for (a in angles) {
			col = colorGeneratingFunction(a + step/2)

			x = c(a, a, a + step, a + step)
			if (a < max(angles)) {
				x[3:4] = x[3:4] + overlap
			}
			y = xylim[2] + c(15, 40, 40, 15) / 360 * (xylim[2] - xylim[1])
			
			graphics::polygon(x, y, col=col, border = FALSE)
			graphics::polygon(y, x, col=col, border = FALSE)
		}
	}
	
	graphics::axis(3, at=at, labels=FALSE)
	graphics::axis(4, at=at, labels=FALSE)
	
	graphics::par(xpd=opar$xpd)
}



#' Posterior Predictive Distribution Plots
#' 
#' This function samples from the posterior predictive distribution for a participant in all conditions and plots the sampled data alongside the participant's actual data.
#' 
#' @param results A results object. The \code{colorGeneratingFunction} element will be used if available.
#' @param pnums The participant number(s) of the participant for whom you want to predict data.
#' @param conditions If not NULL, a vector of the conditions of the experiment to plot. If NULL (the default), all conditions are plotted.
#' @param rowLabels A vector of labels of length equal to the number of conditions. The labels are put on each row of plots. If \code{NULL}, rowLabels are made from the \code{conditions} list in \code{results}.
#' @param xylim A 2-length vector of the xlim and ylim for plotting.
#' @param alpha Transparency for plotted points.
#' @param plotPnum Include the participant number in the plot headers.
#'
#' @return Invisibly, the sampled data in a data frame.
#'
#' @export
#'
posteriorPredictivePlot = function(results, pnums, conditions=NULL, rowLabels=NULL, xylim=c(0, 360), alpha=0.5, plotPnum=FALSE) {

	plotPnum = plotPnum && (length(pnums) == 1) #don't plot pnum if there is more than 1

	tempD = results$data[ results$data$pnum %in% pnums, ]
	
	if (is.null(conditions)) {
		conditions = as.character(sort(unique(tempD$cond)))
	}
	
	
	allSampled = NULL
	
	for (pnum in pnums) {
		for (cond in conditions) {

			thisCondData = tempD[ tempD$cond == cond & tempD$pnum == pnum, ]

			sampled = NULL
			
			for (i in 1:nrow(thisCondData)) {
				
				s = sample(1:results$config$iterations, 1)
				
				param = getTransformedParameters(results, pnum=pnum, cond=cond, iteration=s)
				
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
	
	if (is.null(rowLabels)) {
		rowLabels = paste(results$conditions$type, conditions, sep=" ")
	}
	
	graphics::par(mfrow=c(length(conditions), 2), mar=c(5,4,3,1))
	
	for (condIndex in 1:length(conditions)) {

		cond = conditions[condIndex]
		
		thisCondData = tempD[ tempD$cond == cond, ]
		
		labelBase = paste(ifelse(plotPnum, paste("Part. ", pnums, ", ", sep=""), ""), rowLabels[condIndex], sep="")
		
		scatterplotWithColorBars(thisCondData, 
														 colorGeneratingFunction = results$colorGeneratingFunction, 
														 alpha=alpha, xylim=xylim)
		graphics::mtext(paste("Data - ", labelBase, sep=""), side=3, line=1.5, cex=graphics::par()$cex * 1.3, adj=0)
		
		scatterplotWithColorBars(allSampled[ allSampled$cond == cond, ], 
														 colorGeneratingFunction = results$colorGeneratingFunction, 
														 alpha=alpha, xylim=xylim)
		graphics::mtext(paste("Model - ", labelBase, sep=""), side=3, line=1.5, cex=graphics::par()$cex * 1.3, adj=0)
		
	}
	
	invisible(allSampled)
}





