

#' Plot the Category Weights Function
#' 
#' Makes a plot of the weights function, given in Equation 5 of the Appendix of Hardman, Vergauwe, and Ricker (2017).
#' 
#' @param catMu A vector of category means.
#' @param catSelectivity The categorical selectivity parameter.
#' @param dataType One of `"circular"` or `"linear"`.
#' @param colors A vector of colors for the different categories.
#' @param lty A vector of line types for the different categories.
#' @param lwd A vector of line width for the different categories.
#' @param study The study angles at which the category weights are plotted (i.e. a grid of x-values).
#' @param axes Boolean. If `TRUE`, axes are plotted. If `FALSE`, no axes are plotted.
#' 
#' @return Invisibly, the densities used to make the plot. Each column is related to one catMu.
#'  
#' @seealso [`categoryWeightsFunction`] to get the vector of probabilities for a single study angle.
#' @export
#' 
#' @examples
#' \dontrun{
#' plotWeightsFunction(c(30, 90, 120, 200, 210, 220, 300), 15)
#' }
plotWeightsFunction = function(catMu, catSelectivity, dataType = "circular",
                               colors=grDevices::rainbow(length(catMu)), lty=1:length(catMu),
                               lwd=rep(1, length(catMu)), study = NULL, axes=TRUE) 
{
  if (is.null(study)) {
    if (dataType == "circular") {
      study = seq(0, 360, length.out=100)
    } else {
      stop("If dataType == \"linear\", study must be provided.")
    }
  }
  
  dens = matrix(0, nrow=length(study), ncol=length(catMu))
  
  xlab = if (dataType == "circular") "Study Angle" else "Study Value"
  
  graphics::plot(range(study), c(0,1), type='n', axes=FALSE, 
                 xlab=xlab, ylab="Probability to Select Category")
  graphics::box()
  
  if (axes) {
    graphics::axis(2)
    if (dataType == "circular") {
      graphics::axis(1, at=seq(0, 360, 45))
    } else if (dataType == "linear") {
      graphics::axis(1)
    }
  }
  
  for (i in 1:nrow(dens)) {
    dens[i,] = categoryWeightsFunction(study[i], catMu, catSelectivity, dataType=dataType)
  }
  for (i in 1:ncol(dens)) {
    graphics::lines(study, dens[,i], col=colors[i], lty=lty[i], lwd=lwd[i])
  }
  colnames(dens) = catMu
  rownames(dens) = study
  invisible(dens)
}


#' Calculate Probabilities of Assignment to Categories
#' 
#' This function is the weights function, given in Equation 5 of the Appendix of Hardman, Vergauwe, and Ricker (2017).
#' 
#' @param study A scalar study angle, in degrees.
#' @param catMu A vector of category means, in degrees.
#' @param catSelectivity The categorical selectivity parameter. For circular data, this is a standard deviation in degrees.
#' @param dataType One of `"circular"` or `"linear"`.
#' 
#' @return A vector of the probabilities that the study angle would be assigned to each of the categories centered on the catMu.
#' 
#' @seealso [`plotWeightsFunction`] to plot the weights function for all study angles.
#' 
#' @export
categoryWeightsFunction = function(study, catMu, catSelectivity, dataType = "circular") {
  
  catDensities = NULL
  if (dataType == "circular") {
    catDensities = dvmd(study, catMu, catSelectivity)
  } else if (dataType == "linear") {
    catDensities = stats::dnorm(study, catMu, catSelectivity)
  }
  
  #if the sum of the densities is tiny, give all categories equal weight
  if (sum(catDensities) < 1e-250) {
    catDensities = rep(1, length(catDensities))
  }
  
  catDensities / sum(catDensities)
}




#' Calculate Probabilities of Assignment to Categories
#' 
#' When `distribution` is `"default"`, this function is the weights function, given in Equation 5 of the Appendix of Hardman, Vergauwe, and Ricker (2017).
#' 
#' NOTE TO SELF: This is basically categoryWeightsFunction with a different name.
#' 
#' @param catMu A vector of category means, in degrees.
#' @param catSelectivity The categorical selectivity parameter. For circular data, this is a standard deviation in degrees.
#' @param psWidth The width of the platSpline distribution, in degrees.
#' @param dataType "circular" (default) or "linear".
#' @param distribution "default" (default) or "PlatSpline"
getCategorizationProb = function(catMu, catSelectivity, psWidth,
  dataType = "circular", 
  distribution = "default") 
{
  stop("getCategorizationProb is unimplemented.")
}

