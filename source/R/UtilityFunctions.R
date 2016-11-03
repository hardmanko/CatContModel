
substituteValues = function(x, xvals, subvals) {
	y = rep(subvals[1], length(x)) #to get the type right
	for (i in 1:length(xvals)) {
		y[x == xvals[i]] = subvals[i]
	}
	y
}

getLatexParameterSymbols = function() {
	latex = list()
	
	latex$pMem = "\\pMem"
	latex$pBet = "\\pBet"
	latex$pOB = "\\pOB"
	latex$pOW = "\\pOW"
	latex$pAG = "\\pAG"
	
	latex$contSD = "\\contSD"
	latex$catSD = "\\catSD"
	latex$catSelectivity = "\\catSel"
	
	latex
}

#' Get a list of parameter symbols for plotting
#'
#' @param modelVariant String giving the name of the model variant to use. Can be accessed from \code{results$config$modelVariant}.
#' 
#' @return The list of symbols.
#'
#' @export
#'
getParameterSymbols = function(modelVariant) {
	paramSymbols = list()
	
	paramSymbols$pMem = bquote("P"^"M")
	paramSymbols$pBetween = bquote("P"^"B")
	
	paramSymbols$pContBetween = bquote("P"^"OB")
	if (modelVariant == "betweenItem") {
		paramSymbols$pContBetween = bquote("P"^"O")
	}
	
	paramSymbols$pContWithin = bquote("P"^"OW")
	if (modelVariant == "withinItem") {
		paramSymbols$pContWithin = bquote("P"^"O")
	}
	
	paramSymbols$pCatGuess = bquote("P"^"AG")
	
	paramSymbols$contSD = bquote(sigma^"O")
	paramSymbols$catSD = bquote(sigma^"A")
	paramSymbols$catSelectivity = bquote(sigma^"S")
	
	paramSymbols
}



#' Get Vectors of Parameter Names
#' 
#' @param results The results from the \code{\link{runParameterEstimation}} function.
#' @param filter Filter out parameter not in the current model. If FALSE, \code{results} can be NULL.
#' 
#' @export
getProbParams = function(results, filter=FALSE) {
	pp = c("pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess")
	if (filter) {
		if (results$config$modelVariant == "betweenItem") {
			pp = c("pMem", "pContBetween", "pCatGuess")
		} else if (results$config$modelVariant == "withinItem") {
			pp = c("pMem", "pContWithin", "pCatGuess")
		} else if (results$config$modelVariant == "ZL") {
			pp = c("pMem")
		}
	}
	pp
}

#' @rdname getProbParams
#' @export
getSdParams = function(results, filter=FALSE) {
	pp = c("contSD", "catSelectivity", "catSD")
	if (filter) {
		if (results$config$modelVariant == "ZL") {
			pp = c("contSD")
		}
	}
	pp
}

#' Logit Transformations
#' 
#' The logit transformation is \code{log(p/(1 - p))}.
#' 
#' The inverse logit transformation is \code{exp(f)/(1 + exp(f)))}.
#' 
#' @param p A probability (value in the interval [0,1]) to transform to (-Inf, Inf).
#' @param f A value in the interval \code{(-Inf, Inf)} to transform to [0, 1].
#' 
#' @return The transformed value.
#' 
#' @seealso \link{logitInverse}
#' 
#' @export
logit = function(p) { 
	log(p / (1 - p))
}

#' @rdname logit
#' @export
logitInverse = function(f) { 
	exp(f) / (1 + exp(f))
}

#' Convert Between Degrees and Radians
#'
#' @param rad A numeric vector of angles in radians.
#' @param deg A numeric vector of angles in degrees.
#' 
#' @return The converted angles.
#' 
#' @export
d2r = function(deg) {
	deg * pi / 180
}


#' @rdname d2r
#' @export
r2d = function(rad) {
	rad * 180 / pi
}



#' Von Mises Distribution Functions
#' 
#' These functions are little wrappers that parameterize the Von Mises
#' distribution in terms of degrees. Instead of precision, the von Mises
#' distribution is parameterized in terms of standard deviation (approximately),
#' which is defined as the square root of the inverse of precision.
#' 
#' @param n The number of realizations to draw.
#' @param x A quantile.
#' @param mu The center of the von Mises distribution.
#' @param sdDeg The standard deviation.
#' 
#' @export
rvmd = function(n, mu, sdDeg) {
	sampledRad = CircStats::rvm(n, d2r(mu), sdDeg_to_precRad(sdDeg))
	r2d(sampledRad)
}

#' @rdname rvmd
#' @export
dvmd = function(x, mu, sdDeg) {
	d = CircStats::dvm(d2r(x), d2r(mu), sdDeg_to_precRad(sdDeg))
	d * pi / 180 #scale to degree space
}


#' Convert Between Precision in Radians and Standard Deviation in Degrees
#' 
#' This model is parameterized in terms of standard deviation (square root of the inverse of precision) in degrees.
#' The density functions are all in terms of precision in radians. These functions go back and forth between them.
#' 
#' @param sdDeg A standard deviation in degrees.
#' @param precRad A precision in radians.
#' 
#' @export
#' @examples 
#' sdDeg_to_precRad(10)
#' precRad_to_sdDeg(0.05)
#' precRad_to_sdDeg(sdDeg_to_precRad(10)) == 10
sdDeg_to_precRad = function(sdDeg) {
	sdRad = d2r(sdDeg) #degrees to rad
	1/sdRad^2 #sd to precision
}

#' @rdname sdDeg_to_precRad
#' @export
precRad_to_sdDeg = function(precRad) {
	sdRad = 1 / sqrt(precRad) #precision to sd
	r2d(sdRad) #rad to deg
}

#' Calculate an Optionally Weighted Circular Mean
#' 
#' @param angles A vector of angles. If in radians, set \code{degrees} to FALSE.
#' @param weights A vector of weights.
#' @param degrees Whether the angles are in degrees or radians.
#' 
#' @return The circular mean of the angles.
#' 
#' @export
circMean = function(angles, weights=1, degrees=TRUE) {
	if (degrees) {
		angles = d2r(angles)
	}
	if (length(weights) <= 1) {
		weights = rep(weights, length(angles))
	}
	weights = weights / sum(weights)
	
	cosx = sum(cos(angles) * weights)
	sinx = sum(sin(angles) * weights)
	
	rval = atan2(sinx, cosx) %% (2 * pi)
	if (degrees) {
		rval = r2d(rval)
	}
	rval
}


#' Circular Absolute Distance Between Values
#' 
#' The cicular absolute distance is the smallest angular distance between two values. It is at most 180 degrees.
#' 
#' @param x, y Vectors of values in degrees or radians
#' @param degrees If TRUE, x and y are treated as though they are in degrees. If FALSE, x and y are treated as being in radians.
#' 
#' @export
circAbsDist = function(x, y, degrees=TRUE) {
	offset = 2 * pi
	if (degrees) {
		offset = 360
	}
	
	x = x %% offset
	y = y %% offset
	
	d1 = abs(x - y)
	d2 = abs((x + offset) - y) %% offset
	d3 = abs(x - (y + offset)) %% offset
	
	pmin(d1,d2,d3)
}


#' Get Parameter Transformation Function
#' 
#' Returns a function that will take parameter values in the latent space and convert them to the manifest space, or vice versa if \code{inverse} is TRUE. For probability parameters, the transformation is the inverse logit transformation. For standard deviation parameters, the transformation forces the parameter to be greater than some value, given by \code{results$config$minSD}.
#' 
#' @param param Name of a parameter, e.g. \code{pMem}.
#' @param results A results object.
#' @param inverse If TRUE, the inverse transformation is returned, if possible. Some transformations do not have an inverse.
#' 
#' @return A function of one argument that transforms the argument.
#' 
#' @export
getParameterTransformation = function(param, results, inverse=FALSE) {
	probParams = c("pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess", "pCatActive")
	sdParams = c("contSD", "catSelectivity", "catSD")
	
	if (param %in% probParams) {
		
		if (inverse) {
			transformation = function(x) { logit(x) }
		} else {
			transformation = function(x) { logitInverse(x) }
		}
		
		
	} else if (param %in% sdParams) {
		
		if (inverse) {
			transformation = function(x) { x } #Give warning if x == minSd? Say that you don't know the inverse?
		} else {
			minSd = results$config$minSD
			
			transformation = local( function(x) { pmax(x, minSd); })
		}
		
	} else if (param %in% c("catMu", "catActive")) {
		transformation = function(x) { x }
		
	}	else {
		warning( paste("Transformation not found for parameter: ", param, ".", sep="") )
		transformation = function(x) { x }
	}
	
	transformation
}

#' Calculate Probabilities of Assignment to Categories
#' 
#' This function is the weights function, given in Equation 5 of the Appendix of Hardman, Vergauwe, and Ricker (2016).
#' 
#' @param study A study angle, in degrees.
#' @param catMu A vector of category means, in degrees.
#' @param catSelectivity The categorical selectivity parameter, as standard deviation in degrees.
#' @param dataType One of \code{"circular"} or \code{"linear"}.
#' 
#' @return A vector of the probabilities that the study angle would be assigned to each of the categories centered on the catMu.
#' 
#' @seealso \code{\link{plotWeightsFunction}} to plot the weights function for all study angles.
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

#' Plot the Category Weights Function
#' 
#' Makes a plot of the weights function, given in Equation 5 of the Appendix of Hardman, Vergauwe, and Ricker (2016).
#' 
#' @param catMu A vector of category means, in degrees.
#' @param catSelectivity The categorical selectivity parameter, as standard deviation in degrees.
#' @param dataType One of \code{"circular"} or \code{"linear"}.
#' @param colors A vector of colors for the different categories.
#' @param lty A vector of line types for the different categories.
#' @param lwd A vector of line width for the different categories.
#' @param study The study angles at which the category weights are plotted (i.e. a grid of x-values).
#' @param axes If TRUE, axes are plotted.
#'  
#' @seealso \link{categoryWeightsFunction} to get the vector of probabilities for a single study angle.
#' @export
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
	
	graphics::plot(c(0, 360), c(0,1), type='n', axes=FALSE, 
			 xlab="Study Angle", ylab="Probability to Select Category")
	graphics::box()
	if (axes) {
		graphics::axis(2)
		graphics::axis(1, at=seq(0, 360, 60))
	}
	
	for (i in 1:nrow(dens)) {
		dens[i,] = categoryWeightsFunction(study[i], catMu, catSelectivity, dataType=dataType)
	}
	for (i in 1:ncol(dens)) {
		graphics::lines(study, dens[,i], col=colors[i], lty=lty[i], lwd=lwd)
	}
	dens
}







h_eq20 = function(mu_k, nu_k, mu_kp, nu_kp, catMuPriorSD, dataType) {
	
	if (dataType == "circular") {
		numDens = dvmd(mu_k - mu_kp, 0, catMuPriorSD)
		denDens = dvmd(0, 0, catMuPriorSD)
	} else if (dataType == "linear") {
		numDens = stats::dnorm(mu_k - mu_kp, 0, catMuPriorSD)
		denDens = stats::dnorm(0, 0, catMuPriorSD)
	}
	
	dens = 1 - nu_k * nu_kp * numDens / denDens
	dens^2
}

g_eq19 = function(mus, nus, k, catMuPriorSD, dataType) {
	mu_k = mus[k]
	nu_k = nus[k]
	mus = mus[-k]
	nus = nus[-k]
	
	d2 = h_eq20(mu_k, nu_k, mus, nus, catMuPriorSD, dataType=dataType)
	
	prod(d2)
}

SF_eq21 = function(mus, nus, k, catMuPriorSD, dataType = "circular", steps = 60, muRange = NULL) {
	
	if (dataType == "circular") {
		stepSize = 360 / steps
		
		xs = seq(0, 359.9, stepSize)
	} else if (dataType == "linear") {
		
		stepSize = (muRange[2] - muRange[1]) / steps
		
		xs = seq(muRange[1], muRange[2] - stepSize / 10, stepSize)
		
	}
	
	densSum = 0
	for (i in 1:length(xs)) {
		mus[k] = xs[i]
		densSum = densSum + g_eq19(mus, nus, k, catMuPriorSD, dataType=dataType)
	}
	
	1 / (densSum * stepSize)
}

f_eq18 = function(mus, nus, k, catMuPriorSD, dataType, muRange = NULL) {
	
	sf = SF_eq21(mus, nus, k, catMuPriorSD, dataType=dataType, muRange=muRange)
	
	gv = g_eq19(mus, nus, k, catMuPriorSD, dataType=dataType)
	
	gv * sf
	
}

#' Plot the Prior on catMu
#' 
#' Make a plot of the prior on catMu (and catActive, sort of) given the prior standard deviation and a vector of currently active categories.
#' 
#' @param catMuPriorSD The prior standard deviation on the notched prior on catMu and catActive.
#' @param catMu A vector of category locations, in degrees.
#' @param dataType One of \code{"circular"} or \code{"linear"}.
#' @param catActive A vector of category active parameters. By default, all of the provided catMus are active.
#' @param muRange Required only if \code{dataType == "linear"}. A length 2 vector giving the lower and upper bounds of the catMu parameters. This should generally be the same as the range of the data.
#' 
#' @return Invisibly, a data frame containing the densities in the plot.
#' 
#' @export
plotCatMuPrior = function(catMuPriorSD, catMu, dataType = "circular", catActive = rep(1, length(catMu)), muRange=NULL) {
	
	force(catActive)
	
	usedMus = c(0, catMu)
	catActive = c(1, catActive)
	
	if (dataType == "circular") {
		xs = seq(0, 360, 1)
	} else if (dataType == "linear") {
		xs = seq(muRange[1], muRange[2], 1)
	}
	
	dens = xs
	for (i in 1:length(xs)) {
		usedMus[1] = xs[i]
		dens[i] = f_eq18(usedMus, catActive, k = 1, catMuPriorSD, dataType=dataType, muRange=muRange)
	}
	
	graphics::plot(xs, dens, type='l', xlab="Angle", ylab="Prior Density")
	graphics::points(catMu, rep(0, length(catMu)), pch=16)
	
	invisible(data.frame(catMu=xs, dens=dens))
}




#' Produce Colors from a Warped HSV Color Space
#' 
#' This function provides one way of producing colors from variants of the HSV color space by converting an input angle into an output angle using a linear transformation. See the equations related to producing colors in the Method section for Experiment 1.
#' 
#' @param angle An angle in degrees for which a color should be produced.
#' @param inPoints See description.
#' @param outPoints See description.
#' 
#' @return A color.
#' 
#' @export
#' 
#' @examples
#' #Make a curried function that uses warpedHSVColorGeneratingFunction
#' #to produce colors like those used in Experiment 1 of Hardman, Vergauwe, and Ricker (2016).
#' exp1_colorGeneratingFunction = function(angle) {
#' 	inPoints = c(0, 180, 270, 360)
#' 	outPoints = c(0, 90, 230, 360)
#' 	warpedHSVColorGeneratingFunction(angle, inPoints, outPoints)
#' }
#' exp1_colorGeneratingFunction(90) #Get the color corresponding to 90 degrees.
#' 
warpedHSVColorGeneratingFunction = function(angle, inPoints, outPoints) {
	angle = angle %% 360
	
	interval = which(angle >= inPoints[1:(length(inPoints) - 1)] & angle < inPoints[2:length(inPoints)])
	
	RIL = inPoints[interval]
	RIH = inPoints[interval + 1]
	ROL = outPoints[interval]
	ROH = outPoints[interval + 1]
	
	propIn = (angle - RIL) / (RIH - RIL);
	pointOut = propIn * (ROH - ROL) + ROL;
	
	grDevices::hsv(pointOut / 360, 1, 1)
}

#' Likelihood Function for the Model
#' 
#' For given parameter values, a set of data, and a model variant, this calculates the likelihood for each observation in the data set.
#' 
#' Note that if you don't want to provide any catMu parameters, do not set catMu to NULL. Rather, set it to a zero-length vector with \code{vector(length = 0)}. Do not provide catActive parameters. Rather, filter out catMus that are inactive.
#' 
#' If \code{dataType} is \code{"circular"}, you should provide a value for the minSD argument.
#' 
#' If \code{dataType} is \code{"linear"}, you should provide a value for \code{responseRange}
#' 
#' @param param A list of parameters just like that returned by \code{\link{getTransformedParameters}}.
#' @param data A data.frame of the data for which you want to calculate the likelihood. Only the \code{study} and \code{response} columns are required.
#' @param modelVariant One of betweenItem, withinItem, and ZL.
#' @param dataType One of \code{"circular"} or \code{"linear"}.
#' @param responseRange A length 2 vector giving the theoretical minimum and maximum values of a response. Should be provided if \code{dataType} is \code{"linear"}.
#' @param minSD The minimum standard deviation of the Von Mises or normal distributions.
#' 
#' @return The provided \code{data} data.frame with an additional column named \code{likelihood} giving the likelihood for that observation.
#' 
#' @export
likelihood = function(param, data, modelVariant, dataType = "circular", 
											responseRange=NULL, minSD=NULL) 
{
	allParam = c(getSdParams(NULL), getProbParams(NULL), "catMu")
	if (!all(allParam %in% names(param))) {
		stop("Not all parameters were provided.")
	}
	
	
	config = list(modelVariant = modelVariant, dataType = dataType)
	if (dataType == "circular") {
		if (is.null(minSD)) {
			minSD = 1
			warning(paste("minSD not provided. It was set to ", minSD, ".", sep=""))
		}
		config$minSD = minSD
		
	} else if (dataType == "linear") {
		
		if (is.null(responseRange)) {
			responseRange = range(data$response)
			warning(paste("responseRange not provided. It was set to (", paste(responseRange, collapse=", "), ").", sep=""))
		}

		config$responseRange = responseRange
		config$catMuRange = responseRange
	}
	
	rval = CCM_CPP_likelihoodWrapper(param=param, data=data, config=config)
	
	data$likelihood = rval$likelihoods
	
	data
	
}





