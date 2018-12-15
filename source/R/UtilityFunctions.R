
#' Open CatContModel Package Introduction Manual
#' 
#' Opens the manual for the package. 
#' 
#' @export
CatContModelManual = function() {
	utils::vignette("Introduction", "CatContModel")
}


substituteValues = function(x, xvals, subvals) {
	y = rep(subvals[1], length(x)) #to get the type right
	for (i in 1:length(xvals)) {
		y[x == xvals[i]] = subvals[i]
	}
	y
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
	ps = list()
	
	ps$pMem = bquote("P"^"M")
	ps$pBetween = bquote("P"^"B")
	
	ps$pContBetween = bquote("P"^"OB")
	if (modelVariant == "betweenItem") {
		ps$pContBetween = bquote("P"^"O")
	}
	
	ps$pContWithin = bquote("P"^"OW")
	if (modelVariant == "withinItem") {
		ps$pContWithin = bquote("P"^"O")
	}
	
	ps$pCatGuess = bquote("P"^"AG")
	
	ps$contSD = bquote(sigma^"O")
	ps$catSD = bquote(sigma^"A")
	ps$catSelectivity = bquote(sigma^"S")
	
	ps$catMu = bquote(mu)
	ps$catActive = bquote(nu)
	
	ps
}


#' Logit Transformations
#' 
#' The logit transformation is `log(p/(1 - p))`.
#' 
#' The inverse logit transformation is `exp(f)/(1 + exp(f)))`.
#' 
#' @param p A probability (value in the interval `[0,1]`) to transform to `(-Inf, Inf)`.
#' @param f A value in the interval `(-Inf, Inf)` to transform to `[0,1]`.
#' 
#' @return The transformed value.
#' 
#' @seealso The functions `qlogis` and `plogis` in the `stats` namespace..
#' 
#' @export
logit = function(p) {
	stats::qlogis(p)
}

#' @rdname logit
#' @export
logitInverse = function(f) { 
	stats::plogis(f)
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
#' distribution is parameterized in terms of "standard deviation",
#' the square root of the inverse of precision.
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
#' The Von Mises density function is parameterized in terms of precision in radians. The models in this package are isntead parameterized in terms of standard deviation (square root of the inverse of precision) in degrees. These functions go back and forth between the two representations of precision.
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
#' @param angles A vector of angles. If in radians, set `degrees` to `FALSE`.
#' @param weights A vector of weights.
#' @param degrees If `TRUE`, angles are treated as degrees. If `FALSE`, angles are treated as radians.
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
#' @param x Vector of values in degrees or radians.
#' @param y Vector of values in degrees or radians.
#' @param degrees If `TRUE`, x and y are treated as though they are in degrees. If `FALSE`, x and y are treated as being in radians.
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



#' Calculate Probabilities of Assignment to Categories
#' 
#' This function is the weights function, given in Equation 5 of the Appendix of Hardman, Vergauwe, and Ricker (2017).
#' 
#' @param study A scalar study angle, in degrees.
#' @param catMu A vector of category means, in degrees.
#' @param catSelectivity The categorical selectivity parameter, as standard deviation in degrees.
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





# For the following functions, the equation numbers correspond 
# to equations in the Appendix of Hardman, Vergauwe, and Ricker (2017).
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
#' @param dataType One of `"circular"` or `"linear"`.
#' @param catActive A vector of category active parameters. By default, all of the provided catMus are active.
#' @param muRange Required only if `dataType == "linear"`. A length 2 vector giving the lower and upper bounds of the catMu parameters. This should usually be the same as the range of the data.
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
		xs = seq(muRange[1], muRange[2], length.out = 360)
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
#' This function provides one way of producing colors from variants of the HSV color space by converting an input angle into an output angle using a linear transformation. See the equations related to producing colors in the Method section for Experiment 1 in Hardman, Vergauwe, and Ricker (2017).
#' 
#' @param angle An angle in degrees for which a color should be produced.
#' @param inPoints See description.
#' @param outPoints See description.
#' 
#' @return A color value.
#' 
#' @export
#' 
#' @examples
#' #Make a curried function that uses warpedHSVColorGeneratingFunction
#' #to produce colors like those used in Experiment 1 of Hardman, Vergauwe, and Ricker (2017).
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
#' For given parameter values, a set of data, and a model variant, this calculates the likelihood for each observation in the data set. To be clear, this function calculates likelihood, not log likelihood.
#' 
#' Note that if you don't want to provide any `catMu` parameters, do not set `catMu` to `NULL`. Rather, set it to a zero-length vector with `vector(length = 0)`. Do not provide `catActive` parameters. Rather, do not include inactive `catMu`s.
#' 
#' If `dataType` is `"circular"`, you should provide a value for `minSD`. 
#' If `dataType` is `"linear"`, you should provide a value for `responseRange`.
#' 
#' @param param A list of parameters just like that returned by [`getSingleIterationParameters`].
#' @param data A data.frame of the data for which you want to calculate the likelihood. Only the `study` and `response` columns are required, but the `pnum` and `cond` columns may be included.
#' @param modelVariant One of `"betweenItem"`, `"withinItem"`, and `"ZL"`.
#' @param dataType One of `"circular"` or `"linear"`.
#' @param responseRange A length 2 vector giving the theoretical minimum and maximum values of a response. Should be provided if `dataType` is `"linear"`.
#' @param minSD The minimum standard deviation of the Von Mises or normal distributions.
#' 
#' @return The provided `data` `data.frame` with an additional column named `likelihood` giving the likelihood for each observation.
#' 
#' @export
likelihood = function(param, data, modelVariant, dataType = "circular", 
											responseRange = NULL, minSD = NULL) 
{
	
	tr = list(config = list(modelVariant = modelVariant))
	
	allParam = getAllParams(tr, filter=TRUE)
	if (modelVariant != "ZL") {
		allParam = c(allParam, "catMu")
	}
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


valueIfNull = function(x, value) {
	if (is.null(x)) {
		x = value
	}
	x
}

getDefaultParametersWithConditionEffects = function(modelVariant) {
	pce = c("pMem", "contSD") #ZL and all other models
	
	if (modelVariant == "betweenAndWithin") {
		pce = c(pce, "pBetween", "pContBetween", "pContWithin")
		
	} else if (modelVariant == "betweenItem") {
		pce = c(pce, "pContBetween")
		
	} else if (modelVariant == "withinItem") {
		pce = c(pce, "pContWithin")
	}
	
	pce
}



#' Sample Data from Model with Specific Parameter Values
#' 
#' Samples data from the model given the provided parameter values. This is useful for observing the patterns of data generated by the model and for sampling from the posterior predictive distribution of the data.
#' 
#' @param study A vector of study angles in degrees.
#' @param param A list of parameter values. The values to include are `pMem`, `pBetween`, `pContBetween`, `pContWithin`, `pCatGuess`, `contSD`, `catMu`, `catSelectivity`, `catSD`. `catMu` should be a vector. If there are no categories, `catMu` should be `NULL`.
#' @param modelVariant The model variant, as a string. Should be one of "betweenItem", "withinItem", and "ZL".
#' @param dataType One of `"circular"` or `"linear"`.
#' @param responseRange A length 2 vector giving the theoretical minimum and maximum values of a response. Should be provided if `dataType` is `"linear"`.
#' 
#' @return A data frame containing the `study` angles, the sampled `response` angles, the response `type` (e.g. continuous memory response), and, if the response was categorical in nature, the category it was from (`cat`).
#'  
#' @export
sampleDataFromModel = function(study, param, modelVariant, dataType = "circular", responseRange = NULL) {
	
	trials = length(study)
	
	if (modelVariant == "betweenItem") {
		param$pBetween = 1
	} else if (modelVariant == "withinItem") {
		param$pBetween = 0
	} else if (modelVariant == "ZL") {
		param$pBetween = 1
		param$pContBetween = 1
		param$pCatGuess = 0
		param$catMu = NULL
	}
	
	realzationFunction = NULL
	unifGuessFunction = NULL
	if (dataType == "circular") {
		
		realzationFunction = function(mu, sd) {
			CatContModel::rvmd(1, mu, sd)
		}
		unifGuessFunction = function() {
			stats::runif(1, 0, 360)
		}
		
	} else if (dataType == "linear") {
		
		realzationFunction = function(mu, sd) {
			msm::rtnorm(1, mu, sd, responseRange[1], responseRange[2])
		}
		unifGuessFunction = function() {
			stats::runif(1, responseRange[1], responseRange[2])
		}
	}
	
	
	catCount = length(param$catMu)
	
	
	response = rep(0, trials)
	
	cat = rep(0, trials)
	type = rep("none", trials)
	
	for (i in 1:trials) {
		
		inMemory = (stats::rbinom(1, 1, param$pMem) == 1)
		
		if (!inMemory) {
			
			isCatGuess = (stats::rbinom(1, 1, param$pCatGuess) == 1)
			
			if (isCatGuess && (catCount > 0)) {
				cat[i] = sample(1:catCount, size=1)
				
				response[i] = realzationFunction(param$catMu[cat[i]], param$catSD)
				type[i] = "catGuess"
			} else {
				response[i] = unifGuessFunction()
				type[i] = "unifGuess"
			}
			
		} else {
			
			
			# Pick an error for the continuous representation
			
			contLocation = realzationFunction(study[i], param$contSD)
			
			if (dataType == "circular") {
				contLocation = contLocation %% 360
			}
			
			if (catCount == 0) {
				response[i] = contLocation
				type[i] = "continuous"
				
			} else {
				
				# Pick a category for the color
				catWeights = categoryWeightsFunction(study[i], param$catMu, param$catSelectivity, dataType = dataType)
				catEgory = sample(1:catCount, size=1, prob=catWeights)
				cat[i] = catEgory
				
				
				# Get the categorical response location
				catLocation = realzationFunction(param$catMu[catEgory], param$catSD)
				if (dataType == "circular") {
					catLocation = catLocation %% 360
				}
				
				isBetween = (stats::rbinom(1, 1, param$pBetween) == 1)
				
				if (isBetween) {
					# This is a between response
					isContinuous = (stats::rbinom(1, 1, param$pContBetween) == 1)
					
					if (isContinuous) {
						response[i] = contLocation
						type[i] = "continuous"
					} else {
						response[i] = catLocation
						type[i] = "categorical"
					}
				} else {
					# This is a within response
					
					locMix = NA
					if (dataType == "circular") {
						locMix = circMean(c(contLocation, catLocation), c(param$pContWithin, 1 - param$pContWithin))
						locMix = locMix %% 360
					} else if (dataType == "linear") {
						locMix = param$pContWithin * contLocation + (1 - param$pContWithin) * catLocation
					}
					
					response[i] = locMix
					
					type[i] = "within"
				}
			}
		}
		
	}
	data.frame(study=study, response=response, type=type, cat=cat)
}


# each column of mat corresponds to one row of design.
# designCols controls which columns of design are kept.
reshapeMatrixToDF = function(mat, design, designCols = names(design)) {
	# stringsAsFactors = FALSE is in case mat is character
	df = data.frame(x = as.vector(mat), stringsAsFactors = FALSE)
	for (n in designCols) {
		df[ , n ] = rep(design[ , n ], each=nrow(mat))
	}
	df
}


getSubsampleIterationsToRemove = function(totalIterations, subsamples, subsampleProportion) {
	if (subsamples < 1) {
		stop("You need to use at least one subsample.")
	}
	
	independentSubsamples = FALSE
	if (is.null(subsampleProportion)) {
		independentSubsamples = TRUE
		subsampleProportion = 1 / subsamples
	}
	
	subsampleProportion = min( max(subsampleProportion, 0), 1 )
	
	subsampleIterationsToRemove = list()
	
	if (independentSubsamples) {
		
		shuffledIterations = sample(1:totalIterations, totalIterations, replace=FALSE)
		
		for (sub in 1:subsamples) {
			
			iterationsToUse = floor(subsampleProportion * totalIterations)
			indicesToUse = ((sub - 1) * iterationsToUse + 1):(sub * iterationsToUse)
			subsampleIterationsToRemove[[sub]] = shuffledIterations[-indicesToUse]
			
		}
	} else {
		
		for (sub in 1:subsamples) {
			iterationsToRemove = round((1 - subsampleProportion) * totalIterations, 0)
			subsampleIterationsToRemove[[sub]] = sample(1:totalIterations, iterationsToRemove, replace = FALSE)
		}
		
	}
	
	subsampleIterationsToRemove
}


#' Credible Intervals for Manifest Condition Effect Priors
#' 
#' Which is to say, what prior beliefs do you have about the distribution of condition effects?
#' This function helps deal with the fact that the interpretation of priors on condition effects 
#' 1) depends on the value of the participant parameters and
#' 2) need to be translated from the latent space to the manifest space.
#' 
#' Given a vector of participant parameter values and a prior scale for the condition effect parameter,
#' this:
#' 1) Samples from the condition effect prior (0 centered).
#' 2) Calculates the manifest participant parameter values for each sample from the condition effect prior and
#' 3) Calculates the prior credible interval for the manifest parameter values.
#' This information can help you choose a prior scale value that captures your beliefs about how much the conditions differ from one another (or, more exactly, how much the non-cornerstone conditions differ from the cornerstone condition).
#' 
#' The left plot shows the prior median and the lower and upper credible interval bounds, giving an overall picture of the prior credible interval across the participant parameter values.
#' 
#' The right plot gives the total width of the credible interval and the directional widths from the median to the upper and lower bounds.
#'  
#' @param param The name of the parameter (e.g. `pMem`).
#' @param p_i A vector of manifest participant parameter values. Typically a series of numbers (used for x-axis in plotting). For probability parameters, use something like `seq(0, 1, 0.025)` for the whole range of probabilities. For SD parameters, use something like `seq(0, 40, 1)`.
#' @param ce_scale The scale of the Cauchy prior on the condition effect parameter.
#' @param cip Proportion of the prior inside of the credible interval.
#' @param n Number of samples to take from the prior on the credible interval. Use more for a more accurate approximation.
#' @param minSD Only for standard deviation parameters. The minimum standard deviation (i.e. `config$minSD`).
#' @param plot Whether to make plots.
#' 
#' @return Invisibly, a `data.frame` containing columns:
#' * `p_i`: The participant parameter value (copied from the `p_i` argument).
#' * `lower`, `upper`: Lower and upper bounds of the credible interval.
#' * `median`: The median of the prior.
#' * `lowerW`, `upperW`: The distance from the median to the lower and upper credible intervals.
#' 
#' @export
conditionEffectPriorCredibleInterval = function(param, p_i, ce_scale, cip = 0.95, n = 1e6, minSD = 1, plot = TRUE, doMFRow = TRUE) {
	
	qp = c((1 - cip) / 2, 0.5, (1 + cip) / 2)
	
	# Use same samples for each value of x
	ce = rcauchy(n, 0, ce_scale)
	
	res = list(config = list(minSD = minSD)) # Fake res for getting transformations
	trans = getParameterTransformation(res, param)
	inverse = getParameterTransformation(res, param, inverse = TRUE)
	
	df = data.frame(p_i = p_i)
	
	for (i in 1:length(p_i)) {
		
		manifest = trans(inverse(p_i[i]) + ce)
		
		qs = quantile(manifest, qp)
		df$lower[i] = qs[1]
		df$median[i] = qs[2]
		df$upper[i] = qs[3]
	}
	
	df$lowerW = df$median - df$lower
	df$upperW = df$upper - df$median
	df$totalW = df$upper - df$lower
	
	if (plot) {
		
		lowCol = "red"
		upCol = "blue"
		
		if (doMFRow) {
			par(mfrow=c(1, 2))
		}
		
		ylimL = range(df[ , c("lower", "median", "upper")])
		if (param %in% getProbParams(NULL)) {
			ylimL = c(0, 1)
		}
		
		plot(df$p_i, df$median, ylim=ylimL, type='l', xlab=param, ylab="Manifest Parameter Value")
		lines(df$p_i, df$lower, col=lowCol)
		lines(df$p_i, df$upper, col=upCol)
		
		ylimR = range(df[ , c("lowerW", "totalW", "upperW")])
		
		plot(df$p_i, df$totalW, type = 'l', ylim=ylimR, xlab=param, ylab="Credible Interval Width")
		lines(df$p_i, df$upperW, col=upCol)
		lines(df$p_i, df$lowerW, col=lowCol)
		legend("bottom", legend = c("upper", "total", "lower"), col=c(upCol, "black", lowCol), lty=1)
		
		if (doMFRow) {
			par(mfrow=c(1, 1))
		}
	}
	
	invisible(df)
}

#' Plot Manifest Condition Effect Prior Histogram
#' 
#' @param param A parameter name.
#' @param p_i A manifest participant parameter value. E.g., for probability parameters, give a probability.
#' @param ce_scale The scale of the Cauchy prior on the condition effect parameter.
#' @param sdCutoff For standard deviation parameters, extremely large sample values are common, so for plotting purposes the plot has to be cut off somewhere. `sdCutoff` sets the cutoff.
#' @param n Number of samples to take from the prior on the credible interval. Use more for a more accurate approximation.
#' @param minSD Only for standard deviation parameters. The minimum standard deviation (i.e. `config$minSD`).
#' 
#' @return Invisibly, the vector of manifest parameter samples from the condition effect prior that were used for plotting (so some samples are cut off for SD parameters; see the `sdCutoff` argument).
conditionEffectPriorHist = function(param, p_i, ce_scale, sdCutoff=50, n=1e6, minSD=1) {
	
	ce = rcauchy(n, 0, ce_scale)
	
	res = list(config = list(minSD = minSD))
	trans = getParameterTransformation(res, param)
	inverse = getParameterTransformation(res, param, inverse=TRUE)
	
	manifest = trans(inverse(p_i) + ce)
	if (param %in% getSdParams(NULL)) {
		manifest = manifest[manifest < sdCutoff]
	}
	
	main = paste0(param, " = ", p_i, ", scale = ", ce_scale)
	hist(manifest, 
			 xlab=paste0("Manifest ", param), 
			 ylab="Prior Density", 
			 main=main, prob=TRUE)
	
	invisible(manifest)
}

