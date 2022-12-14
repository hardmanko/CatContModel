
#' Open CatContModel Package Introduction Manual
#' 
#' Opens the manual for the package. 
#' 
#' @export
CatContModelManual = function() {
	utils::vignette("Introduction", "CatContModel")
}


valueIfNull = function(x, value) {
  if (is.null(x)) {
    x = value
  }
  x
}


substituteValues = function(x, xvals, subvals) {
	y = rep(subvals[1], length(x)) #to get the type right
	for (i in 1:length(xvals)) {
		y[x == xvals[i]] = subvals[i]
	}
	y
}

#' Get a list of parameter symbols for mathy plotting
#'
#' @param modelVariant String giving the name of the model variant to use. See [`makeModelConfig`] for choices.
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
#' The inverse logit transformation is `exp(q)/(1 + exp(q)))`.
#' 
#' @param p Probabilities in the interval `[0,1]` to transform to `(-Inf, Inf)`.
#' @param q Quantiles in the interval `(-Inf, Inf)` to transform to `[0,1]`.
#' 
#' @return The transformed values.
#' 
#' @seealso The functions `qlogis` and `plogis` in the `stats` namespace..
#' 
#' @export
logit = function(p) {
	stats::qlogis(p)
}

#' @rdname logit
#' @export
logitInverse = function(q) { 
	stats::plogis(q)
}

#' Convert Between Degrees and Radians
#' 
#' Degrees to radians is `deg * pi / 180`. Radians to degrees is `rad * 180 / pi`.
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

#' Convert Between Precision in Radians and Standard Deviation in Degrees
#' 
#' The Von Mises density function is typically parameterized in terms of precision in radians. 
#' The models in this package are instead parameterized in terms of standard deviation (square root of the inverse of precision) in degrees. 
#' These functions go back and forth between the two representations of precision.
#' 
#' If precision is `kappa` (radians) and `sd` is standard deviation (degrees):
#' + `precRad_to_sdDeg`: `sd = 1 / sqrt(kappa) * 180 / pi`. 
#' + `sdDeg_to_precRad`: `kappa = 1 / (sd * pi / 180)^2`
#' 
#' @param sdDeg A standard deviation in degrees.
#' @param precRad A precision in radians (i.e. `kappa`).
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


#' Von Mises Distribution Functions
#' 
#' These functions are wrappers that parameterize the von Mises distribution in terms of degrees. 
#' 
#' Instead of precision, the von Mises distribution is parameterized in terms of "standard deviation",
#' the square root of the inverse of precision.
#' See [`sdDeg_to_precRad`] and [`precRad_to_sdDeg`].
#' 
#' For Von Mises density and sampling functions parameterized in terms 
#' of radians and precision, see the CircStats package, in particular 
#' the functions `dvm` and `rvm`, which are used internally by `dVonMises` and `rVonMises`.
#' Or set the `degrees` argument to `FALSE`.
#' 
#' @param n The number of realizations to draw.
#' @param xs A vector of quantiles.
#' @param mu The center of the von Mises distribution.
#' @param sd The standard deviation in degrees. If `degrees=FALSE`, `sd` is precision in radians (`kappa`).
#' @param log If true, the log density is returned.
#' @param degrees If `FALSE`, functions use radians instead of degrees. That means that `xs` and `mu` are radians and `sd` is precision (`kappa`).
#' @param impl Von Mises density implementation. The first value is used. LUT and noLUT use the C++ implementation.
#' 
#' @return A realization (`rVonMises`) or a density (`dVonMises`).
#' 
#' @export
dVonMises = function(xs, mu, sd, log=FALSE, degrees=TRUE, impl=c("CircStats", "LUT", "noLUT")) {
  
  impl = impl[1]
  
  if (impl == "CircStats") {
    
    if (degrees) {
      xs = d2r(xs)
      mu = d2r(mu)
      sd = sdDeg_to_precRad(sd)
    }
    
    dens = CircStats::dvm(xs, mu, sd)
    
    if (degrees) {
      # dvm() gives a density for arguments in radians. The width of
      # the space is 2*PI radians but 360 degrees. This means that the
      # same Von Mises distribution will have higher density in radians 
      # than degrees because the numeric circumference is lower for radians
      # than degrees. Rescale the density to be correct for degrees by multiplying
      # by 2*pi/360 = pi/180.
      dens = dens * pi / 180
    }
    
    if (log) {
      dens = log(dens)
    }
    
  } else if (impl == "LUT") {
    
    dens = CCM_CPP_dVonMises(xs, mu, sd, log=log, degrees=degrees, useLUT=TRUE)
    
  } else if (impl == "noLUT") {
    
    dens = CCM_CPP_dVonMises(xs, mu, sd, log=log, degrees=degrees, useLUT=FALSE)
    
  }
  
  dens
}

#' @rdname dVonMises
#' @export
rVonMises = function(n, mu, sd, degrees=TRUE) {
  
  if (degrees) {
    mu = d2r(mu)
    kappa = sdDeg_to_precRad(sd)
  } else {
    kappa = sd
  }
  
  sampled = CircStats::rvm(n, mu, kappa)
  
  if (degrees) {
    sampled = r2d(sampled)
  }
  
  sampled
}



#' Calculate an Optionally Weighted Circular Mean
#' 
#' @param angles A vector of angles. If in radians, set `degrees` to `FALSE`.
#' @param weights A vector of weights. If weights is length 0 or 1, equal weights are used.
#' @param degrees If `TRUE`, angles are treated as degrees. If `FALSE`, angles are treated as radians.
#' 
#' @return The circular mean of the angles.
#' 
#' @export
#' 
circMean = function(angles, weights=1, degrees=TRUE) {
  if (length(weights) <= 1) {
    weights = rep(1, length(angles))
  }
  weights = weights / sum(weights)
  
  CatContModel:::CCM_CPP_circMean(angles, weights, degrees)
}






# For the following functions, the equation numbers correspond 
# to equations in the Appendix of Hardman, Vergauwe, and Ricker (2017).
# k is a vector of active category indices.
h_eq20 = function(mu_k, nu_k, mu_kp, nu_kp, catMuPriorSD, dataType) {
	
	if (dataType == "circular") {
		#numDens = dvmd(mu_k - mu_kp, 0, catMuPriorSD)
		numDens = dVonMises(mu_k - mu_kp, 0, catMuPriorSD, degrees=TRUE)
		
		#denDens = dvmd(0, 0, catMuPriorSD)
		denDens = dVonMises(0, 0, catMuPriorSD, degrees=TRUE)
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
	
	mvUsedParam = getParamNames(modelVariant, types=c("prob", "sd"))
	if (modelVariant != "ZL") {
	  mvUsedParam = c(mvUsedParam, "catMu")
	}
	
	missingParamNames = mvUsedParam[ !(mvUsedParam %in% names(param)) ]
	if (length(missingParamNames) > 0) {
		stop(paste0("Some needed parameters were missing: ", paste(missingParamNames, collapse=", ")))
	}
	
	# Add NA parameters not used by the modelVariant to the list
	unusedParam = getParamNames(NULL)
	unusedParam = unusedParam[ unusedParam != "catActive"] # Exclude catActive
	
	unusedParam = unusedParam[ !(unusedParam %in% mvUsedParam) ]
	for (up in unusedParam) {
	  param[[up]] = NA
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



# Used to convert catMu between degrees and radians (and other operations)
convertCatMuUnits = function(results, conversionFun) {
  
  for (p in unique(results$data$pnum)) {
    for (i in 1:results$config$maxCategories) {
      cmName = paste("catMu[", p, ",", i-1, "]", sep="")
      results$posteriors[[ cmName ]] = conversionFun(results$posteriors[[ cmName ]])
    }
  }
  
  results
}


#' Plot Data by Condition and Participant
#' 
#' Make scatterplots and (optionally) histograms of data.
#' 
#' @param data Data to plot. Formatted with `pnum`, `cond`, `study`, and `response` columns.
#' @param byCond Make plots for each condition (collapsing across participants).
#' @param byPnum Make plots for each participant (collapsing across conditions).
#' @param pnumByCond Make plots for each pnum * cond cell.
#' @param pdfFile Optional. Name of pdf file to plot to.
#' @param pdfPanelWH A length 2 vector giving width and height in inches of panels in pdf file.
#' @param histograms If `TRUE`, response histograms are plotted along with scatterplots.
#' 
#' @family plotting functions
#' @export
plotData = function(data, allData=TRUE, byCond=TRUE, byPnum=TRUE, pnumByCond=TRUE, plotRange=NULL, pdfFile=NULL, pdfPanelWH=c(6,6), histograms=TRUE) {
  
  
  #whichPlots = c("all", "cond", "pnum", "pnum*cond")
  
  plotToFile = !is.null(pdfFile) && (pdfFile != "")
  
  if (histograms) {
    pdfPanelWH[1] = pdfPanelWH[1] * 2
  }
  
  if (plotToFile) {
    grDevices::pdf(pdfFile, width=pdfPanelWH[1], height=pdfPanelWH[2])
  }
  
  if (is.null(plotRange)) {
    plotRange = range(c(data$study, data$response))
    
    # If almost circular, make circular
    if (abs(plotRange[1] - 0) < 2 && abs(plotRange[2] - 360) < 2) {
      plotRange = c(0, 360)
    }
  }
  
  alphaFun = function(nObs) {
    alpha = 6 / sqrt(nObs)
    max(0, min(alpha, 1))
  }
  
  if (histograms) {
    graphics::par(mfrow=c(1,2))
  }
  
  if (allData) {
    scatterplotWithColorBars(data, alpha=alphaFun(nrow(data)), xlim=plotRange, ylim=plotRange)
    graphics::title("All cond and pnum")
    
    if (histograms) {
      graphics::hist(data$response, breaks=60, main="response")
    }
  }
  
  if (byCond) {
    for (cond in unique(data$cond)) {
      condData = data[ data$cond == cond, ]
      scatterplotWithColorBars(condData, alpha=alphaFun(nrow(condData)), xlim=plotRange, ylim=plotRange)
      graphics::title(paste0("cond=", cond, " (all pnums)"))
      
      if (histograms) {
        graphics::hist(condData$response, breaks=60, main="response")
      }
    }
  }
  
  if (byPnum) {
    for (pnum in unique(data$pnum)) {
      pnumData = data[ data$pnum == pnum, ]
      scatterplotWithColorBars(pnumData, alpha=alphaFun(nrow(pnumData)), xlim=plotRange, ylim=plotRange)
      graphics::title(paste0("pnum=", pnum, " (all conds)"))
      
      if (histograms) {
        graphics::hist(pnumData$response, breaks=60, main="response")
      }
    }
  }
  
  if (pnumByCond) {
    
    for (pnum in unique(data$pnum)) {
      for (cond in unique(data$cond)) {
        pcData = data[ data$pnum == pnum & data$cond == cond, ]
        scatterplotWithColorBars(pcData, alpha=alphaFun(nrow(pcData)), xlim=plotRange, ylim=plotRange)
        graphics::title(paste0("pnum*cond=", pnum, "*", cond))
        
        if (histograms) {
          graphics::hist(pcData$response, breaks=60, main="response")
        }
      }
    }
  }
  
  if (plotToFile) {
    grDevices::dev.off()
  }
  
}


