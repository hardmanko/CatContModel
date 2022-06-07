
# paramGroup like pMem, contSD, catMu, etc.
getOrderedMatchingParameterNames = function(res, paramGroup) {
  
  paramNames = names(res$posteriors)
  
  # Find the parameters and put them in order
  matchingParam = grepl(paste0("^", paramGroup), paramNames)
  matchingParam = paramNames[ matchingParam ]
  
  # put .mu, .var, and _cond first
  muLoc = which(grepl("\\.mu$", matchingParam))
  varLoc = which(grepl("\\.var$", matchingParam))
  condLocs = which(grepl("_cond", matchingParam))
  paramOrder = c(muLoc, varLoc, condLocs)
  
  allLocs = 1:length(matchingParam)
  participantLocs = allLocs[ !(allLocs %in% paramOrder) ]
  
  paramOrder = c(paramOrder, participantLocs)
  
  matchingParam[ paramOrder ]
  
}

#' Plots of Posterior Parameter Chains
#' 
#' Makes plots, preferably to a pdf file, of posterior parameter chains for purposes of visual convergence diagnostics.
#' There are other uses for these plots, including model validation (are the parameters behaving reasonably?).
#' 
#' @param res A results object.
#' @param filename The name of a pdf file to plot to.
#' @param whichIterations Which iterations to plot.
convergencePlots = function(res, filename=NULL, whichIterations=NULL) {
  convergencePlots_parallel(list(res), filename=filename, whichIterations=whichIterations)
}



# mega function to do all parameters in one file
convergencePlots_parallel = function(resList, filename=NULL, whichIterations=NULL) {
  
  titlePagePlot = function(text, x=0.5, y=0.5, cex=3) {
    par(mfrow=c(1,1))
    plot(c(0,1), c(0,1), type='n', axes=FALSE, xlab="", ylab="")
    text(x, y, text, cex=cex)
  }
  
  usedPnums = resList[[1]]$pnums
  
  if (!is.null(filename)) {
    grDevices::pdf(filename, width = 8, height = 4)
  }

  paramNameGroups = getAllParams(resList[[1]], filter=TRUE)
  
  # Main title page
  origMar = par()$mar
  par(xpd=TRUE, mar=c(0,0,0,0))
  titlePagePlot("Convergence Diagnostic Plots", y=0.9)
  
  printedParamNameGroups = paramNameGroups[1:min(length(paramNameGroups), 3)]
  if (length(paramNameGroups) > 3) {
    printedParamNameGroups = c(printedParamNameGroups, "etc.")
  }
  
  labels = c("For all standard parameters, posterior chains and density are plotted. ", 
             "   Population mean and variance: param.mu, param.var",
             "   Effects of task condition: param_cond[cond]",
             "   Participant level parameters: param[pnum]",
             paste0("Standard params: ", paste(printedParamNameGroups, collapse=", "))
             )
  
  if (resList[[1]]$config$modelVariant != "ZL") {
    labels = c(labels, "Category parameters (catMu, catActive) follow with a different format.")
  }
  
  text(x=0.05, y=seq(0.7, 0.1, length.out=length(labels)), labels=labels, cex=1.2, adj=0)
  par(xpd=FALSE, mar=origMar)
  
  # standard params
  for (i in 1:length(paramNameGroups)) {
    
    titlePagePlot(paramNameGroups[i])
    
    # TODO: Select participants
    matchingParam = getOrderedMatchingParameterNames(resList[[1]], paramNameGroups[i])
    
    #graphics::par(mfrow=c(1,2))
    for (i in 1:length(matchingParam)) {
      convergencePlotStandard_parallel(resList, matchingParam[i], type=c("chain", "density"), 
                                       whichIterations=whichIterations, combinedDens=FALSE, panelize=TRUE)
    }
    
  }
  
  

  if (resList[[1]]$config$modelVariant != "ZL") {
  
    # category params are done for each results object separately so participants have their chains together
    
    # title page
    if (length(resList) > 1) {
      titlePagePlot("Category Parameters", y=0.75)
      labels = c("For each participant:", "Multiple chains: catMu density, catActive iterations", "Each chain: catMu scatter, catActive iterations")
      text(x=0.5, y=c(0.45, 0.3, 0.15), labels=labels, cex=1.5)
    } else {
      titlePagePlot("Category Parameters", y=0.5)
    }
    
    # For each participant, plot the multi chain results, then each individual chain
    for (i in 1:length(usedPnums)) {
      
      if (length(resList) > 1) {
        #convergencePlotCategory_parallel(resList, pnums=usedPnums[i], panelize=TRUE)
        par(mfrow=c(1,1))
        convergencePlot_catMuDensity_parallel(resList, pnums=usedPnums[i], whichIterations=whichIterations)
        convergencePlot_catActive_parallel(resList, pnums=usedPnums[i], whichIterations=whichIterations)
      }
      
      for (resInd in 1:length(resList)) {
        #convergencePlotCategory(resList[[resInd]], pnums=usedPnums[i], panelize=TRUE)
        
        par(mfrow=c(1,2))
        convergencePlot_catMuScatter(resList[[resInd]], pnum=usedPnums[i], whichIterations=whichIterations)
        convergencePlot_catActive_parallel(list(resList[[resInd]]), pnums=usedPnums[i], whichIterations=whichIterations)
        
      }
    }
  }
    
  if (!is.null(filename)) {
    dev.off()
  }
  
}

convergencePlotStandard = function(res, param, type=c("chain", "density"), whichIterations=NULL, panelize=TRUE) {
  convergencePlotStandard_parallel(list(res), param=param, type=type, whichIterations=whichIterations, combinedDens=FALSE, panelize=panelize)
}


# For non-category parameters (all but catMu and catActive)
# type %in% chain, density
convergencePlotStandard_parallel = function(resList, param, type=c("chain", "density"), whichIterations=NULL, combinedDens=FALSE, panelize=TRUE) {
  
  # check param name
  if (!(param %in% names(resList[[1]]$posteriors)) || 
      grepl("catMu", param, fixed=TRUE) || 
      grepl("catActive", param, fixed=TRUE)) 
  {
    warning(paste0("Invalid parameter name: \"", param, "\""))
    return()
  }
  
  # shared setup
  mm = matrix(nrow=resList[[1]]$config$iterations, ncol=length(resList))
  for (i in 1:length(resList)) {
    mm[,i] = resList[[i]]$posteriors[[param]]
  }
  
  whichIterations = valueIfNull(whichIterations, 1:resList[[1]]$config$iterations)
  whichIterations = whichIterations[ whichIterations %in% 1:nrow(mm) ]
  mm = mm[whichIterations,,drop=FALSE]
  
  valRange = range(mm)
  
  col = rainbow(ncol(mm))
  
  
  # Clean type
  type = unique(type[ type %in% c("chain", "density") ])
  if (panelize && length(type) == 2) {
    graphics::par(mfrow=c(1,2))
  }
  
  # chain
  if ("chain" %in% type) {
    base::plot(x=range(whichIterations), y=valRange, type='n', xlab="Iteration", ylab=param)
    for (i in 1:ncol(mm)) {
      graphics::lines(x=whichIterations, y=mm[,i], col=col[i])
    }
    graphics::title(paste0("Chain for ", param))
  }
  
  
  # density
  paramIsConstant = valRange[1] == valRange[2]
  if ("density" %in% type) {
    require(polspline)
    
    # Don't try to do density for constant parameters
    if (paramIsConstant) {
      
      densXs = valRange[1] + c(-1, -0.01, 0, 0.01, 1)
      densMat = matrix(0, nrow=length(densXs), ncol=ncol(mm))
      densMat[3,] = 1
      
    } else {
      densXs = seq(valRange[1], valRange[2], length.out = 100)
      
      densMat = matrix(nrow=length(densXs), ncol=ncol(mm))
      for (i in 1:ncol(mm)) {
        ls = logspline(mm[,i])
        densMat[,i] = dlogspline(densXs, ls)
      }
    }
    
    base::plot(x=valRange, y=c(0,max(densMat)), type='n', xlab=param, ylab="Density")
    for (i in 1:ncol(mm)) {
      graphics::lines(x=densXs, y=densMat[,i], col=col[i])
    }
    graphics::title(paste0("Density for ", param))
    
    # Add combined density line
    if (combinedDens && !paramIsConstant) {
      totalDensLS = polspline::logspline(as.vector(mm))
      graphics::lines(x=densXs, y=polspline::dlogspline(densXs, totalDensLS), col="black", lwd=2, lty=2)
    }
  }
}




# single chain, single participant
convergencePlot_catMuScatter = function(res, pnum, whichIterations=NULL) {
  
  whichIterations = valueIfNull(whichIterations, 1:res$config$iterations)
  
  catCol = rainbow(res$config$maxCategories)
  alphaHex = gettextf("%X", 51)
  catCol = paste0(catCol, alphaHex)
  
  pMats = convertPosteriorsToMatrices(res, param=c("catMu", "catActive"))
  
  pInd = which(pnum == res$pnums)
  
  cmi = pMats$catMu[pInd,,]
  cai = pMats$catActive[pInd,,]
  
  
  ylim = c(0, 360)
  if (res$config$dataType == "linear") {
    ylim = res$config$catMuRange
  }
  
  base::plot(x=range(whichIterations), y=ylim, axes=FALSE, type='n', xlab="Iteration", ylab="Category Location")
  graphics::box()
  graphics::axis(1)
  graphics::axis(2, at=seq(0, 360, 45))
  graphics::title(paste0("Active catMu for pnum:", pnum))
  
  for (j in 1:dim(cmi)[1]) {
    cmij = cmi[j,whichIterations]
    caij = cai[j,whichIterations]
    
    cmij = CatContModel::clampAngle(cmij, pm180 = FALSE, degrees = TRUE)
    
    cmij[ caij == 0 ] = -100 # This is such a gross hack
    
    graphics::points(x=whichIterations, y=cmij, pch=18, col=catCol[j]) # pch=20
  }
  
}


# All participants in one density plot.
# Only active categories
convergencePlot_catMuDensity_parallel = function(resList, pnums=NULL, catMuPrecision=2, whichIterations=NULL) {
  nRes = length(resList)
  pnums = valueIfNull(pnums, resList[[1]]$pnums)
  chainColors = grDevices::rainbow(nRes)
  
  whichIterations = valueIfNull(whichIterations, 1:resList[[1]]$config$iterations)
  
  allPlotData = NULL
  

  for (i in 1:nRes) {
    post = convertPosteriorsToMatrices(resList[[i]], param = c("catMu", "catActive"))
    
    
    includedPnumInd = which(dimnames(post$catMu)[[1]] %in% pnums)
    post$catMu = post$catMu[includedPnumInd,,whichIterations]
    
    includedPnumInd = which(dimnames(post$catActive)[[1]] %in% pnums)
    post$catActive = post$catActive[includedPnumInd,,whichIterations]
    
    
    plotData = catMu_getPlotData(post$catMu, post$catActive, 
                                 precision = catMuPrecision, 
                                 dataType = resList[[i]]$config$dataType, 
                                 responseRange = resList[[i]]$config$responseRange, 
                                 colorGeneratingFunction = resList[[i]]$config$colorGeneratingFunction)
    
    plotData$chain = i
    allPlotData = rbind(allPlotData, plotData)
    
  }
  
  catMu_plotSetup(allPlotData, resList[[1]]$config$dataType) #, ylim_1 = colorBarYlim[2])
  
  if (length(pnums) > 1) {
    graphics::title("Active catMu for multiple pnums")
  } else {
    graphics::title(paste0("Active catMu for pnum:", pnums))
  }
  
  for (i in 1:nRes) {
    plotData = allPlotData[ allPlotData$chain == i, ]
    catMu_plotData(plotData$x, plotData$y, col = chainColors[i], type = "line")
  }
  
  invisible(allPlotData)
}


# All participants in one chain plot
# plots lines, no density
# type = c("line", "point")
convergencePlot_catActive_parallel = function(resList, pnums=NULL, type="line", whichIterations=NULL) {
  
  nRes = length(resList)
  #nIterations = resList[[1]]$config$iterations
  pnums = valueIfNull(pnums, resList[[1]]$pnums)
  chainColors = grDevices::rainbow(nRes)
  
  whichIterations = valueIfNull(whichIterations, 1:resList[[1]]$config$iterations)
  
  # per chain, count of active categories averaged across participants
  catActivePlotData = matrix(0, nrow=length(whichIterations), ncol=nRes)
  
  for (i in 1:nRes) {
    post = convertPosteriorsToMatrices(resList[[i]], param = "catActive")
    
    includedPnumInd = which(dimnames(post$catActive)[[1]] %in% pnums)
    post$catActive = post$catActive[includedPnumInd,,whichIterations]
    
    
    if (length(pnums) > 1) {
      # sum for each participant
      caSum = apply(post$catActive, c(1,3), sum)
      # mean of participants
      catActivePlotData[,i] = apply(caSum, 2, mean)
    } else {
      # sum for the single participant
      catActivePlotData[,i] = apply(post$catActive, 2, sum)
    }
    
  }
  
  yRange = range(catActivePlotData)

  base::plot(range(whichIterations), yRange, type='n', xlab="Iteration", ylab="Number of active categories")
  
  if (length(pnums) > 1) {
    graphics::title("Mean of sum of catActive for multiple pnums")
  } else {
    graphics::title(paste0("Sum of catActive for pnum:", pnums))
  }
  
  for (i in 1:nRes) {
    if (type == "line") {
      graphics::lines(x=whichIterations, y=catActivePlotData[,i], col=chainColors[i], lwd=1)
    } else if (type == "point") {
      graphics::points(x=whichIterations, y=catActivePlotData[,i], col=chainColors[i], pch=18)
    }
  }
  
  invisible(catActivePlotData)
}


