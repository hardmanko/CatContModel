
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

titlePagePlot = function(text, x=0.5, y=0.5, cex=3, setPar=TRUE, resetPar=TRUE) {
  
  opar = graphics::par()
  
  if (setPar) {
    graphics::par(xpd=TRUE, mar=c(0,0,0,0), mfrow=c(1,1))
    #graphics::par(mfrow=c(1,1))
  } else {
    resetPar = FALSE
  }
  
  graphics::plot(c(0,1), c(0,1), type='n', axes=FALSE, xlab="", ylab="")
  graphics::text(x, y, text, cex=cex)
  
  if (resetPar) {
    graphics::par(xpd=opar$xpd, mar=opar$mar, mfrow=opar$mfrow)
  }
  
  opar
}

#' Plots of Posterior Parameter Chains
#' 
#' Makes plots, preferably to a pdf file, of posterior parameter chains for purposes of visual convergence diagnostics.
#' There are other uses for these plots, including model validation (are the parameters behaving reasonably?).
#' Title pages with
#' 
#' @param res A results object. Both WP and Parallel types are supported.
#' @param pdfFile The name of a pdf file to plot to. If the file extension is not ".pdf", the ".pdf" extension will be appended.
#' @param whichIterations Numeric vector. Which iterations to plot. If `NULL`, all iterations are plotted.
#' 
#' @export
convergencePlots = function(res, pdfFile=NULL, whichIterations=NULL) {
  
  if (resultIsType(res, "WP")) {
    # Make fake Parallel object.
    res = list(chains = list(res), parallelConfig = makeParallelConfig(1))
  } else if (resultIsType(res, "BP")) {
    # TODO: What about making separate files for each group?
    stop("BP results objects are not supported.")
  }

  convergencePlots_parallel(res, pdfFile=pdfFile, whichIterations=whichIterations)
}



# mega function to do all parameters in one file
convergencePlots_parallel = function(parRes, pdfFile=NULL, whichIterations=NULL) {
  
  usedPnums = parRes$chains[[1]]$pnums
  
  if (!is.null(pdfFile)) {
    if (!grepl("\\.pdf$", pdfFile)) {
      pdfFile = paste0(pdfFile, ".pdf")
    }
    
    grDevices::pdf(pdfFile, width = 8, height = 4)
  }

  standardParams = getParamNames(parRes$chains[[1]]$config$modelVariant, types=c("prob", "sd"))
  
  # Main title page
  opar = titlePagePlot("Convergence Diagnostic Plots", y=0.9, setPar=TRUE, resetPar = FALSE)
  
  printedParamNames = standardParams[1:min(length(standardParams), 3)]
  if (length(standardParams) > 3) {
    printedParamNames = c(printedParamNames, "etc.")
  }
  
  labels = c("For standard parameters, posterior chains and density are plotted. ", 
             ">> Population mean and variance: parName.mu, parName.var",
             ">> Effect of task condition: parName_cond[cond]",
             ">> Participant level parameters: parName[pnum]",
             paste0("Standard parameters: ", paste(printedParamNames, collapse=", "))
             )
  
  if (parRes$chains[[1]]$config$modelVariant != "ZL") {
    labels = c(labels, "Category parameters (catMu, catActive) follow with a different format.")
  }
  
  graphics::text(x=0.05, y=seq(0.7, 0.1, length.out=length(labels)), labels=labels, cex=1.2, adj=0)
  graphics::par(xpd=FALSE, mar=opar$mar) # reset
  
  # standard params
  for (i in 1:length(standardParams)) {
    
    titlePagePlot(standardParams[i])
    
    # TODO: Select participants
    matchingParam = getOrderedMatchingParameterNames(parRes$chains[[1]], standardParams[i])
    
    #graphics::par(mfrow=c(1,2))
    for (i in 1:length(matchingParam)) {
      convergencePlotStandard_parallel(parRes, matchingParam[i], type=c("chain", "density"), 
                                       whichIterations=whichIterations, combinedDens=FALSE, panelize=TRUE)
    }
    
  }
  
  # Category parameters
  if (parRes$chains[[1]]$config$modelVariant != "ZL") {
    
    # category params are done for each results object separately so participants have their chains together
    
    catMuScatter_colMaxIter = 50
    
    # title page
    
    explanationText = c("For each participant there are 2 all-chain plots plus a pair of plots for each chain.", 
               "All chains (wide plots): Density for active catMu; sum of catActive per iteration.", 
               "Each chain (pair of plots): catMu scatterplot; catActive sum",
               "catMu scatterplot colors indicate how many iterations the category has been active for.",
               paste0(">> Red = 1 iteration, blue = ", catMuScatter_colMaxIter, " or more iterations."))
  
    if (length(parRes$chains) == 1) {
      explanationText = explanationText[4:5]
    }

    opar = titlePagePlot("Category Parameters", y=0.85, resetPar=FALSE)

    ys = 0.65 - 0.12 * (1:length(explanationText) - 1)
    graphics::text(x=0, y=ys, labels=explanationText, cex=1, adj=0)

    graphics::par(xpd=FALSE, mar=opar$mar) # reset
    
    
    # For each participant, plot the multi chain results, then each individual chain
    for (i in 1:length(usedPnums)) {
      
      if (length(parRes$chains) > 1) {
        #convergencePlotCategory_parallel(resList, pnums=usedPnums[i], panelize=TRUE)
        graphics::par(mfrow=c(1,1))
        convergencePlot_catMuDensity_parallel(parRes, pnums=usedPnums[i], whichIterations=whichIterations)
        convergencePlot_catActive_parallel(parRes, pnums=usedPnums[i], whichIterations=whichIterations)
      }
      
      for (resInd in 1:length(parRes$chains)) {
        #convergencePlotCategory(resList[[resInd]], pnums=usedPnums[i], panelize=TRUE)
        
        graphics::par(mfrow=c(1,2))
        #convergencePlot_catMuScatter(parRes$chains[[resInd]], pnum=usedPnums[i], whichIterations=whichIterations)
        convergencePlot_catMuScatter_activeCumSum(parRes$chains[[resInd]], pnum=usedPnums[i], 
                                                  whichIterations=whichIterations, colMaxIter = catMuScatter_colMaxIter)
        
        fakeParRes = list(chains = list(parRes$chains[[resInd]]), parallelConfig = makeParallelConfig(1))
        convergencePlot_catActive_parallel(fakeParRes, pnums=usedPnums[i], whichIterations=whichIterations)
        
      }
    }
  }
    
  if (!is.null(pdfFile)) {
    grDevices::dev.off()
  }
  
}

convergencePlotStandard = function(res, parName, type=c("chain", "density"), whichIterations=NULL, panelize=TRUE) {
  parRes = list(chains = list(res), parallelConfig = makeParallelConfig(1))
  convergencePlotStandard_parallel(parRes, parName=parName, type=type, whichIterations=whichIterations, combinedDens=FALSE, panelize=panelize)
}


# For non-category parameters (all but catMu and catActive)
# type %in% chain, density
convergencePlotStandard_parallel = function(parRes, parName, type=c("chain", "density"), whichIterations=NULL, combinedDens=FALSE, panelize=TRUE) {
  
  # check for invalid parameters
  if (!(parName %in% names(parRes$chains[[1]]$posteriors)) || 
      grepl("catMu", parName, fixed=TRUE) || 
      grepl("catActive", parName, fixed=TRUE)) 
  {
    warning(paste0("Invalid parameter name: \"", parName, "\""))
    return()
  }
  
  # shared setup
  mm = matrix(nrow=parRes$chains[[1]]$runConfig$iterations, ncol=length(parRes$chains))
  for (i in 1:length(parRes$chains)) {
    mm[,i] = parRes$chains[[i]]$posteriors[[parName]]
  }
  
  whichIterations = valueIfNull(whichIterations, 1:parRes$chains[[1]]$runConfig$iterations)
  whichIterations = whichIterations[ whichIterations %in% 1:nrow(mm) ]
  mm = mm[whichIterations,,drop=FALSE]
  
  valRange = range(mm)
  
  col = grDevices::rainbow(ncol(mm))
  
  
  # Clean type
  type = unique(type[ type %in% c("chain", "density") ])
  if (panelize && length(type) == 2) {
    graphics::par(mfrow=c(1,2))
  }
  
  # chain
  if ("chain" %in% type) {
    base::plot(x=range(whichIterations), y=valRange, type='n', xlab="Iteration", ylab=parName)
    for (i in 1:ncol(mm)) {
      graphics::lines(x=whichIterations, y=mm[,i], col=col[i])
    }
    graphics::title(paste0("Chain for ", parName))
  }
  
  
  # density
  paramIsConstant = valRange[1] == valRange[2]
  if ("density" %in% type) {
    requireNamespace("polspline")
    
    polsplineError = ""
    errorCatcher = function(e) {
      polsplineError <<- paste0("Error from polspline: ", e$message)
    }
    
    # Don't try to do density for constant parameters
    if (paramIsConstant) {
      
      densXs = valRange[1] + c(-1, -0.01, 0, 0.01, 1)
      densMat = matrix(0, nrow=length(densXs), ncol=ncol(mm))
      densMat[3,] = 1
      
    } else {
      densXs = seq(valRange[1], valRange[2], length.out = 100)
      
      densMat = matrix(nrow=length(densXs), ncol=ncol(mm))
      for (i in 1:ncol(mm)) {
        
        tryCatch({
          ls = polspline::logspline(mm[,i])
          densMat[,i] = polspline::dlogspline(densXs, ls)
        }, error = errorCatcher)
        
        if (polsplineError != "") {
          errToPrint = paste0("Error estimating density for parameter = ", parName, ". Density set to 0. ", polsplineError)
          warning(errToPrint)
          polsplineError = ""
          densMat[,i] = 0
        }
        
      }
    }
    
    base::plot(x=valRange, y=c(0,max(densMat)), type='n', xlab=parName, ylab="Density")
    for (i in 1:ncol(mm)) {
      graphics::lines(x=densXs, y=densMat[,i], col=col[i])
    }
    graphics::title(paste0("Density for ", parName))
    
    # Add combined density line
    if (combinedDens && !paramIsConstant) {
      
      tryCatch({
        totalDensLS = polspline::logspline(as.vector(mm))
        totalDensY = polspline::dlogspline(densXs, totalDensLS)
        graphics::lines(x=densXs, y=totalDensY, col="black", lwd=2, lty=2)
      }, error = errorCatcher)
      
      if (polsplineError != "") {
        errToPrint = paste0("Error estimating combined density for parameter=", parName, ". ", polsplineError)
        warning(errToPrint)
        polsplineError = ""
      }
      
    }
  }
}




# single chain, single participant
# Colors are per category index
convergencePlot_catMuScatter = function(res, pnum, whichIterations=NULL) {
  
  whichIterations = valueIfNull(whichIterations, 1:res$runConfig$iterations)
  
  catCol = grDevices::rainbow(res$config$maxCategories)
  alphaHex = gettextf("%X", 51)
  catCol = paste0(catCol, alphaHex)
  
  pMats = convertPosteriorsToMatrices(res, parNames=c("catMu", "catActive"))
  
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
    
    cmij = clampAngle(cmij, pm180 = FALSE, degrees = TRUE)
    
    cmij[ caij == 0 ] = -100 # This is such a gross hack
    
    graphics::points(x=whichIterations, y=cmij, pch=18, col=catCol[j]) # pch=20
  }
  
}

convergencePlot_catMuScatter_activeCumSum = function(res, pnum, whichIterations=NULL, colMaxIter=50) {
  
  if (colMaxIter <= 1) {
    colMaxIter = 1
    colors = "black"
  } else {
    colors = rev(grDevices::hcl.colors(colMaxIter, palette = "Zissou 1")) # Fall is also ok
    
    # Test the palette
    # plot(1:colMaxIter, y=rep(1, colMaxIter), col=colors, pch=16)
    
    alphaHex = gettextf("%X", 51)
    colors = paste0(colors, alphaHex)
  }

  
  whichIterations = valueIfNull(whichIterations, 1:res$runConfig$iterations)
  
  pMats = convertPosteriorsToMatrices(res, parNames=c("catMu", "catActive"))
  
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
    
    cumActive = rep(0, length(caij))
    cumActive[1] = caij[1]
    for (i in 2:length(caij)) {
      if (caij[i] == 1) {
        cumActive[i] = cumActive[i-1] + 1
      } else {
        cumActive[i] = 0
      }
    }
    # Any cumulative active counts get the color for the highest count with a color.
    cumActive[ cumActive > colMaxIter ] = colMaxIter
    
    cmij = clampAngle(cmij, pm180 = FALSE, degrees = TRUE)
    
    cmij[ caij == 0 ] = ylim[1] - 1000 # TODO: This is such a gross hack
    
    activeCounts = sort(unique(cumActive))
    for (ac in activeCounts) {
      whichCount = which(cumActive == ac)
      graphics::points(x=whichIterations[whichCount], y=cmij[whichCount], pch=18, col=colors[ac]) # pch=20
    }
  }
  
}


# All participants in one density plot.
# Only active categories
convergencePlot_catMuDensity_parallel = function(parRes, pnums=NULL, catMuPrecision=2, whichIterations=NULL) {
  nRes = length(parRes$chains)
  pnums = valueIfNull(pnums, parRes$chains[[1]]$pnums)
  chainColors = grDevices::rainbow(nRes)
  
  whichIterations = valueIfNull(whichIterations, 1:parRes$chains[[1]]$runConfig$iterations)
  
  allPlotData = NULL

  for (i in 1:nRes) {
    post = convertPosteriorsToMatrices(parRes$chains[[i]], parNames = c("catMu", "catActive"))
    
    
    includedPnumInd = which(dimnames(post$catMu)[[1]] %in% pnums)
    post$catMu = post$catMu[includedPnumInd,,whichIterations]
    
    includedPnumInd = which(dimnames(post$catActive)[[1]] %in% pnums)
    post$catActive = post$catActive[includedPnumInd,,whichIterations]
    
    
    plotData = catMu_getPlotData(post$catMu, post$catActive, 
                                 precision = catMuPrecision, 
                                 dataType = parRes$chains[[i]]$config$dataType, 
                                 responseRange = parRes$chains[[i]]$config$responseRange, 
                                 colorGeneratingFunction = parRes$chains[[i]]$config$colorGeneratingFunction)
    
    plotData$chain = i
    allPlotData = rbind(allPlotData, plotData)
    
  }
  
  catMu_plotSetup(allPlotData, parRes$chains[[1]]$config$dataType) #, ylim_1 = colorBarYlim[2])
  
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
convergencePlot_catActive_parallel = function(parRes, pnums=NULL, type="line", whichIterations=NULL) {
  
  nRes = length(parRes$chains)
  #nIterations = resList[[1]]$config$iterations
  pnums = valueIfNull(pnums, parRes$chains[[1]]$pnums)
  chainColors = grDevices::rainbow(nRes)
  
  whichIterations = valueIfNull(whichIterations, 1:parRes$chains[[1]]$runConfig$iterations)
  
  # per chain, count of active categories averaged across participants
  catActivePlotData = matrix(0, nrow=length(whichIterations), ncol=nRes)
  
  for (i in 1:nRes) {
    post = convertPosteriorsToMatrices(parRes$chains[[i]], parNames = "catActive")
    
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



# TODO: Work on these catActive diagnostics.

getCatActiveSequential = function(res) {
  catActive = convertPosteriorsToMatrices(res, "catActive")$catActive
  
  seqActive = catActive
  seqActive[,,] = 0
  
  for (i in 1:dim(catActive)[1]) {
    for (j in 1:dim(catActive)[2]) {
      
      caij = catActive[i,j,]
      
      seqActive[i,j,1] = caij[1]
      for (iter in 2:length(caij)) {
        
        if (caij[iter] == 1) {
          seqActive[i,j,iter] = seqActive[i,j,iter-1] + 1
        } else {
          seqActive[i,j,iter] = 0
        }
        
      }
    }
  }
  
  seqActive
}

getCatActiveAlternation = function(res) {
  
  catActive = convertPosteriorsToMatrices(res, "catActive")$catActive
  
  caAlt = catActive
  caAlt[,,] = 0
  
  nIter = dim(catActive)[3]
  
  for (i in 1:dim(catActive)[1]) {
    for (j in 1:dim(catActive)[2]) {
      
      caij = catActive[i,j,]
      
      caAlt[i,j,2:nIter] = caij[1:(nIter-1)] != caij[2:nIter]
    }
  }
  
  caAlt
}


