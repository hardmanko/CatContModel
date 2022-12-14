
runAnalysisDemo = function(modelVariant="betweenItem") {
  
  demoData = SIM_sampleTestData(modelVariant)
  
  modCfg = makeModelConfig(demoData, modelVariant)
  
  rval = basicAnalyses(demoData, modCfg, iterations=500, burnIn=0, outputDir="analysisDemo", mhOptim=FALSE, parallel = NULL)
  
  invisible(rval)
}



#' Fit Model and do Standard Analyses
#' 
#' Runs parameter estimation and does several basic analyses. This function is like a demo. 
#' There are many more possible analyses implemented in the CatContModel package than are run by this function.
#' 
#' @details 
#' If `loadExistingResults == TRUE` and `iterations` is greater than the number of completed iterations, 
#' `iterations` is treated as the desired total number of iterations and [`continueSampling`] is called.
#' 
#' Learn more about the analyses by reading documentation for the analysis functions:
#' + `dataPlots`: [`plotData`].
#' + 
#' + `WAIC`: [`calculateWAIC`].
#' + `testFactors`: [`testMainEffectsAndInteractions`].
#'
#' @param data The data to analyze.
#' @param config A model configuration. See [`makeModelConfig`].
#' @param iterations Number of iterations to sample.
#' @param parallel Passed to `parallel` argument of [`runParameterEstimation`].
#' @param mhOptim Passed to `mhOptim` argument of [`runParameterEstimation`].
#' @param burnIn Number of burn-in iterations to remove from the beginning of chains. See [`removeBurnIn`].
#' @param whichAnalyses A vector of strings naming the analyses to do. See details.
#' @param outputDir Directory for output files. Relative to working directory.
#' @param loadExistingResults If `TRUE`, raw results are loaded from file and analyses are rerun. See details.
#' 
#' @return 
#' 
runStandardAnalyses = function(data, config, iterations,
                         parallel = 3,
                         mhOptim = TRUE,
                         burnIn = 500,
                         whichAnalyses = c("dataPlots", "convergence", "WAIC", "plotParameterSummary", "testConditionEffects", "testFactors", "categoricalResponding"),
                         outputDir = "standardAnalyses",
                         loadExistingResults = FALSE
                         )
{
  
  # TODO: Analysis names. Not just here, but maybe the main functions should be renamed.
  # Should ANOVA be FactorTests? OmnibusTests? testFactors?
  # testConditionEffects -> testConditionPairs?
  
  
  if (dir.exists(outputDir)) {
    # TODO: Stuff about overwriting.
  } else {
    dir.create(outputDir)
    if (!dir.exists(outputDir)) {
      stop(paste0('Could not create the output directory "', getwd(), "/", outputDir, '" Try creating it manually before running analyses.'))
    }
  }
  
  
  
  oldWD = getwd()
  setwd(outputDir)
  
  message(paste0("Output files will be written to: ", getwd()))
  
  ####
  if ("dataPlots" %in% whichAnalyses) {
    tempFN = "dataPlots.pdf"
    message(paste0("\nPlotting raw data: ", tempFN))
    plotData(data, pdfFile=tempFN)
  }
  
  ####
  # Load existing results
  rawResults = NULL
  rawResFilename = "rawResults.RDS"
  
  if (loadExistingResults) {
    if (file.exists(rawResFilename)) {
      rawResults = readRDS(rawResFilename)
      if (!resultIsType(rawResults, c("WP", "Parallel"))) {
        warning(paste0("Existing results object (", rawResFilename, ") was not of a valid type."))
        rawResults = NULL
      } else {
        message(paste0("\nExisting results loaded: ", rawResFilename))
      }
    } else {
      warning(paste0("Loading results object (", rawResFilename, ") failed: File not found. Rerunning fit."))
    }
  }
  
  ####
  # Fit model or continue sampling
  
  saveResults = FALSE
  if (is.null(rawResults)) {
    message("\nFitting the model.")
    rawResults = runParameterEstimation(data, config, iterations=iterations, parallel=parallel, mhOptim=mhOptim)
    saveResults = TRUE
    
  } else if (getCompletedIterations(rawResults) < iterations) {
    
    message(paste0("\nContinuing sampling to a total of ", iterations, " iterations."))
    rawResults = continueSampling(rawResults, iterations, totalIterations = TRUE)
    saveResults = TRUE
  }
  
  if (saveResults) {
    message(paste0("Parameter estimation complete. Saving raw results: ", rawResFilename))
    saveRDS(rawResults, rawResFilename)
  }
  
  
  

  
  ####
  # Convergence
  if ("convergence" %in% whichAnalyses) {
    
    tempFN = "convergence.pdf"
    message(paste0("\nMaking convergence plots: ", tempFN))
    
    convergencePlots(rawResults, filename=tempFN)
  }
  
  ####
  # Remove burn-in iterations and combine parallel chains into one chain
  if (burnIn >= getCompletedIterations(rawResults)) {
    message("burnIn >= iterations. burnIn set to 0.")
    burnIn = 0
  }
  
  message(paste0("Removing ", burnIn, " burn in iterations (from each chain)."))
  cleanRes = removeBurnIn(rawResults, burnIn)
  if (resultIsType(cleanRes, "Parallel")) {
    cleanRes = combineResults(resList = cleanRes$chains)
  }
  
  ###
  cleanResFilename = "cleanResults.RDS"
  message(paste0("Saving clean results: ", cleanResFilename))
  saveRDS(cleanRes, cleanResFilename)
  
  rval = list(
    rawResults = rawResults,
    cleanResults = cleanRes
  )
  
  message(paste0("Following analyses performed on clean results with ", getCompletedIterations(cleanRes), " iterations."))
  
  ####
  # Analyze clean results
  
  if ("categoricalResponding" %in% whichAnalyses && cleanRes$config$modelVariant != "ZL") {
    
    tempFN = "testCategoricalResponding.csv"
    
    message(paste0("\nTesting for evidence of categorical responding: ", tempFN))
    
    catResp = testCategoricalResponding(cleanRes)
    
    utils::write.csv(catResp, tempFN, row.names = FALSE)
    
    rval$categoricalResponding = catResp
  }
  
  if ("plotParameterSummary" %in% whichAnalyses) {
    tempFN = "parameterSummary.pdf"
    message(paste0("\nPlotting parameter summary: ", tempFN))
    
    plotParameterSummary(cleanRes, pdfFile=tempFN, openPdf = FALSE)
    
    # TODO: Also do numerical summary?
    #PMCI = posteriorMeansAndCredibleIntervals(cleanRes)
  }
  
  ###
  if ("testConditionEffects" %in% whichAnalyses) {
    tempFN = "testConditionEffects.csv"
    
    message(paste0("\nTesting pairwise differences between task conditions: ", tempFN))
    
    condEff = testConditionEffects(cleanRes)
    utils::write.csv(condEff, tempFN, row.names = FALSE)
    
    rval$conditionEffects = condEff
  }
  
  ###
  if ("testFactors" %in% whichAnalyses) {
    
    fn = "testFactors.csv"
    
    message(paste0("\nANOVA-like hypothesis tests on factors: ", fn))
    
    mei = testMainEffectsAndInteractions(cleanRes)
    
    utils::write.csv(mei, fn, row.names = FALSE)
    
    rval$MEI = mei
  }
  
  if ("WAIC" %in% whichAnalyses) {
    tempFN = "WAIC.csv"
    message(paste0("\nCalculating WAIC: ", tempFN))
    
    waic = calculateWAIC(cleanRes)
    utils::write.csv(waic, tempFN, row.names=FALSE)
    
    rval$WAIC = waic
  }
  
  setwd(oldWD)
  
  invisible(rval)
}



# The idea is to split basicAnalyses into fitting and analyzing, but ordering is a pain.
# This function may not work right because of ordering, just use basic basicAnalyses
standardAnalyses = function(res, burnIn=0, 
                            whichAnalyses = c("dataPlots", "convergence", "WAIC", "plotParameterSummary", "testConditionEffects", "ANOVA")) 
{
  if (burnIn > 0) {
    res = removeBurnIn(res, burnIn)
  }
  
  # Collapse parallel chains
  if (resultIsType(res, "Parallel")) {
    res = combineResults(resList = res$results)
  }
  
  # TODO: BP?
  
  
  
}


exploratoryAnalysis = function(data, modelVariants, outputDir, iterations=2500, burnIn=500, mhOptim=TRUE, parallel=1) {
  

  for (mv in modelVariants) {
    config = makeModelConfig(data, modelVariant=mv, iterations=iterations)
    
    if (length(modelVariants) > 1) {
      outputDir = paste0(outputDir, "/", mv)
    }
    
    basicAnalyses(data = data, config = config, outputDir=outputDir, 
                  iterations = iterations, burnIn=burnIn, 
                  parallel = parallel, mhOptim = mhOptim)
  }
  
  
}


