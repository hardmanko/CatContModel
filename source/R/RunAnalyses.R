



#' Package Demo
#' 
#' Samples simulated data, fits the model, does standard analyses, and writes output. 
#' See [`runStandardAnalyses`] if you have data to fit but otherwise want a demo.
#' 
#' @param modelVariant Name of a model variant. See `modelVariant` argument of [`makeModelConfig`] for names.
#' @param outputDir Where to save output (relative to working directory).
#' 
#' @export
#' @return Invisibly, the same value as [`runStandardAnalyses`].
runAnalysisDemo = function(modelVariant="betweenItem", outputDir="analysisDemo") {
  
  demoData = SIM_sampleTestData(modelVariant)
  
  modCfg = makeModelConfig(demoData, modelVariant)
  
  rval = runStandardAnalyses(demoData, config=modCfg, iterations=500, burnIn=0, outputDir=outputDir, mhOptim=FALSE, parallel = NULL)
  
  invisible(rval)
}



SA_checkOutputDir = function(outputDir, overwrite) {
	success = FALSE
	
	if (dir.exists(outputDir)) {
		if (!loadExistingResults) {
			if (is.na(overwrite)) {
				if (interactive()) {
					resp = askYesNo(paste0('Directory "', outputDir, '"  already exists. Overwrite analysis files?'))
					overwrite = is.logical(resp) && resp == TRUE
				} else {
					overwrite = FALSE
				}
			}
			
			if (overwrite) {
				logMsg('Analysis files in "', outputDir, '" will be overwritten.')
				success = TRUE
			} else {
				logMsg('Directory "', outputDir, '" will not be overwritten.')
				# return()
				success = FALSE
			}
		}
	} else {
		dir.create(outputDir, recursive = TRUE)
		if (!dir.exists(outputDir)) {
			stop(paste0('Could not create the output directory "', outputDir, '" Try creating it manually before running analyses.'))
		}
		success = TRUE
	}
	
	success
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
#' + `dataPlots`: Plots of raw data. See [`plotData`].
#' + `convergence`: Plots of model parameter convergence. See [`convergencePlots`].
#' + `parameterSummary`: Means and credible intervals for model parameters. For plots, see [`plotParameterSummary`]. For numerical summary, see [`posteriorMeansAndCredibleIntervals`].
#' + `testConditionEffects`: Comparisons between conditions (`cond`). See [`testConditionEffects`].
#' + `testFactors`: ANOVA-like tests of factors. See [`testMainEffectsAndInteractions`] for tests, [`makeFactors`] to map conditions to factor levels, and `factors` argument of [`makeModelConfig`] to use your factors.
#' + `categoricalResponding`: Tests for presence of categorical responding (for models with categorical responding). See [`testCategoricalResponding`].
#' + `WAIC`: Calculates model with using WAIC. See [`calculateWAIC`]. Only useful if you fit more than one model variant to the data.
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
#' @param overwrite `NA`: The user is prompted to overwrite. `TRUE`: Analysis files are automatically overwritten. `FALSE`: No files will be overwritten.
#'
#' 
#' @return A list with results objects and analysis results.
#' 
#' @export
#' 
runStandardAnalyses = function(data, config, iterations,
                         parallel = 3,
                         mhOptim = TRUE,
                         burnIn = 500,
                         whichAnalyses = c("dataPlots", "convergence", "parameterSummary", "testConditionEffects", "testFactors", "categoricalResponding", "WAIC"),
                         outputDir = "standardAnalyses",
                         loadExistingResults = FALSE,
												overwrite = NA
                         )
{

	# TODO: Rename analyses
	
	# @param analysisDescription Optional description that will be logged.
	
	if (!SA_checkOutputDir(outputDir, overwrite)) {
		return()
	}
  
	# How to handle stops within called functions? TryCatch with sink outside?
	# sink("Analysis Log.txt", append=FALSE, type = "output", split=TRUE)
  
  oldWD = getwd()
  setwd(outputDir)
  
  logMsg("Output files will be written to: ", getwd())
  
  ####
  # TODO: Should data plots be before running chains? Probably? But so much easier if run after.
  if ("dataPlots" %in% whichAnalyses) {
    tempFN = "dataPlots.pdf"
    logMsg("\nPlotting raw data: ", tempFN)
    plotData(data, pdfFile=tempFN)
  }
  
  ####
  # Load existing results
  rawResults = NULL
  rawResFilename = "rawResults.RDS"
  
  # TODO: On load fail, should you warn and re-run or stop?
  if (loadExistingResults) {
    if (file.exists(rawResFilename)) {
      rawResults = readRDS(rawResFilename)
      if (!resultIsType(rawResults, c("WP", "Parallel"))) {
        #warning(paste0("Existing results object (", rawResFilename, ") was not of a valid type."))
      	logWarning("Existing results object (", rawResFilename, ") was not of a valid type.")
        rawResults = NULL
      } else {
        logMsg("\nExisting results loaded: ", rawResFilename)
      }
    } else {
      #warning(paste0("Loading existing results (", rawResFilename, ") failed: File not found. Rerunning parameter estimation."))
    	logWarning("Loading existing results (", rawResFilename, ") failed: File not found. Rerunning parameter estimation.")
    }
  }
  
  ####
  # Fit model or continue sampling
  
  saveResults = FALSE
  if (is.null(rawResults)) {
    logMsg("\nFitting the model.")
  	
    rawResults = runParameterEstimation(data, config, iterations=iterations, parallel=parallel, mhOptim=mhOptim)
    saveResults = TRUE
    
  } else if (getCompletedIterations(rawResults) < iterations) {
    
  	logMsg("\nContinuing sampling to a total of ", iterations, " iterations.")
    rawResults = continueSampling(rawResults, iterations, totalIterations = TRUE)
    saveResults = TRUE
  }
  
  if (saveResults) {
  	logMsg("Parameter estimation complete. Saving raw results: ", rawResFilename)
    saveRDS(rawResults, rawResFilename)
  }
  
  
  ####
  # Convergence
  if ("convergence" %in% whichAnalyses) {
    
    tempFN = "convergence.pdf"
    logMsg("\nMaking convergence plots: ", tempFN)
    
    convergencePlots(rawResults, pdfFile=tempFN)
  }
  
  ####
  # Remove burn-in iterations and combine parallel chains into one chain
  if (burnIn >= getCompletedIterations(rawResults)) {
  	logWarning("burnIn >= iterations. burnIn set to 0.")
    burnIn = 0
  }
  
  logMsg("Removing ", burnIn, " burn in iterations (from each chain).")
  cleanRes = removeBurnIn(rawResults, burnIn)
  if (resultIsType(cleanRes, "Parallel")) {
    cleanRes = combineResults(resList = cleanRes$chains)
  }
  
  ###
  cleanResFilename = "cleanResults.RDS"
  logMsg("Saving clean results: ", cleanResFilename)
  saveRDS(cleanRes, cleanResFilename)
  
  rval = list(
    rawResults = rawResults,
    cleanResults = cleanRes
  )
  
  logMsg("Following analyses performed on clean results with ", getCompletedIterations(cleanRes), " iterations.")
  
  ####
  # Analyze clean results
  
  if ("categoricalResponding" %in% whichAnalyses && cleanRes$config$modelVariant != "ZL") {
    
    tempFN = "testCategoricalResponding.csv"
    
    logMsg("\nTesting for evidence of categorical responding: ", tempFN)
    
    catResp = testCategoricalResponding(cleanRes)
    
    utils::write.csv(catResp, tempFN, row.names = FALSE)
    
    rval$categoricalResponding = catResp
  }
  
  if ("parameterSummary" %in% whichAnalyses) {
    
    # Plot summary
    tempFN = "parameterSummary.pdf"
    logMsg("\nPlotting parameter summary: ", tempFN)
    
    plotParameterSummary(cleanRes, pdfFile=tempFN, openPdf = FALSE)
    
    # Numerical summary
    tempFN = "parameterSummary.csv"
    logMsg("\nWriting numerical parameter summary: ", tempFN)
    PMCI = posteriorMeansAndCredibleIntervals(cleanRes)
    
    utils::write.csv(PMCI, tempFN, row.names=FALSE)
  }
  
  ###
  if ("testConditionEffects" %in% whichAnalyses) {
    tempFN = "testConditionEffects.csv"
    
    logMsg("\nTesting pairwise differences between task conditions: ", tempFN)
    
    condEff = testConditionEffects(cleanRes)
    utils::write.csv(condEff, tempFN, row.names = FALSE)
    
    rval$conditionEffects = condEff
  }
  
  ###
  if ("testFactors" %in% whichAnalyses) {
    
    fn = "testFactors.csv"
    
    logMsg("\nANOVA-like hypothesis tests on factors: ", fn)
    
    mei = testMainEffectsAndInteractions(cleanRes)
    
    utils::write.csv(mei, fn, row.names = FALSE)
    
    rval$MEI = mei
  }
  
  if ("WAIC" %in% whichAnalyses) {
    tempFN = "WAIC.csv"
    logMsg("\nCalculating WAIC: ", tempFN)
    
    waic = calculateWAIC(cleanRes)
    utils::write.csv(waic, tempFN, row.names=FALSE)
    
    rval$WAIC = waic
  }
  
  setwd(oldWD)
  
  invisible(rval)
}

SA_preEstimationAnalyses = function(whichAnalyses, data) {
	
	if ("dataPlots" %in% whichAnalyses) {
		tempFN = "dataPlots.pdf"
		logMsg("\nPlotting raw data: ", tempFN)
		plotData(data, pdfFile=tempFN)
	}
	
}

SA_postEstimationAnalyses = function(whichAnalyses, rawRes, cleanRes) {
	
	rval = list()
	
	####
	# Convergence
	if ("convergence" %in% whichAnalyses) {
		
		tempFN = "convergence.pdf"
		logMsg("\nMaking convergence plots: ", tempFN)
		
		convergencePlots(rawRes, pdfFile=tempFN)
	}
	
	logMsg("Following analyses performed on clean results with ", getCompletedIterations(cleanRes), " iterations.")
	
	####
	# Analyze clean results
	
	if ("categoricalResponding" %in% whichAnalyses && cleanRes$config$modelVariant != "ZL") {
		
		tempFN = "testCategoricalResponding.csv"
		
		logMsg("\nTesting for evidence of categorical responding: ", tempFN)
		
		catResp = testCategoricalResponding(cleanRes)
		
		utils::write.csv(catResp, tempFN, row.names = FALSE)
		
		rval$categoricalResponding = catResp
	}
	
	if ("parameterSummary" %in% whichAnalyses) {
		
		# Plot summary
		tempFN = "parameterSummary.pdf"
		logMsg("\nPlotting parameter summary: ", tempFN)
		
		plotParameterSummary(cleanRes, pdfFile=tempFN, openPdf = FALSE)
		
		# Numerical summary
		tempFN = "parameterSummary.csv"
		logMsg("\nWriting numerical parameter summary: ", tempFN)
		PMCI = posteriorMeansAndCredibleIntervals(cleanRes)
		
		utils::write.csv(PMCI, tempFN, row.names=FALSE)
	}
	
	###
	if ("testConditionEffects" %in% whichAnalyses) {
		tempFN = "testConditionEffects.csv"
		
		logMsg("\nTesting pairwise differences between task conditions: ", tempFN)
		
		condEff = testConditionEffects(cleanRes)
		utils::write.csv(condEff, tempFN, row.names = FALSE)
		
		rval$conditionEffects = condEff
	}
	
	###
	if ("testFactors" %in% whichAnalyses) {
		
		fn = "testFactors.csv"
		
		logMsg("\nANOVA-like hypothesis tests on factors: ", fn)
		
		mei = testMainEffectsAndInteractions(cleanRes)
		
		utils::write.csv(mei, fn, row.names = FALSE)
		
		rval$MEI = mei
	}
	
	if ("WAIC" %in% whichAnalyses) {
		tempFN = "WAIC.csv"
		logMsg("\nCalculating WAIC: ", tempFN)
		
		waic = calculateWAIC(cleanRes)
		utils::write.csv(waic, tempFN, row.names=FALSE)
		
		rval$WAIC = waic
	}
	
	rval
}

SA_run_internal = function(data, config, iterations, parallel, mhOptim, burnIn, loadExistingResults) {
	
	####
	# Load existing results
	rawResults = NULL
	rawResFilename = "rawResults.RDS"
	
	# TODO: On load fail, should you warn and re-run or stop?
	if (loadExistingResults) {
		if (file.exists(rawResFilename)) {
			rawResults = readRDS(rawResFilename)
			if (!resultIsType(rawResults, c("WP", "Parallel"))) {
				logWarning("Existing results object (", rawResFilename, ") was not of a valid type.")
				rawResults = NULL
			} else {
				logMsg("\nExisting results loaded: ", rawResFilename)
			}
		} else {
			logWarning("Loading existing results (", rawResFilename, ") failed: File not found. Rerunning parameter estimation.")
		}
	}
	
	####
	# Fit model or continue sampling
	
	saveResults = FALSE
	if (is.null(rawResults)) {
		logMsg("\nRunning parameter estimation.")
		
		rawResults = runParameterEstimation(data, config, iterations=iterations, parallel=parallel, mhOptim=mhOptim)
		saveResults = TRUE
		
	} else if (getCompletedIterations(rawResults) < iterations) {
		
		logMsg("\nContinuing sampling to a total of ", iterations, " iterations.")
		rawResults = continueSampling(rawResults, iterations, totalIterations = TRUE)
		saveResults = TRUE
	}
	
	if (saveResults) {
		logMsg("Parameter estimation complete. Saving raw results: ", rawResFilename)
		saveRDS(rawResults, rawResFilename)
	}
	
	
	####
	# Remove burn-in iterations and combine parallel chains into one chain
	if (burnIn >= getCompletedIterations(rawResults)) {
		logWarning("burnIn >= iterations. burnIn set to 0.")
		burnIn = 0
	}
	
	logMsg("Removing ", burnIn, " burn in iterations (from each chain).")
	cleanRes = removeBurnIn(rawResults, burnIn)
	if (resultIsType(cleanRes, "Parallel")) {
		cleanRes = combineResults(resList = cleanRes$chains)
	}
	
	###
	cleanResFilename = "cleanResults.RDS"
	logMsg("Saving clean results: ", cleanResFilename)
	saveRDS(cleanRes, cleanResFilename)
	

	list(rawResults = rawResults, cleanResults = cleanRes)
}



runStandardAnalyses_2 = function(data, config, iterations,
															 parallel = 3,
															 mhOptim = TRUE,
															 burnIn = 500,
															 whichAnalyses = c("dataPlots", "convergence", "parameterSummary", "testConditionEffects", "testFactors", "categoricalResponding", "WAIC"),
															 outputDir = "standardAnalyses",
															 loadExistingResults = FALSE,
															 overwrite = NA
)
{
	
	
	# 
	if (!SA_checkOutputDir(outputDir, overwrite)) {
		return()
	}
	
	oldWD = getwd()
	setwd(outputDir)
	
	logMsg("Output files will be written to: ", getwd())
	
	# How to handle stops within called functions? TryCatch with sink outside?
	sink("Analysis Log.txt", append=FALSE, type = "output", split=TRUE)
	
	trySuccess = FALSE
	
	tryValue = try({
		SA_preEstimationAnalyses(whichAnalyses, data)
		
		rawRes = SA_run_internal(data, config, iterations, parallel, mhOptim, burnIn, loadExistingResults)
		
		SA_postEstimationAnalyses(whichAnalyses, rawRes, cleanRes)
		
		trySuccess = TRUE
	}, silent=FALSE)
	
	if (!trySuccess) {
		# print the error
		logWarning(tryValue) # ???
	}

	# TODO: Close all open files? How?
	# Reset WD
	
}




