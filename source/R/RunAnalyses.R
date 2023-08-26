

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
runAnalysisDemo = function(modelVariant="ZL", outputDir="StandardAnalysesDemo") {
  
  demoData = SIM_sampleTestData(modelVariant)
  
  modCfg = makeModelConfig(demoData, modelVariant)
  
  rval = runStandardAnalyses(demoData, modCfg=modCfg, iterations=500, burnIn=0, outputDir=outputDir, mhOptim=FALSE, parallel = NULL)
  
  invisible(rval)
}




#' Fit Model and do Standard Analyses
#' 
#' Runs parameter estimation and does several basic analyses. 
#' Output files are saved to a directory.
#' There are many more possible analyses than are run by this function.
#' 
#' @details 
#' While the analyses are running, output is logged to the console. 
#' Output is also logged to a file named "Analysis Log.txt" and, if using parallel parameter estimation, to additional per-chain log files.
#' 
#' If `loadExistingResults == TRUE` and `iterations` is greater than the number of completed iterations, 
#' `iterations` is treated as the desired total number of iterations and [`continueSampling`] is called to sample more iterations.
#' Note that when sampling more iterations, `iterations` includes MH optimization iterations.
#' When first sampling, however, `iterations` does not including MH optimization iterations (see `mhOptim` argument).
#' 
#'
#' @param data The data to analyze.
#' @param modCfg A model configuration. See [`makeModelConfig`].
#' @param iterations Number of iterations to sample (after MH optimization is complete).
#' @param parallel Passed to `parallel` argument of [`runParameterEstimation`].
#' @param mhOptim Passed to `mhOptim` argument of [`runParameterEstimation`].
#' @param burnIn Number of burn-in iterations to remove from the beginning of chains. See [`removeBurnIn`].
#' @param outputDir Directory for output files. Relative to working directory.
#' @param loadExistingResults If `TRUE`, raw results are loaded from file and analyses are rerun. See details.
#' @param overwrite `NA`: The user is prompted to overwrite. `TRUE`: Analysis files are automatically overwritten. `FALSE`: No files will be overwritten.
#' @param shortName Short name for the analysis that will be prepended to filenames.
#' 
#' @return A list with results objects and analysis results.
#' 
#' @export
#' 
runStandardAnalyses = function(data, modCfg, iterations,
															 parallel = 3,
															 mhOptim = TRUE,
															 burnIn = 1000,
															 outputDir = "standardAnalyses",
															 loadExistingResults = FALSE,
															 overwrite = NA,
															 shortName = NULL
)
{
	
  # @param whichAnalyses A vector of strings naming the analyses to do. See details.
  #whichAnalyses = c("dataPlots", "convergence", "parameterSummary", "testConditionEffects", "testFactors", "categoricalResponding", "WAIC"),
  
	# TODO: Rename analyses

  # Analysis Stages:
  # 1) Pre-model: data plots, NEW response error analysis?, participant characteristics
  # 2) Parameter Estimation: RPE, burn in
  # 3) Model validation and selection: convergence, categorical responding, WAIC
  # 4) Winning model analysis: parameterSummary, pairwise HT, factorial HT
	
	
	# 
	if (!SA_checkOutputDir(outputDir, overwrite, loadExistingResults)) {
		return()
	}
	
	oldWD = getwd()
	setwd(outputDir)
	

	# Route output to a logfile
	sink("Analysis Log.txt", append=FALSE, type = "output", split=TRUE)
	
	#####
	# Startup message
	logMsg("# CatContModel Standard Analyses", "\n")
	
	if (!is.null(shortName)) {
	  logMsg("Short name: ", shortName)
	}
	logMsg("Output directory: ", getwd())
	logMsg("Analyses started: ", date(), "\n")
	
	logMsg("## General Information\n")
	logMsg("These standard analyses are typically useful as part of analyzing continuous report data. There are many more possible analyses.", "\n")
	logMsg("Each analysis has a name, brief description, files produced by the analysis, functions that were used for the analysis, and messages produced during the analysis.", "\n\n", 
	       "The analysis functions are part of the CatContModel package. The help pages for these functions have information about the analysis that was performed and how to read the output from the analysis.", "\n")
	

	# Run the analyses
	tryValue = try({
		preEA = SA_preEstimationAnalyses(data, modCfg, shortName)
		
		res = SA_run_internal(data, modCfg, iterations, parallel, mhOptim, burnIn, loadExistingResults, shortName)
		
		postEA = SA_postEstimationAnalyses(res$rawResults, res$cleanResults, shortName)
		
		c(res, preEA, postEA)
	}, silent=FALSE)
	
	
	# Cleanup messages
	logMsg("\n\n", "# Analyses Completed", "\n\n", "All analyses completed: ", date())
	
	if ("try-error" %in% class(tryValue)) {
	  logMsg("Analyses completed with error: ", tryValue)
	}
	
	
	sink() # Close output file
	setwd(oldWD) # Reset WD
	
	invisible(tryValue)
}



#############################################

SA_checkOutputDir = function(outputDir, overwrite, loadExistingResults) {
  success = FALSE
  
  if (dir.exists(outputDir)) {
    if (!loadExistingResults) {
      if (is.na(overwrite)) {
        if (interactive()) {
          resp = utils::askYesNo(paste0('Directory "', outputDir, '"  already exists. Overwrite analysis files?'))
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


#############################################


# hlev must include its own space at end
SA_logAnalysisInfo = function(name, brief, files, functions, messages=NULL, hlev="", log=TRUE) {
  
  endl = "\n"
  
  str = paste0(hlev, name, endl, endl, 
               brief, endl, endl)
  
  if (!is.null(files)) {
    str = paste0(str, "Files: ", paste(paste0('"', files, '"'), collapse=", "), endl)
  }
  
  if (!is.null(functions)) {
    str = paste0(str, "Functions: ", paste(functions, collapse=", "), endl)
  }
  
  if (!is.null(messages)) {
    str = paste0(str, endl, "Messages: ", paste(messages, collapse=endl))
  }
  
  logMsg(str, disable=!log)
  
  invisible(str)
}

# Expr should return TRUE or FALSE
SA_runSingleAnalysis = function(expr, name, brief, files, functions, hlev="") {
  
  SA_logAnalysisInfo(name=name, hlev=hlev, brief=brief, files=files, functions=functions, messages="")
  
  exprValue = try(expr, silent=FALSE)
  
  if ("try-error" %in% class(exprValue)) {
    logMsg("Analysis completed with errors.")
  } else {
    logMsg("Analysis completed.")
  }
  #} else if (is.logical(exprSuccess)) {
  #	if (exprSuccess) {
  #		logMsg("Analysis completed.")
  #	} else {
  #		logMsg("Analysis completed with possible errors.")
  #	}
  #}
  
  exprValue
}


SA_chancePerformanceAnalysis = function(data, modCfg, perPnum = TRUE) {
  
  getChanceLevelPerformance = function(nObsPerPnum, modCfg, sampleSize=10000) {
    
    # Calculate a chance level by bootstrapping pure guessing participants using actual observation counts per participant.
    # Assume uniform guessing across response space.
    if (modCfg$dataType == "circular") {
      sampleFun = function(n) {
        study = stats::runif(n, 0, 360)
        response = stats::runif(n, 0, 360)
        circDist(study, response, absDist = TRUE, degrees = TRUE)
      }
    } else if (modCfg$dataType == "linear") {
      sampleFun = function(n) {
        study = stats::runif(n, modCfg$responseRange[1], modCfg$responseRange[2])
        response = stats::runif(n, modCfg$responseRange[1], modCfg$responseRange[2])
        
        abs(response - study)
      }
    }
    
    sampleMeans = rep(0, sampleSize)
    
    for (i in 1:sampleSize) {
      samp = sampleFun(nObsPerPnum)
      sampleMeans[i] = mean(samp)
    }
    
    #hist(sampleMeans)
    
    rval = list(mean = mean(sampleMeans), sd = stats::sd(sampleMeans))
    
    # minus 2 SD or quantile methods produce similar results
    rval$chanceMeanMinus2Sd =  rval$mean - 2 * rval$sd
    rval$chanceQuantile = stats::quantile(sampleMeans, 0.05)
    
    rval
  }
  
  
  obsPerParticipant = stats::aggregate(study ~ pnum, data, length)
  names(obsPerParticipant) = c("pnum", "nObs")
  
  obsPerParticipant$chance_MeanMinus2SD = NA
  obsPerParticipant$chance_5thPercentile = NA
  
  if (perPnum) {
    for (uoc in unique(obsPerParticipant$nObs)) {
      
      whichRows = which( obsPerParticipant$nObs == uoc )
      
      cp = getChanceLevelPerformance(uoc, modCfg)
      
      obsPerParticipant$chance_MeanMinus2SD[ whichRows ] = cp$chanceMeanMinus2Sd
      obsPerParticipant$chance_5thPercentile[ whichRows ] = cp$chanceQuantile
      
    }
  } else {
    meanObsPerParticipant = mean(stats::aggregate(study ~ pnum, data, length)$study)
    
    cp = getChanceLevelPerformance(meanObsPerParticipant, modCfg)
    
    obsPerParticipant$chance_MeanMinus2SD = cp$chanceMeanMinus2Sd
    obsPerParticipant$chance_5thPercentile = cp$chanceQuantile
  }
  

  # Look for participants near chance level
  errAgg = stats::aggregate(absErr ~ pnum, data, mean)
  
  # errAgg is aggregated the same as obsPerParticipant, so columns can be copied directly
  errAgg$nObs = obsPerParticipant$nObs
  errAgg$chance_MeanMinus2SD = obsPerParticipant$chance_MeanMinus2SD
  errAgg$chance_5thPercentile = obsPerParticipant$chance_5thPercentile
  
  errAgg = errAgg[ order(errAgg$absErr, decreasing = TRUE), ]
  
  errAgg
}


#############################################

# Pre estimation analyses do use any modeling results.
# 1: Plot raw data
# 2: Observation counts by participant * condition.
# 3: Factorial analysis of response errors. Plots? BayesFactor ANOVA?


SA_preEstimationAnalyses = function(data, modCfg, shortName = NULL) {
  
  # TODO: This is temporary
  #whichAnalyses = c("dataPlots", "observationCounts", "chancePerformance", "factorialErrorAnalysis")
  
  
  # For some of these analyses, absolute error is needed
  if (modCfg$dataType == "circular") {
    data$absErr = circDist(data$study, data$response, absDist = TRUE, degrees = TRUE)
  } else if (modCfg$dataType == "linear") {
    data$absErr = abs(data$response - data$study)
  }
  
  
  rval = list()

  #########################################################
  logMsg("\n\n", "## Analysis Stage 1: Pre-Model Data Exploration", "\n\n", 
         "The first stage of analysis is to examine the raw data.")

  pfx = "1 "
  if (!is.null(shortName)) {
    pfx = paste0(shortName, " ", pfx)
  }
  hlev = "\n\n### "
  
  ###################################
  # Plot raw data
  #if ("dataPlots" %in% whichAnalyses) {
    tempFN = paste0(pfx, "Data Scatterplots.pdf")
    
    rval$dataPlots = SA_runSingleAnalysis({
      plotData(data, pdfFile=tempFN);
      
      "Data plots sucessfully plotted."
    }, 
    name = "Data Plots", hlev=hlev,
    brief = "Scatterplots and histograms of the raw data.", 
    files = tempFN, 
    functions = "plotData()")
  #}
  
  
  ###################################
  # pnum * cond matrix of observation counts
  #if ("observationCounts" %in% whichAnalyses) {
    
    tempFN = paste0(pfx, "Observation Counts Participant by Condition.csv")

    rval$observationCounts = SA_runSingleAnalysis({
      obsCounts = stats::aggregate(study ~ pnum * cond, data, length)
      
      u_pnum = sort(unique(obsCounts$pnum))
      u_cond = unique(obsCounts$cond)
      
      
      m = matrix(0, nrow=length(u_pnum) + 1, ncol=length(u_cond) + 1)
      dimnames(m) = list(c(u_pnum, "cond sum"), c(u_cond, "pnum sum"))
      
      for (i in 1:(nrow(m) - 1)) {
        for (j in 1:(ncol(m) - 1)) {
          pnum = dimnames(m)[[1]][i]
          cond = dimnames(m)[[2]][j]
          m[pnum,cond] = obsCounts$study[ obsCounts$pnum == pnum & obsCounts$cond == cond ]
        }
      }
      
      m[ , length(u_cond) + 1] = apply(m, 1, sum)
      m[length(u_pnum) + 1,] = apply(m, 2, sum)
      
      utils::write.csv(m, tempFN)
      
      m
    }, 
    name = "Observation Counts", hlev=hlev,
    brief = "Number of observations (study/response pairs) by participant and condition.", 
    files = tempFN, 
    functions = "No documented functions used for this analysis.")
    
  #}
  
  
  ###################################
  # Check for participants at chance level
  #if ("chancePerformance" %in% whichAnalyses) {
    
    tempFN = paste0(pfx, "Chance Performance.csv")

    rval$chancePerformance = SA_runSingleAnalysis({
      chanceP = SA_chancePerformanceAnalysis(data, modCfg, perPnum=TRUE)
      utils::write.csv(chanceP, tempFN)
      chanceP
    }, 
    name = "Chance Performance", hlev=hlev,
    brief = "Outputs average response error by participant. Uses a simulated participant who is guessing at random to estimate chance-level cutoffs that may be optionally used by a researcher to exclude participants with average response error above the cutoff.", 
    files = tempFN, 
    functions = "No documented functions used for this analysis.")
    
  #}
  
  ###################################
  # Factorial line chart.
  #if ("factorialErrorAnalysis" %in% whichAnalyses) {
    
    tempFN = paste0(pfx, "Factorial Response Error Plot.pdf")

    rval$factorialErrorAnalysis = SA_runSingleAnalysis({
      wpFactors = getFactorTypeToName(modCfg$factors)$wp
      
      tempData = data
      
      # Copy factors into data
      for (cond in unique(tempData$cond)) {
        rows = which(tempData$cond == cond)
        for (wpf in wpFactors) {
          tempData[ rows , wpf ] = modCfg$factors[ modCfg$factors$cond == cond, wpf ]
        }
      }
      
      form = stats::formula(paste0("absErr ~ ", paste(wpFactors, collapse=" *")))
      
      # TODO: Maybe a table is more useful?
      #agg = aggregate(form, data, mean)
      
      grDevices::pdf(tempFN, width=6, height=6)
      
      plotDf = LineChart::lineChart(form, tempData, errBarType = "Cred95", replicate = "pnum", ylab="Absolute Response Error")
      
      grDevices::dev.off()
      
      plotDf
    }, 
    name = "Factorial Error Analysis", hlev=hlev,
    brief = "Plots absolute response error by the factors in the design. Error bars are 95% credible intervals", 
    files = tempFN, 
    functions = "No documented functions used for this analysis.")
  #}
  
  rval
}


#############################################

SA_run_internal = function(data, modCfg, iterations, parallel, mhOptim, burnIn, loadExistingResults, shortName) {
  
  
  parallelConfig = checkVariableParallelConfig(parallel) #, "RunProgress_")
  
  
  rawResFilename = "2 Raw Results.RDS"
  cleanResFilename = "2 Clean Results.RDS"
  
  if (!is.null(shortName)) {
    rawResFilename = paste0(shortName, " ", rawResFilename)
    cleanResFilename = paste0(shortName, " ", cleanResFilename)
  }
  
  # Startup message

  SA_logAnalysisInfo("\n\n## Analysis Stage 2: Parameter Estimation", 
                     brief = "Parameter estimation can take a while. While waiting, look at the output files from Analysis Stage 1.", 
                     files = c(rawResFilename, cleanResFilename),
                     functions = c("runParameterEstimation()", "continueSampling()", "cleanResults()"),
                     messages = "")
  

  
  ####
  # Load existing results
  rawResults = NULL
  
  # TODO: On load fail, should you warn and re-run or stop?
  if (loadExistingResults) {
    if (file.exists(rawResFilename)) {
      rawResults = readRDS(rawResFilename)
      if (!resultIsType(rawResults, c("WP", "Parallel"))) {
        logWarning("Existing results object (", rawResFilename, ") was not of a valid type.")
        rawResults = NULL
      } else {
        logMsg("Existing results loaded: ", rawResFilename)
      }
    } else {
      logWarning("Loading existing results (", rawResFilename, ") failed: File not found. Rerunning parameter estimation.")
    }
  }
  
  ####
  # Fit model or continue sampling
  
  saveResults = FALSE
  if (is.null(rawResults)) {
    logMsg("Running parameter estimation.")
    
    rawResults = runParameterEstimation(data, modCfg, iterations=iterations, parallel=parallel, mhOptim=mhOptim)
    saveResults = TRUE
    
  } else if (getCompletedIterations(rawResults) < iterations) {
    
    logMsg("Continuing sampling to a total of ", iterations, " iterations.")
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
  cleanRes = cleanResults(rawResults, burnIn)
  
  logMsg("Saving clean results: ", cleanResFilename)
  saveRDS(cleanRes, cleanResFilename)
  
  logMsg("Most of the following analyses use clean results with ", getCompletedIterations(cleanRes), " total iterations.")
  
  list(rawResults = rawResults, cleanResults = cleanRes)
}



SA_postEstimationAnalyses = function(rawRes, cleanRes, shortName=NULL) {
  
  hlev = "\n\n### "
  
  rval = list()
  
  
  ################################################################
  
  logMsg("\n\n## Analysis Stage 3: Model Validation and Selection", "\n\n", 
         "After parameter estimation, the model fits must be validated and a best model selected (if more than one model is being used).")
  
  pfx = "3 "
  if (!is.null(shortName)) {
    pfx = paste0(shortName, " ", pfx)
  }
  
  ####
  # Convergence
  #if ("convergence" %in% whichAnalyses) {
    
    tempFN = paste0(pfx, "Convergence.pdf")
    
    SA_runSingleAnalysis({
      convergencePlots(rawRes, pdfFile=tempFN);
      
      "Convergence plots sucessfully plotted."
    }, 
    name = "Convergence", hlev=hlev,
    brief = "Plots can help identify problems with model convergence. All iterations are plotted, including burn-in iterations that are removed for follow-up analyses. Warnings for this analysis are typically related to density estimation (often mentioning knots) and can be safely ignored.", 
    files=tempFN, 
    functions="convergencePlots()")
    
  #}
  
  ####
  # Analyze clean results
  
  ###
  #if ("categoricalResponding" %in% whichAnalyses) {
    
    if (cleanRes$config$modelVariant == "ZL") {
      
      SA_logAnalysisInfo("Categorical Responding", hlev=hlev,
                      brief="This analysis is not performed for the ZL model variant as there are no categories in that variant.", 
                      files=NULL, functions=NULL, messages=NULL)
    } else {
      
      tempFN = paste0(pfx, "Categorical Responding.csv")
      
      rval$categoricalResponding = SA_runSingleAnalysis({
        catResp = testCategoricalResponding(cleanRes)
        
        utils::write.csv(catResp, tempFN, row.names = FALSE)
        
        catResp
      }, 
      name = "Categorical Responding", hlev=hlev,
      brief = "Hypothesis tests of whether participants' responses are entirely categorical or entirely continuous.", 
      files = tempFN, 
      functions = "testCategoricalResponding()")
    }
    
  #}
  
  ###
  #if ("WAIC" %in% whichAnalyses) {
    
    tempFN = paste0(pfx, "WAIC.csv")
    
    rval$WAIC = SA_runSingleAnalysis({
      waic = calculateWAIC(cleanRes)
      
      utils::write.csv(waic, tempFN, row.names=FALSE)
      
      waic
    }, 
    name = "WAIC", hlev=hlev,
    brief = "After running multiple models on the same data, the best fitting model can be selected using WAIC. If you only fit one model, WAIC is not informative.", 
    files = tempFN, 
    functions = "calculateWAIC()")
    
  #}
  
  #####################################################################
  
  logMsg("\n\n## Analysis Stage 4: Analysis of Selected Model", "\n\n", 
         "Once a model has been selected, its parameters can be analyzed to answer substantive research equations.")
  
  pfx = "4 "
  if (!is.null(shortName)) {
    pfx = paste0(shortName, " ", pfx)
  }
  
  ###
  #if ("parameterSummary" %in% whichAnalyses) {
    
    
    tempFN_plot = paste0(pfx, "Parameter Summary Plots.pdf")
    tempFN_table = paste0(pfx, "Parameter Summary Table.csv")
    
    
    rval$parameterSummary = SA_runSingleAnalysis({
      
      # Plot summary
      plotParameterSummary(cleanRes, pdfFile=tempFN_plot, openPdf = FALSE)
      
      # Numeric summary
      PMCI = posteriorMeansAndCredibleIntervals(cleanRes)
      
      utils::write.csv(PMCI, tempFN_table, row.names=FALSE)
      
      PMCI
    }, 
    name = "Parameter Summaries", hlev=hlev,
    brief = "Model parameters are summarized in plots and tables.", 
    files = c(tempFN_plot, tempFN_table), 
    functions = c("plotParameterSummary()", "posteriorMeansAndCredibleIntervals()"))

  #}
  
  ###
  #if ("testConditionEffects" %in% whichAnalyses) {
    tempFN = paste0(pfx, "Pairwise Comparisons of Conditions.csv")
    
    rval$conditionEffects = SA_runSingleAnalysis({
      condEff = testConditionEffects(cleanRes)
      utils::write.csv(condEff, tempFN, row.names = FALSE)
      
      condEff
    }, 
    name = "Pairwise Comparisons of Conditions", hlev=hlev,
    brief = "Pairwise differences between task conditions are tested with Bayes factors and credible intervals.", 
    files = tempFN, 
    functions = "testConditionEffects()")
  #}
  
  ###
  #if ("testFactors" %in% whichAnalyses) {
    
    tempFN = paste0(pfx, "Factorial Tests.csv")
    
    rval$factorialTests = SA_runSingleAnalysis({
      mei = testMainEffectsAndInteractions(cleanRes)
      
      utils::write.csv(mei, tempFN, row.names = FALSE)
      
      mei
    }, 
    name = "Hypothesis Tests on Factors", hlev=hlev,
    brief = "ANOVA-like hypothesis tests on factors in the design.", 
    files = tempFN, 
    functions = "testMainEffectsAndInteractions()")

  #}
  
  rval
}



