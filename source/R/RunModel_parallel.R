

#' Configure Parallel Parameter Estimation
#' 
#' Creates a configuration list that can be passed as the `parallelConfig` argument of [`runParameterEstimation`].
#' 
#' @param numNodes The number of cluster nodes to use (on a typical computer, each node is 1 CPU core). Each node will run one chain at a time.
#' @param numChains The number of parameter estimation chains to run in parallel. Can be greater than `numNodes`.
#' @param outputFile A text file name (e.g. "progress.txt") or an incomplete file name without .txt extension (e.g. "progress_"). See details.
#' @param cluster A cluster created with `parallel::makeCluster()`. If `NULL`, when parameter estimation is run, a cluster will be created with `numNodes` nodes.
#' 
#' @details 
#' The parallel parameter estimation functions use the R package `parallel` to run multiple 
#' parallel parameter chains, not to speed estimation of one parameter chain.
#' Once the parallel chains have been sampled, 
#' 1) each chain must have burn-in iterations removed (see [`removeBurnIn`]) and 
#' 2) the chains must be combined into a single results object (see [`combineResults`]).
#' 
#' There are two major advantages to parallel estimation: Speed and convergence diagnostics.
#' Speed: By running chains in parallel, more iterations of the parameter chains can be sampled more quickly.
#' Convergence: If multiple chains, each with different starting values, all converge to the same
#' distributions, that is a strong convergence diagnostic. (See [`convergencePlots`], which works with parallel chains.)
#' 
#' If `cluster` is not provided, a cluster will be created on the local machine using `numNodes` parallel R processes.
#' For best performance (shortest runtime), do not exceed the number of physical CPU cores on your machine. 
#' However, you may want more chains than you have cores, in which case you simply request the number of desired chains and wait longer.
#' Note: use `parallel::detectCores(logical=FALSE)` to detect number of physical cores.
#' 
#' Argument `outputFile` specifies how nodes should output their progress. 
#' Progress is normally printed to the console, but that doesn't work with parallel nodes.
#' 1) If `outputFile` is a text file with extension ("ParallelProgress.txt"), output from only the first chain is printed to that file. Not recommended if `numChains > numNodes`.
#' 2) If `outputFile` is an incomplete file name ("ParallelProgress_"), one file will be created for each chain (e.g. "ParallelProgress_1.txt").
#' 3) If `outputFile` is the empty string ("") or `NULL`, no output file is created.
#' File contents are overwritten each time a new run starts.
#' Files are made for each chain, not each node.
#' 
#' Each parallel chain uses a RNG seed that is sampled in the main process. This means that using `set.seed` before
#' starting parallel parameter estimation will correctly produce repeatable chains.
#' 
#' For between-participants designs, it would be helpful to be able to run different participant groups in parallel,
#' but the parallel parameter estimation functions do not help with that. You will still need to run estimation for
#' each group separately. You can use the parallel functions to run parallel chains for each group, however, which can
#' speed parameter estimation.
#' 
#' @return A configuration list that can be passed as the `parallelConfig` argument of [`runParameterEstimation`].
#' 
#' @export
makeParallelConfig = function(numNodes,
                              numChains = numNodes,
                              outputFile = "ParallelProgress_",
                              cluster = NULL) 
{

  if (numNodes <= 0) {
    stop("Must have at least 1 node.")
  }

  parallelConfig = list(numNodes = numNodes,
                        numChains = numChains,
                        outputFile = outputFile,
                        cluster = cluster
  )
  
  parallelConfig
}

# Handles single number or list

checkVariableParallelConfig = function(parallel, outputFile) {
	
	useOutputFile = !missing(outputFile)
	
	parallelConfig = NULL
	
	if (!is.null(parallel)) {
		if (is.list(parallel)) {
			parallelConfig = parallel
			if (useOutputFile) {
				parallelConfig$outputFile = outputFile
			}
		} else if (is.numeric(parallel) && length(parallel) == 1 && parallel > 1) {
			
			if (useOutputFile) {
				parallelConfig = makeParallelConfig(numNodes = parallel, outputFile = outputFile)
			} else {
				parallelConfig = makeParallelConfig(numNodes = parallel)
			}
		}
	}
	
	parallelConfig
}



runParallelJobs = function(pc, jobs, jobFun=NULL, startupMsg=NULL) {
  
  requireNamespace("parallel")
  
  for (i in 1:length(jobs)) {
    # Sample seeds for the jobs
    if (is.null(jobs[[i]]$seed)) {
      jobs[[i]]$seed = sample.int(.Machine$integer.max, 1)
    }
    
    # If jobFun is provided, set it per job
    if (!is.null(jobFun)) {
      jobs[[i]]$jobFun = jobFun
    }
    
    if (is.null(jobs[[i]]$jobFun)) {
      stop("No job function provided.")
    }
  }
  
  # Make the cluster if needed
  createCluster = is.null(pc$cluster)
  if (createCluster) {
    pc$cluster = parallel::makeCluster(pc$numNodes)
  }
  
  # Startup message
  if (is.null(startupMsg) || !is.character(startupMsg)) {
    logMsg("Starting ", length(jobs), " parallel jobs on ", pc$numNodes, " nodes.")
  } else {
  	logMsg(startupMsg)
  }
  
  
  # Deal with output filenames
  if (!is.null(pc$outputFile) && pc$outputFile != "") {
    
    outputFiles = rep("", length(jobs))
    
    if (endsWith(pc$outputFile, ".txt")) {
      outputFiles[1] = pc$outputFile
      
      logMsg("Parallel progress for the first job will be written to: ", pc$outputFile)
    } else {

      for (i in 1:length(jobs)) {
        outputFiles[i] = paste0(pc$outputFile, i, ".txt")
      }
      
    	logMsg("Parallel progress for each job will be written to: ", paste(outputFiles, collapse=", "))
    }
    
    for (i in 1:length(jobs)) {
      jobs[[i]]$outputFile = outputFiles[i]
    }
    
  } else {
  	logMsg("No parallel progress file will be written. See outputFile argument of makeParallelConfig().")
  }
  
 
  jobRunnerWapper = function(job) {
    require(CatContModel)
    
    set.seed(job$seed) # if NULL, seeds randomly
    
    if (!is.null(job$outputFile) && job$outputFile != "") {
      sink(job$outputFile, append=FALSE)
    }
    
    rval = job$jobFun(job)
    
    rval
  }
  
  # Run the jobs
  rval = parallel::clusterApply(pc$cluster, jobs, jobRunnerWapper)
  
  # Done: Stop the cluster
  if (createCluster) {
    parallel::stopCluster(pc$cluster)
    pc$cluster = NULL
  }
  
  rval
}


# This assumes that the whole job is passed to jobFun (including jobFun)
makeParallelJob = function(seed=NULL, outputFile=NULL, jobFun=NULL, ...) {
  
  job = list(...)
  job$jobFun = jobFun
  job$seed = seed
  job$outputFile = outputFile
  
  job
  
}



RPE_parallel = function(data, modCfg, runConfig, parallelConfig, mhOptim) {

  RPE_parallelWrapper = function(job) {
    
    rval = RPE_mhOptim(data = job$data, modCfg = job$modCfg, runConfig=job$runConfig, mhOptim = job$mhOptim)
    
    rval$runConfig$seed = job$seed # TODO: This probably should have a special place
    
    rval
  }
  
  jobs = list()
  for (i in 1:parallelConfig$numChains) {
    jobs[[i]] = makeParallelJob(data = data, modCfg = modCfg, runConfig = runConfig, mhOptim = mhOptim)
  }
  
  startupMsg = paste0("Starting parallel parameter estimation of ", parallelConfig$numChains, " chains on ", parallelConfig$numNodes, " nodes.")
  
  jobRes = runParallelJobs(parallelConfig, jobs, jobFun=RPE_parallelWrapper, startupMsg = startupMsg)
  
  # Collect output and create CCM_Parallel object
  rval = list(chains = jobRes, parallelConfig = parallelConfig)
  class(rval) = c(class(rval), "CCM_Parallel")
  
  rval
  
}


continueSampling_parallel = function(parRes, iterations) {

  continueSampling_parallelWrapper = function(job) {
    
    rval = continueSampling(results=job$results, 
                            iterations=job$iterations, 
                            combinedOnly=TRUE)
    
    rval$runConfig$seed = job$seed
    
    rval
  }
  
  pc = parRes$parallelConfig
  
  # Set number of chains to number of results objects (if for some reason different)
  pc$numChains = length(parRes$chains)
  
  jobs = list()
  for (i in 1:pc$numChains) {
    jobs[[i]] = makeParallelJob(results=parRes$chains[[i]], iterations=iterations)
  }
  
  startupMsg = paste0("Continuing parallel parameter estimation of ", pc$numChains, " chains on ", pc$numNodes, " nodes.")
  
  jobRes = runParallelJobs(pc, jobs, continueSampling_parallelWrapper, startupMsg = startupMsg)
  
  # Collect output and create CCM_Parallel object
  rval = list(chains = jobRes, parallelConfig = pc)
  class(rval) = c(class(rval), "CCM_Parallel")
  
  rval
}

