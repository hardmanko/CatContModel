

configureOutputFile = function(clust, numNodes, outputFile) {
  
  if (!is.null(outputFile) && outputFile != "") {
    outExt = substr(outputFile, nchar(outputFile) - 3, nchar(outputFile))
    
    if (outExt == ".txt") {
      # If full file name, use single output file only for the first node.
      clusterApply(clust, x=as.list(1:numNodes), fun=function(x, ...) {
        if (x == 1) {
          sink(outputFile, append = FALSE)
        }
      }, outputFile=outputFile)
      
    } else {
      # If no .txt extension, use individual output files for each node.
      clusterApply(clust, x=as.list(1:numNodes), fun=function(x, ...) {
        sink(paste0(list(...)$outputPrefix, x, ".txt"), append = FALSE)
      }, outputPrefix=outputFile)
    }
  }
  
}

verifyOrSampleSeeds = function(seeds, numNodes) {
  if (!is.null(seeds) && length(seeds) != numNodes) {
    warning("Wrong number of seeds provided. They will be ignored.")
    seeds = NULL
  }
  if (is.null(seeds)) {
    seeds = sample.int(.Machine$integer.max, numNodes, replace=FALSE)
  }
  as.list(seeds)
}


runParameterEstimation_parallelWrapper = function(seed, ...) {
  
  require(CatContModel)
  
  set.seed(seed)
  
  dots = list(...)
  for (listArgs in c("mhTuningOverrides", "priorOverrides", "startingValueOverrides", "constantValueOverrides")) {
    if (is.null(dots[[listArgs]])) {
      dots[[listArgs]] = list()
    }
  }
  
  rval = runParameterEstimation(config = dots$config,
                                data = dots$data,
                                mhTuningOverrides = dots$mhTuningOverrides,
                                priorOverrides = dots$priorOverrides,
                                startingValueOverrides = dots$startingValueOverrides,
                                constantValueOverrides = dots$constantValueOverrides)
  
  rval
}


#' Run Parameter Estimation in Parallel
#' 
#' There are two major advantages to parallel estimation: Speed and convergence diagnostics.
#' Speed: By running chains in parallel, more iterations of the parameter chains can be sampled more quickly.
#' Convergence: If multiple chains, each with different starting values, all converge to the same
#' distributions, that is a strong convergence diagnostic.
#' 
#' @param numNodes The number of cluster nodes to use. Each node will run one chain.
#' @param config See [`runParameterEstimation`].
#' @param data See [`runParameterEstimation`].
#' @param mhTuningOverrides See [`runParameterEstimation`].
#' @param priorOverrides See [`runParameterEstimation`].
#' @param startingValueOverrides See [`runParameterEstimation`].
#' @param constantValueOverrides See [`runParameterEstimation`].
#' @param seeds A vector of length `numNodes` of RNG seeds that are used to initialize each node's RNG with `set.seed`. If `NULL`, seeds will be sampled.
#' @param nodeSetupFun A function that is called on each node before running parameter estimation. Not usually needed.
#' @param outputFile A text file name (e.g. "progress.txt") or an incomplete file name without .txt extension (e.g. "progress_"). See details.
#' @param cluster A cluster created with `parallel::makeCluster()`. If `NULL`, a cluster will be created with `numNodes` nodes.
#' 
#' @details The parallel parameter estimation functions use the R package `parallel` to run multiple 
#' parallel parameter chains, not to speed estimation of one parameter chain.
#' Once the parallel chains have been sampled, 
#' 1) each chain must have burn-in iterations removed (see [`removeBurnIn`]) and 
#' 2) the chains must be merged into a single results object (see [`mergeResults`]).
#' 
#' If `cluster` is not provided, a cluster will be created on the local machine using `numNodes` parallel R processes.
#' For best performance (shortest runtime), do not exceed the number of physical CPU cores on your machine. 
#' Note: `parallel::detectCores` returns logical cores, not physical cores.
#' 
#' Argument `outputFile` specifies how nodes should output their progress. 
#' Progress is normally printed to the console, but that doesn't work with multiple nodes.
#' If `outputFile` is a text file with extension ("runProgress.txt"), output from only one of the nodes is printed to that file.
#' If `outputFile` is an incomplete file name ("runProgress_"), one file will be created for each node (e.g. "runProgress_1.txt").
#' If `outputFile` is empty ("") or NULL, no output file is created.
#' File contents are overwritten each time a new run starts.
#' 
#' For between-participants designs, it would be helpful to be able to run different participant groups in parallel,
#' but the parallel parameter estimation functions do not help with that. You will still need to run estimation for
#' each group separately. You can use the parallel functions to run parallel chains for each group, however, which can
#' speed parameter estimation.
#' 
#' @seealso [`continueSampling_parallel`]
#' 
#' @return A list with `numNodes` elements, each element being the return value from [`runParameterEstimation`].
#' 
#' @export
runParameterEstimation_parallel = function(numNodes, 
                                           config, 
                                           data, 
                                           mhTuningOverrides = list(),
                                           priorOverrides = list(),
                                           startingValueOverrides = list(),
                                           constantValueOverrides = list(), 
                                           seeds = NULL,
                                           nodeSetupFun = NULL,
                                           outputFile = NULL,
                                           cluster = NULL)
{
  
  require(parallel)
  
  seedsList = verifyOrSampleSeeds(seeds, numNodes)
  
  # Make the cluster if needed
  createCluster = is.null(cluster)
  if (createCluster) {
    cluster = makeCluster(numNodes)
  }
  
  configureOutputFile(cluster, numNodes, outputFile)
  
  # Set up each cluster instance in the same way
  if (!is.null(nodeSetupFun)) {
    clusterCall(cluster, nodeSetupFun)
  }
  
  # Run the parameter estimation
  clustResults = clusterApply(cluster, seedsList, runParameterEstimation_parallelWrapper, 
                              config=config, data=data, 
                              mhTuningOverrides=mhTuningOverrides,
                              priorOverrides=priorOverrides, 
                              startingValueOverrides=startingValueOverrides, 
                              constantValueOverrides=constantValueOverrides)
  
  # Done: Stop the cluster
  if (createCluster) {
    stopCluster(cluster)
  }
  
  # Add seed to the results
  for (i in 1:numNodes) {
    clustResults[[i]]$info$randomSeed = seedsList[[i]]
  }
  
  clustResults
}


continueSampling_parallelWrapper = function(nodeConfig, ...) {
  
  require(CatContModel)
  
  set.seed(nodeConfig$seed)
  
  rval = continueSampling(results=nodeConfig$results, 
                          iterations=nodeConfig$iterations, 
                          combinedOnly=TRUE)
  
  rval
}


#' Continue Sampling from Parallel Parameter Chains
#' 
#' This is the parallel version of [`continueSampling`].
#' 
#' There is a difference between the return value of [`continueSampling`] and `continueSampling_parallel`:
#' `continueSampling` returns old results, new results, and merged results while `continueSampling_parallel` only returns merged results. 
#' You can easily get parallel old results and new results: `continueSampling_parallel` takes the old results as its argument and the new results can be gotten by removing the old results from the merged results with [`removeBurnIn`].
#' 
#' @param resList A list of results objects (like returned by [`runParameterEstimation_parallel`]).
#' @param iterations The number of new iterations to sample for each chain.
#' @param seeds See [`runParameterEstimation_parallel`].
#' @param nodeSetupFun See [`runParameterEstimation_parallel`].
#' @param outputFile See [`runParameterEstimation_parallel`].
#' @param cluster See [`runParameterEstimation_parallel`].
#' 
#' @return A list of results objects (like `resList`) with additional iterations sampled.
#' 
#' @export
continueSampling_parallel = function(resList, iterations, 
                                       seeds=NULL, 
                                       nodeSetupFun=NULL, 
                                       outputFile = NULL,
                                       cluster = NULL) 
{
  require(parallel)
  
  numNodes = length(resList)
  
  seedsList = verifyOrSampleSeeds(seeds, numNodes)
  
  # Make the cluster
  createCluster = is.null(cluster)
  if (createCluster) {
    cluster = makeCluster(numNodes)
  }
  
  configureOutputFile(cluster, numNodes, outputFile)
  
  # Set up each cluster instance in the same way
  if (!is.null(nodeSetupFun)) {
    clusterCall(cluster, nodeSetupFun)
  }
  
  # Continue sampling
  nodeConfig = list()
  for (i in 1:numNodes) {
    nodeConfig[[i]] = list(seed=seedsList[[i]], results=resList[[i]], iterations=iterations)
  }
  
  clustResults = clusterApply(cluster, nodeConfig, continueSampling_parallelWrapper)
  
  # Done: Stop the cluster
  if (createCluster) {
    stopCluster(cluster)
  }
  
  # Add seed to the results
  for (i in 1:numNodes) {
    clustResults[[i]]$info$randomSeed = seedsList[[i]]
  }
  
  clustResults
  
}
