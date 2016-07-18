#' CatContModel: A package for delayed estimation tasks with circular, continuous stimuli.
#' 
#' CatContModel is used for estimating the parameters of the model used by Hardman, Vergauwe, and Ricker (In press). This is a model for a continuous estimation task in which the stimuli exist in a circular space. The model is difficult to implement for a number of reasons, primarily that it has a large number of parameters, most of which do not have conjugate priors. This makes parameter estimation slow because the Metropolis-Hastings algorithm is used for estimating those parameters, which requires that the likelihood function be evaluated many times. Each evaluation of the likelihood function is costly, due to the structure of the model and the data. To improve parameter estimation speed, the core parameter estimation code is written in C++.
#' 
#' The primary function for this package is \link{runParameterEstimation}, which does the main parameter estimation. Other useful functions include \link{examineMHAcceptance}, \link{continueSampling}, \link{removeBurnIn}, \link{plotParameterSummary}, \link{mergeResults}, \link{testConditionEffects}, among others.
#' 
#'
#' @name CatContModel
#' @docType package
#' 
#' @useDynLib CatContModel
#' 
#' @importFrom Rcpp evalCpp 
NULL