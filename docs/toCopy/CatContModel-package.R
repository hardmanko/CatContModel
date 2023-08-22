#' CatContModel: A package for delayed estimation tasks with circular, continuous stimuli.
#' 
#' CatContModel is used for estimating the parameters of the model used by Hardman, Vergauwe, and Ricker (2017). This is a model for a continuous estimation task in which the stimuli exist in a circular space. The model is difficult to implement for a number of reasons, primarily that it has a large number of parameters, most of which do not have conjugate priors. This makes parameter estimation slow because the Metropolis-Hastings algorithm is used for estimating those parameters, which requires that the likelihood function be evaluated many times. Each evaluation of the likelihood function is costly, due to the structure of the model and the data. To improve parameter estimation speed, the core parameter estimation code is written in C++.
#' 
#' Be sure to read the manual, which can be opened with `CatContModelManual()`.
#' 
#' Some examples of how to use the package can be found in the manual. Additional examples can be found in the \href{https://github.com/hardmanko/CatContModel}{GitHub repository} for this package.
#' 
#' See also the [`Glossary`] of common terms.
#' 
#' The primary function for this package is [`runParameterEstimation`], which does the main parameter estimation. Other useful functions include [`examineMHAcceptance`],  [`removeBurnIn`], [`continueSampling`], [`combineResults`], [`plotParameterSummary`], [`testConditionEffects`], [`testMainEffectsAndInteractions`], among others. 
#' 
#'
#' @name CatContModel
#' @docType package
#' @md
#' 
#' @useDynLib CatContModel, .registration = TRUE
#' 
#' @importFrom Rcpp evalCpp 
NULL