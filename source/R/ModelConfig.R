



#' Make Model Configuration List
#' 
#' This function mostly exists as a place to document model configuration options.
#' The object returned by this function can be modified before calling [`runParameterEstimation`].
#' 
#' @param data The data to which you will fit the model. (Data is not stored in the model config, but is used for setting some defaults.)
#' @param modelVariant One of "betweenItem", "withinItem", "ZL" (Zhang & Luck 2008), or "betweenAndWithin" (not recommended).
#' @param dataType One of "circular" or "linear".
#' @param maxCategories For modelVariants with category parameters (catMu and catActive), the number of categories to use.
#' @param minSD The smallest value of standard deviation parameters. Use `NULL` for defaults that depend on `dataType`. See details.
#' @param catMuPriorApproximationPrecision The joint prior on `catMu` and `catActive` is not easily normalized so that it integrates to 1. 
#' The prior is rescaled to integrate to 1 by approximating the integral using a grid with `catMuPriorApproximationPrecision` points.
#' @param factors A data.frame explained in details. Use `NULL` for a one-factor design.
#' @param conditionEffects A named list where the names are parameters (e.g. "pMem") and the values are vectors of names of factors in `factors`. See details.
#' @param cornerstoneConditionName The name of the condition that will be used as the cornerstone condition. Defaults to the condition with the most observations in the data set.
#' @param responseRange Only used if `dataType=="linear"`. The possible range of response values as a length 2 numeric
#' vector, where the first value in the vector is the lower limit and the second the upper limit. 
#' By default, this is taken to be the range of the response values from the data set.
#' @param catMuRange The allowed range of `catMu` parameters. Defaults to `responseRange`, but there are reasons why you might want to allow categories to be outside of the response range.
#' @param calculateParticipantLikelihoods Whether (log) likelihoods should be calculated on each iteration for each participant. Setting to `TRUE` enables the use of [`calculateInappropriateFitStatistics`]. Recommend `FALSE`.
#' @param priorOverrides User-provided values for fixed priors. If only some of the priors are provided, the rest will be set to defaults. See [`getDefaultPriors()`] for naming/format of list.
#' @param mhTuningOverrides A list of overrides of the default tuning parameters for the Metropolis-Hastings algorithm. 
#' See the Tuning Metropolis-Hastings Acceptance Rates section of the manual for more information. 
#' See also [`getDefaultMHTuning()`] for naming/format of list.
#' @param constantParamValues A named list of constant values for parameters. See [`setConstantParameterValue`] and [`setConstantCategoryParameters`].
#' @param startingParamValues A named list of starting values for parameters with the same format as `constantParamValues`. 
#' By default, starting values will be randomly sampled and overdispersed so they can find their way to a converged distribution.
#' @param checkConfig If `TRUE`, [`checkModelConfig`] will be called.
#' @param ... Undocumented features. Do not use.
#' 
#' @return A list that can be passed as the `config` argument of [`runParameterEstimation`].
#' 
#' @details  `minSD` sets the smallest allowable standard deviation for standard deviation parameters (contSD, catSD, catSelectivity).
#' For circular data, smaller values increase the setup time and memory usage of the Von Mises look up table exponentially.
#' The default `minSD` for circular data is 1 degree and for linear data is 0.001 units.
#'  
#' The `factors` data.frame describes how the conditions in your data can be arranged into factors. 
#' See examples (for this function) for a case with 4 conditions arranged into a 2x2 factorial design.
#'  
#' The `conditionEffects` list gives the factors along which parameters are allowed to vary. Thus, it can be used to
#' create models that differ in terms of flexibility, enabling nested model comparison with [`calculateWAIC`].
#' Each element is the name of a parameter (e.g. "pMem") and the value of each element is a vector of factor names (which must be in `factors`).
#' The special value "all" means that the parameter varies freely with all factors while "none" means that the parameter 
#' has the same value regardless of factor levels.
#' See examples.
#' 
#' @section Prior Overrides:
#' Some things about the priors can be adjusted with the `priorOverrides` argument. 
#' The prior distributions are specified in Hardman, Vergauwe, and Ricker (2017) and the default fixed values controlling those prior distributions can be gotten with [`getDefaultPriors`]. 
#' For all of the parameters, the priors are set in the latent space. See the Priors section of the manual for more information.
#' 
#' Most participant level parameters have a hierarchical Normal prior with estimated mean, `mu`, and variance, `var`. 
#' The prior on `mu` is Normal with fixed mean and variance. 
#' The prior on `var` is Inverse Gamma with fixed alpha and beta. 
#' To set the priors, use the `priorOverrides` argument. 
#' For example, to set the priors on the distribution of `pMem`, you would set `pMem.mu.mu`, `pMem.mu.var`, `pMem.var.a`, and/or `pMem.var.b`. 
#' 
#' The condition effect parameters have Cauchy priors with constant location and scale. 
#' For example, for `pMem`, the location and scale priors are `pMem_cond.loc` and `pMem_cond.scale`. 
#' The locations all default to 0. The locations should almost certainly always be 0 unless you have a compelling argument to use a different value.
#' The scales are different for different parameter types (probability parameters vs standard deviation parameters). 
#'
#' 
#' @export
#' @seealso [`checkModelConfig`]
#' 
#' @examples 
#' myConds = c("A1", "A2", "B1", "B2")
#' factors = makeFactors(cond=myConds, 
#'     letterFact=c("A", "A", "B", "B"), 
#'     numberFact=c("1", "2", "1", "2"))
#' conditionEffects = list(pMem="all", contSD="letterFact", catSD="none")
#'
#' testData = SIM_sampleTestData("betweenItem", conds=myConds)
#' config = makeModelConfig(testData, "betweenItem", 
#'     factors=factors, conditionEffects=conditionEffects)
#' 
makeModelConfig = function(data, modelVariant, 
                           dataType = "circular", 
                           maxCategories = 16, minSD = NULL, catMuPriorApproximationPrecision = 60,
                           factors = NULL, conditionEffects = NULL, 
                           cornerstoneConditionName = NULL,
                           responseRange=NULL, catMuRange = NULL,
                           calculateParticipantLikelihoods = FALSE,
                           priorOverrides = list(),
                           mhTuningOverrides = list(),
                           constantParamValues = list(),
                           startingParamValues = list(),
                           checkConfig = TRUE,
                           ...) 
{
  
  cfg = list(modelVariant = modelVariant,
             dataType = dataType,
             maxCategories = maxCategories, 
             minSD = minSD, 
             catMuPriorApproximationPrecision = catMuPriorApproximationPrecision,
             factors = factors,
             conditionEffects = conditionEffects, 
             cornerstoneConditionName = cornerstoneConditionName,
             responseRange = responseRange, 
             catMuRange = catMuRange,
             calculateParticipantLikelihoods=calculateParticipantLikelihoods,
             priorOverrides = priorOverrides,
             mhTuningOverrides = mhTuningOverrides,
             constantParamValues = constantParamValues,
             startingParamValues = startingParamValues)
  
  # Load hidden configuration options
  dots = list(...)
  for (n in names(dots)) {
    cfg[[n]] = dots[[n]]
  }
  
  if (checkConfig) {
    cfg = checkModelConfig(data, cfg)
  }
  
  cfg
}




#' Check Model Configuration Values
#' 
#' This function is used to verify the correctness of a configuration list and to add in default values for all missing 
#' elements. This function is used internally by, e.g., [`runParameterEstimation`].
#' 
#' @param data The data that will be used as the `data` argument of [`runParameterEstimation`].
#' @param config A configuration to be used as the `config` argument of [`runParameterEstimation`].
#' @param immediateWarnings If `TRUE`, warnings will be printed immediately. Regardless of the value, warnings will also be stored for later access with `warnings()`.
#' @param checkData If `TRUE`, check the data as well as the model `config`.
#' 
#' @return An updated configuration list, possibly with additional items added or existing items modified.
#' 
#' @export
#' 
#' @examples 
#' config = list(modelVariant = "betweenItem")
#' 
#' config$factors = data.frame(
#'   cond = c("a1", "a2", "a3", "b1", "b2", "b3"),
#'   letters = rep(c("a", "b"), each = 3),
#'   numbers = rep(c("1", "2", "3"), 2)
#' )
#' 
#' config$conditionEffects = list(pMem = c("letters", "numbers"),
#'   pContBetween = "all", #shortcut to listing them all
#'   contSD = "letters",
#'   pCatGuess = "numbers")
#' # All unspecified parameters will be set to "none".
#' 
#' # Make fake data
#' data = data.frame(pnum = 1, cond = config$factors$cond,
#'   study = 0, response = 0)
#' 
#' config = checkModelConfig(data, config)
checkModelConfig = function(data, config, immediateWarnings = FALSE, checkData = TRUE) {
  
  usedConfigKeys = names(config)
  
  allAllowedConfigKeys = c("modelVariant", 
                           "cornerstoneConditionName", "conditionEffects", "factors",
                           "maxCategories", "minSD", "catMuPriorApproximationPrecision", 
                           "dataType", "responseRange", "catMuRange", 
                           "calculateParticipantLikelihoods", 
                           "priorOverrides", "mhTuningOverrides", "constantParamValues", "startingParamValues",
                           "sharedParameters",
                           "privateConfig")
  
  disallowedConfigKeys = usedConfigKeys[ !(usedConfigKeys %in% allAllowedConfigKeys) ]
  
  if (length(disallowedConfigKeys) > 0) {
    msg = paste0("The following configuration settings were used, but are not allowed: ", 
                paste(disallowedConfigKeys, collapse=", "))
    stop(msg)
  }


  
  ########################
  # modelVariant
  possibleModelVariants = c("betweenAndWithin", "betweenItem", "withinItem", "ZL")
  visibleModelVariants = c("betweenItem", "withinItem", "ZL")
  if (is.null(config$modelVariant)) {
    stop(paste0("config$modelVariant not set. Choose from one of: ", paste(visibleModelVariants, collapse = ", "), "."))
  }
  if (!(config$modelVariant %in% possibleModelVariants)) {
    stop(paste0("Invalid model variant '", config$modelVariant, "' selected. Choose from one of: ", 
                paste(visibleModelVariants, collapse = ", "), "."))
  }
  
  ####################
  # dataType
  
  if (is.null(config$dataType)) {
    config$dataType = "circular"
    logMsg("Note: config$dataType set to ", config$dataType)
  } else if (!(config$dataType %in% c("circular", "linear"))) {
    stop("Invalid value for config$dataType. Choose from one of \"circular\" or \"linear\".")
  }
  
  if (config$dataType == "linear") {
    
    if (!is.null(config$responseRange) && length(config$responseRange) != 2) {
      stop("responseRange provided, but it is not valid. It must be a length 2 vector of numeric values.")
    }
    
    if (is.null(config$responseRange)) {
      config$responseRange = range(data$response)
      logMsg("Note: config$responseRange set to ", paste(config$responseRange, collapse=", "))
    }
    
    if (config$responseRange[1] >= config$responseRange[2]) {
      stop("config$responseRange[1] >= config$responseRange[2]. The lower end of the response range must be lower!")
    }
    
    if (any(data$response < config$responseRange[1]) || any(data$response > config$responseRange[2])) {
      stop("Some responses are outside of responseRange. responseRange must contain all responses.")
    }
    
    if (is.null(config$catMuRange)) {
      config$catMuRange = config$responseRange
    }
    
    if (config$catMuRange[1] >= config$catMuRange[2]) {
      stop("config$catMuRange[1] >= config$catMuRange[2]. The lower end of the catMu range must be lower!")
    }
    
  }
  
  ####################
  # Data checks
  # Must come after dataType
  if (checkData) {
    requiredColumns = c("pnum", "cond", "study", "response")
    if (!all(requiredColumns %in% names(data))) {
      stop( paste0("The required columns are not in the data set. These columns are: ", paste(requiredColumns, collapse=", ")) )
    }
    
    if (config$dataType == "circular" && (all(data$study < 10) || all(data$response < 10))) {
      msg = "Data appear to be in radians rather than degrees. The data should be in degrees. See CatContModel::r2d()"
      #warning(msg, immediate. = TRUE)
      warning(msg)
    }
  }
  

  ##################################
  # cornerstoneConditionName
  if (is.null(config$cornerstoneConditionName)) {
    # By default, set cornerstone condition to the condition with the most data.
    dataCounts = stats::aggregate(study ~ cond, data, length) # This is not actually aggregating study, it's just for length
    config$cornerstoneConditionName = dataCounts$cond[which.max(dataCounts$study)]
    logMsg("Note: config$cornerstoneConditionName set to \"", config$cornerstoneConditionName, "\".")
  }
  if (!is.character(config$cornerstoneConditionName)) {
    config$cornerstoneConditionName = as.character(config$cornerstoneConditionName)
  }
  
  ######################
  # maxCategories
  if (config$modelVariant == "ZL") {
    if (!is.null(config$maxCategories) && config$maxCategories != 0) {
      logMsg("Note: config$maxCategories has been set to 0 because you are using the ZL modelVariant.")
    }
    config$maxCategories = 0
  }
  if (is.null(config$maxCategories)) {
    config$maxCategories = 16
    logMsg("Note: config$maxCategories set to default of ", config$maxCategories, ".")
  }
  
  ####################
  # minSD
  if (is.null(config$minSD)) {
    if (config$dataType == "circular") {
      config$minSD = 1
      msg = paste0("Note: config$minSD set to default of ", config$minSD, " degree(s)")
    } else if (config$dataType == "linear") {
      config$minSD = 0.001
      msg = paste0("Note: config$minSD set to default of ", config$minSD, " units")
    }
  	logMsg(msg)
  }
  
  ##########################################
  # catMuPriorApproximationPrecision
  if (is.null(config$catMuPriorApproximationPrecision)) {
    config$catMuPriorApproximationPrecision = 60
    logMsg("Note: config$catMuPriorApproximationPrecision not set. Set to ", 
              config$catMuPriorApproximationPrecision, " points at which the prior is evaluated.")
  }
  
  ###########################################
  # calculateParticipantLikelihoods
  if (is.null(config$calculateParticipantLikelihoods)) {
    config$calculateParticipantLikelihoods = FALSE
    logMsg("Note: config$calculateParticipantLikelihoods set to ", config$calculateParticipantLikelihoods, ".")
  }
  
  #####################
  # Priors
  # Missing priors are set in c++ so the config list is not full of priors.
  defaultPriors = getDefaultPriors()
  if (is.null(config$priorOverrides)) {
    config$priorOverrides = list()
  } else {
    for (n in names(config$priorOverrides)) {
      if (!(n %in% names(defaultPriors))) {
      	logWarning("Invalid prior name: ", n)
      }
    }
  }

  #########
  # TODO
  #config$priorOverrides = linearizeDefaultPriors(config$responseRange, config$priorOverrides)
  
  
  ###########################################
  # MH
  defaultMhTuning = getDefaultMHTuning()
  if (is.null(config$mhTuningOverrides)) {
    config$mhTuningOverrides = list()
  } else {
    for (n in names(config$mhTuningOverrides)) {
      if (!(n %in% names(defaultMhTuning))) {
        stop(paste0("Invalid MH tuning override name: ", n))
      }
    }
  }
  
  #######################
  # factors

  if (is.null(config$factors)) {
    config$factors = makeFactors(cond=sort(unique(data$cond)), check=FALSE)
    logMsg("Note: factors not provided. A one-factor design is assumed.")
  }
  config$factors = checkFactors(data, config$factors, factorsToCharacter = TRUE)
  

  ######################
  # Check starting and constant values
  
  for (n in names(config$startingParamValues)) {
    if (grepl("catActive", n, fixed=TRUE)) {
      if (!(config$startingParamValues[[n]] %in% c(0, 1))) {
        stop("Starting values of catActive must be either 0 or 1.")
      }
    }
  }
  
  for (n in names(config$constantParamValues)) {
    if (grepl("catActive", n, fixed=TRUE)) {
      if (!(config$constantParamValues[[n]] %in% c(0, 1))) {
        stop("Constant values of catActive must be either 0 or 1.")
      }
    }
  }
  
  
  ######################
  # conditionEffects
  config$conditionEffects = verifyConditionEffects(config, immediateWarnings=immediateWarnings)
  
  
  ######################
  # Single Participant
  if (length(unique(data$pnum)) == 1) {
    logMsg("Single participant design detected. Removing hierarchical priors. See documentation for singleParticipantConfiguration() for more information.")
    
    config = singleParticipantConfiguration(config)
  }
  
  ######################
  # Shared Participant Parameters
  if (is.null(config$sharedParameters)) {
    config$sharedParameters = character(0)
  } else {
    
    # TODO: Clean this whole section up. config$model????
    
    # This is for MemFirst
    allowedShared = c("catMu", "catActive")
    
    
    # If CatFirst
    if (!is.null(config$privateConfig$CatFirst)) {
      allowedShared = c(getParamNames(), "catType")
    }
    
    disallowedUsed = config$sharedParameters[ !(config$sharedParameters %in% allowedShared) ]
    if (length(disallowedUsed) > 0) {
      stop(paste0("The following parameters are not allowed to be shared: ", paste0(disallowedUsed, collapse=", ")))
    }
    
    # Sharing dependencies:
    catMuShared = "catMu" %in% config$sharedParameters
    catActiveShared = "catActive" %in% config$sharedParameters
    catTypeShared = "catType" %in% config$sharedParameters
    if (!catMuShared && (catActiveShared || catTypeShared)) {
      if (catActiveShared) {
        stop("catActive cannot be shared if catMu is not shared.")
      }
      if (catTypeShared) {
        stop("catType cannot be shared if catMu is not shared.")
      }
    }
  }
  
  config
}



#' Define Factors Levels for Conditions
#' 
#' Convenience function for defining how conditions (`cond`) map to factor levels.
#' See `factors` argument of [`makeModelConfig`].
#' 
#' @param cond Names of conditions in the `cond` column of `data`.
#' @param ... Named arguments, each a character vector with the same length as `condNames`. The name is taken to be the name of a factor. 
#' The values in the vector are levels of the factor that correspond to levels of `cond`. See examples.
#' @param data The data you will fit. Used to check that the conditions are in the data. Can be `NULL` or missing to skip checks.
#' @param check If `TRUE`, the factors will be checked.
#' 
#' @return A `data.frame` that can be passed as the `factors` argument of [`makeModelConfig`].
#' 
#' @export
#' @examples 
#' factors = makeFactors(cond=c("A1", "A2", "B1", "B2"), 
#'   letterFact=c("A", "A", "B", "B"), 
#'   numberFact=c("1", "2", "1", "2"))
makeFactors = function(cond, ..., data=NULL, check=TRUE) {
  
  factorMap = list(...)
  
  if (length(factorMap) == 0) {
    factorMap$DefaultFactor = cond
    #logMsg("Note: config$factors not provided. A one-factor design is assumed.")
  }
  
  # Copy over factors
  rval = data.frame(cond=cond, stringsAsFactors = FALSE)
  
  if (any(names(factorMap) == "")) {
    stop("Factors must have a non-empty name.")
  }

  if (any(names(factorMap) == "cond")) {
    stop('Factors cannot be named "cond".')
  }
  
  for (factName in names(factorMap)) {
    # Check length
    if (length(factorMap[[factName]]) != length(cond)) {
      stop(paste0("The length of factor ", factName, " is different from the length of cond."))
    }
    
    rval[,factName] = factorMap[[factName]]
  }
  
  if (check) {
    rval = checkFactors(data=data, factors=rval)
  }
  
  rval
}


checkFactors = function(data, factors, factorsToCharacter=TRUE) {
  
  # C1: Check condition names
  if (!is.null(data)) {
    dataConds = unique(data$cond)
    allConds = union( dataConds, factors$cond )
    
    for (cond in allConds) {
      if (!(cond %in% factors$cond)) {
        stop(paste0('The condition "', cond, '" is in the data but not in factors.'))
      }
      if (!(cond %in% dataConds)) {
        stop(paste0('The condition "', cond, '" is in factors but not in the data.'))
      }
    }
  }
  
  # C2: Convert factors to character (in case they were R factors)
  if (factorsToCharacter) {
    for (n in names(factors)) {
      factors[ , n ] = as.character(factors[ , n ])
    }
  }
  
  # C3: Check for illegal factor names and factor level names
  if (any(c("key", "group") %in% names(factors))) {
    stop('Factors cannot be named either "key" or "group".')
  }
  
  for (n in names(factors)) {
    nameBad = grepl("[:.]+", n)
    levelsBad = any(grepl("[:.]+", factors[,n]))
    if (nameBad || levelsBad) {
      stop('Factor names and factor levels may not contain colon (":") or period (".").')
    }
  }
  
  # C4: Factors contains cond
  if (!("cond" %in% names(factors))) {
    stop('Factors must have a "cond" column.')
  }
  
  # C5: Elements in cond are unique
  if (length(factors$cond) != length(unique(factors$cond))) {
    stop("Duplicate cond in factors.")
  }
  
  factors
}


# Checks the correctness of the condition effects in 
# config$conditionEffects, including possibly modifying them.
# This function is used by checkModelConfig
verifyConditionEffects = function(config, immediateWarnings = FALSE) {
  
  parametersWithPossibleConditionEffects = getParamNames(config$modelVariant, types=c("prob", "sd"))
  
  if (is.null(config$conditionEffects)) {
    
    pceToUse = getDefaultParametersWithConditionEffects(config$modelVariant)
    
    config$conditionEffects = list()
    for (parName in pceToUse) {
      config$conditionEffects[[ parName ]] = "all"
    }
    
    logMsg("Note: config$conditionEffects not set. It was set to use all factors for ", paste(pceToUse, collapse = ", "), ".")
    
  } else { 
    
    for (n in names(config$conditionEffects)) {
      if (!(n %in% parametersWithPossibleConditionEffects)) {
        config$conditionEffects[[n]] = NULL
        msg = paste0("In config$conditionEffects, parameter \"", n, "\" was included, but it is not a parameter that can have condition effects (possibly because it is not used by the current modelVariant). Its condition effects have been ignored.")
        if (immediateWarnings) {
          warning(msg, immediate. = TRUE)
        }
        warning(msg)
      }
    }
  }
  
  #set all unmentioned condition effects to "none"
  unmentionedCE = c()
  for (n in parametersWithPossibleConditionEffects) {
    if (!(n %in% names(config$conditionEffects))) {
      config$conditionEffects[[n]] = "none"
      unmentionedCE = c(unmentionedCE, n)
    }
  }
  if (length(unmentionedCE) > 0) {
    logMsg('In conditionEffects, the following parameters were set to the default value of "none": ', 
    			 paste(unmentionedCE, collapse=", "))
  }
  
  #Double check that condition effect factor names are in factors
  factNames = getAllFactorNames(config$factors)
  for (n in names(config$conditionEffects)) {
    ce = config$conditionEffects[[n]]
    if (sum(ce %in% c("all", "none")) >= 2) {
      stop(paste0("config$conditionEffects$", n, " contains more than one instance of \"all\" or \"none\"."))
    }
    if (length(ce) > 1 || !all(ce %in% c("all", "none"))) {
      notIn = ce[ !(ce %in% factNames) ]
      if (length(notIn) > 0) {
        stop( paste0("config$conditionEffects$", n, " contains factor names not found in config$factors. The bad factor name(s): ", paste(notIn, collapse=", "), ".") )
      }
    }
  }
  
  # Check that config$conditionEffects are correct, given the constant value overrides that are being used.
  # This was checkConditionEffectsGivenConstantParameters
  for (parName in names(config$conditionEffects)) {
    
    if (length(config$conditionEffects[[parName]]) > 1 || config$conditionEffects[[parName]] != "none") {
      thisCPN = paste(parName, "_cond[", config$factors$cond, "]", sep="")
      
      inCPV = thisCPN %in% names(config$constantParamValues)
      
      if (all(inCPV)) {
        logMsg("Note: Parameter ", parName, " has constant value overrides on all condition effect parameters. It has been noted to have no condition parameters.")
        
        config$conditionEffects[[parName]] = "none"
        
      } else if (any(inCPV)) {
        msg = paste0("Parameter ", parName, " has constant value overrides on some condition effect parameters. It will have condition effects estimated, but some follow-up tests may not work correctly. After parameter estimation is complete, consider setting results$config$conditionEffects$", parName, " to \"none\".")
        if (immediateWarnings) {
          warning( msg, immediate. = TRUE )
        }
        warning( msg )
      }
    }
    
  }
  
  config$conditionEffects
}



# Converts units of starting and constant catMu, typically from degrees to radians
convertModelConfigUnits = function(config, convFun) {
  
  for (n in names(config$startingParamValues)) {
    if (grepl("catMu", n, fixed=TRUE)) {
      config$startingParamValues[[n]] = convFun(config$startingParamValues[[n]])
    }
  }
  
  for (n in names(config$constantParamValues)) {
    if (grepl("catMu", n, fixed=TRUE)) {
      config$constantParamValues[[n]] = convFun(config$constantParamValues[[n]])
    }
  }

  config
}


#' Configure Model for Single Participant Design
#' 
#' If you only have data for one participant, the hierarchical priors need to be removed. 
#' This function does that by setting the normally-estimated parameters of the hierarchical priors to constant values. See details.
#' 
#' @param config A model configuration. See [`makeModelConfig`].
#' 
#' @details
#' All of the standard parameters (all but catMu or catActive) have hierarchical priors.
#' For a participant parameter `P_i`: `P_i ~ Normal(P.mu, P.var)`. 
#' 
#' Normally, `P.mu` and `P.var` are estimated based on the participant parameters and have prior distributions controlled by constant values.
#' + P.mu ~ Normal(`P.mu.mu`, `P.mu.var`)
#' + P.var ~ InverseGamma(`P.var.a`, `P.var.b`).
#' Call [`getDefaultPriors`] to see default values for `P.mu.mu`, etc.
#' 
#' This function sets `P.mu` and `P.var` to constant values (see `constantParamValues` in the returned list).
#' The values for `P.mu` and `P.var` are set equal to the values for `P.mu.mu` and `P.mu.var` for all standard parameters.
#' The values of `P.mu.mu` and `P.mu.var` are gotten from `config$priorOverrides`. 
#' If there is no prior override for a parameter, the value is gotten from [`getDefaultPriors`].
#' 
#' @return An updated model configuration.
#' 
#' @export
singleParticipantConfiguration = function(config) {
  
  pNames = getParamNames(config$modelVariant, types=c("prob", "sd"))
  
  priors = getDefaultPriors()
  
  for (pon in names(config$priorOverrides)) {
    priors[[ pon ]] = config$priorOverrides[[ pon ]]
  }
  
  for (pn in pNames) {
    # Do not overwrite constant values.
    
    muPN = paste0(pn, ".mu")
    if (!(muPN %in% names(config$constantParamValues))) {
      config$constantParamValues[[ muPN ]] = priors[[ paste0(pn, ".mu.mu") ]]
    }
    
    varPN = paste0(pn, ".var")
    if (!(varPN %in% names(config$constantParamValues))) {
      config$constantParamValues[[ varPN ]] = priors[[ paste0(pn, ".mu.var") ]]
    }

  }
  
  config
}




