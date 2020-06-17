
#' Set Parameters to Constant Values
#' 
#' This function helps with setting parameters to constant values. It returns a
#' list mapping from parameter name to parameter value. It sets all participant
#' parameters to the same value. It also sets the hierarchical, population
#' level parameters to constant values (this doesn't really have any effect:
#' if all of the participant level parameters are constant, the population level
#' parameters do nothing).
#' 
#' If `zeroConditionEffects` is `TRUE`, the condition effects are all set to 0, which
#' means that the participants will have the same parameter values in all conditions.
#' Note that this overrides the `conditionEffects` setting in the
#' configuration for [`runParameterEstimation`]. If you allow varying condition
#' effects with `config$conditionEffects` but set `zeroConditionEffects` to `TRUE`
#' you will get no condition effects. The reverse is not true: Setting 
#' `zeroConditionEffects` to `FALSE` does not create varying condition effects
#' in `config$conditionEffects`.
#' 
#' If you need to call this function several times to get multiple lists of fixed parameter
#' values, you can combine the lists with the concatenate function, `c()`.
#' 
#' @param data The data you will use, in the same format as required by [`runParameterEstimation`].
#' @param param The parameter to set to a constant value, e.g. `"pMem"`.
#' @param value The constant value of the paramter. See also the `transformValueToLatent` argument.
#' @param zeroConditionEffects Whether condition effects should all be set to 0. (See [`Glossary`].)
#' @param transformValueToLatent Whether `value` should be transformed to the latent space. If `TRUE` (the default) you should provide parameter values in the manifest space (e.g. probabilities should be between 0 and 1). If `FALSE`, you should provide parameter values in the latent space.
#' 
#' @return A list mapping from parameter name to parameter value. 
#' This list can be passed as the `constantValueOverrides` argument of [`runParameterEstimation`].
#' 
#' @seealso [`setConstantCategoryParameters`] for setting catMu and catActive to constant values.
#' 
#' @examples
#' data = data.frame(pnum=rep(1:5, each=2), cond=rep(1:2, each=5))
#' constParam = setConstantParameterValue(data, param = "pContBetween", value = 0.5)
#' 
#' @export
setConstantParameterValue = function(data, param, value, zeroConditionEffects = TRUE, transformValueToLatent = TRUE) {
	
	if (transformValueToLatent) {
		trans = getParameterTransformation(NULL, param, inverse=TRUE)
		value = trans(value)
		if (param %in% getProbParams(NULL)) {
			value = min(max(value, -100), 100) #clamp probability parameters to be finite.
		}
	}
	
	pnums = sort(unique(data$pnum))
	conds = sort(unique(data$cond))
	
	rval = list()
	for (pnum in pnums) {
		partParam = paste( param, "[", pnum, "]", sep="")
		rval[[ partParam ]] = value
	}
	
	rval[[ paste( param, ".mu", sep="" ) ]] = 0
	rval[[ paste( param, ".var", sep="" ) ]] = 1
	
	if (zeroConditionEffects) {
		for (cond in conds) {
			condParam = paste( param, "_cond[", cond, "]", sep="")
			rval[[ condParam ]] = 0
		}
	}
	
	rval
}

#' Set catMu and catActive Parameters to Constant Values
#' 
#' 
#' @param data Your data set, in the same format as required by [`runParameterEstimation`].
#' @param catParam A `data.frame` with constant values of `catMu` and `catActive` parameters. See details.
#' @param maxCategories The maximum number of categories that will be allowed when running parameter estimation.
#' @param activateConstantCats If `TRUE`, categories for which `catMu` values are given are forced to be active.
#' @param deactivateUnspecifiedCats If `TRUE`, `catActive` is set to 0 (inactive) for all categories for which constant values of both `catMu` and `catActive` are not specified. That is, if `maxCategories` is `5` but you only list 3 categories per participant, the additional 2 categories are set to be inactive. If `FALSE`, `catActive` and `catMu` are left as free parameters for the extra, unspecified categories.
#' 
#' @return A list mapping from parameter name to parameter value which can be passed as the `constantValueOverrides` argument of `runParameterEstimation`.
#' 
#' @details The `catParam` `data.frame` can have three columns: `pnum`, `catMu`, and `catActive`, where the `catActive` column is optional. Each row specifies the setting for a `catMu` and/or `catActive` for the given participant. Each participant can have multiple rows, each row specifying constant parameters for a different category. If you want to allow a parameter to be estimated, set the value for that parameter to `NA`. 
#' 
#' Because either or both of `catMu` and `catActive` can be specified, there are 4 possibilities per category: 
#' 1) both `catMu` and `catActive` set to constant values. This forces the model to use (or not use) the category at a fixed location.
#' 2) `catMu` set to a constant value but `catActive` freely estimated. In this case, what `catActive` will tell you is mow much of the time a category at that fixed location was used. You might do this if you have candidates for category locations that should be used, but don't know how much each category might be used.
#' 3) `catActive` set to a constant value, but `catMu` freely estimated. You should not set `catActive` to 0 in this case, because `catMu` cannot be meaningfully estimated for inactive categories. Thus, by setting `catActive` to 1 but leaving `catMu` unspecified, it forces the model to use a category, but allows it to put that category anywhere.
#' 4) Neither `catMu` nor `catActive` specified. This is the default behavior of the model, where both are freely estimated.
#' 
#' Each of these 4 possibilities is set per category. Thus, you can do complex things like force each participant to have 3 active categories in fixed locations, 2 more active categories in freely-estimated locations, and some additional number of fully freely-estimated categories.
#' 
#' @examples 
#' catParam = data.frame(pnum = rep(1, 4),
#' 	catMu = c(60, 120, NA, NA),
#' 	catActive = c(1, NA, 1, NA))
#' data = data.frame(pnum=1) #This is just for testing
#' constCat = setConstantCategoryParameters(data, catParam, maxCategories = 5, 
#' 	activateConstantCats = TRUE, deactivateUnspecifiedCats = TRUE)
#' 
#' @export
setConstantCategoryParameters = function(data, catParam, maxCategories, activateConstantCats = TRUE, 
																				 deactivateUnspecifiedCats = TRUE) 
{
	
	allPnums = sort(unique(data$pnum))
	
	if (is.null(catParam$catActive)) {
		catParam$catActive = NA
	}
	
	if (activateConstantCats) {
		catParam[ !is.na(catParam$catMu), "catActive" ] = 1.0
	}
	
	if (!all(catParam$catActive %in% c(0.0, 1.0, NA))) {
		stop("All catActive parameters must be either 0 or 1, or NA if not constant.")
	}
	
	rval = list()
	for (pnum in allPnums) {
		
		#mus = catParam$catMu[ catParam$pnum == pnum ]
		catMu = catParam$catMu[ catParam$pnum == pnum ]
		catActive = catParam$catActive[ catParam$pnum == pnum ]
		
		if (length(catMu) > maxCategories) {
			warning( paste("More categories provided for participant ", pnum, " than the maximum number of categories.", sep="") )
		}
		
		for (i in 1:maxCategories) {
			
			catInd = i - 1 #zero-indexed
			
			catMuName = paste( "catMu[", pnum, ",", catInd, "]", sep="")
			catActiveName = paste( "catActive[", pnum, ",", catInd, "]", sep="")
			
			if (!is.na(catMu[i])) {
				rval[[ catMuName ]] = catMu[i]
			}
			
			if (!is.na(catActive[i])) {
				rval[[ catActiveName ]] = catActive[i]
			}
			
			if (deactivateUnspecifiedCats && is.na(catActive[i]) && is.na(catActive[i])) {
				rval[[ catMuName ]] = 0.0
				rval[[ catActiveName ]] = 0.0
			}
			
		}
	}
	
	rval
}

