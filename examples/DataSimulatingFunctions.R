
library(CatContModel)
library(msm) #for rtnorm


makeCategories = function(nCat, maxCat, noiseScale = 1, catRange = c(0, 360)) {
	if (nCat == 0) {
		return(NULL)
	}
	if (nCat > maxCat) {
		nCat = maxCat
		warning("nCat > maxCat")
	}
	
	rangeWidth = catRange[2] - catRange[1]
	
	noiseRange = noiseScale * rangeWidth / maxCat / 2

	baseLocations = seq(catRange[1], catRange[2] * (maxCat - 1)/maxCat, length.out=maxCat)
	
	catMu = sample(baseLocations, nCat)
	
	catMu = catMu + runif(nCat, -noiseRange, noiseRange)

	#get into range
	catMu = ((catMu - catRange[1]) %% rangeWidth) + catRange[1]
	
	catMu
}


sampleSimulatedData = function(conditionEffects, partParam, trialsPerCondition, 
															 modelVariant="betweenAndWithin", dataType = "circular", 
															 minSd = 1, studyRange = NULL, responseRange = NULL) {
	
	nCond = length(conditionEffects$pMem)
	nPart = length(partParam$pMem)
	
	ppList = list()
	
	if (modelVariant == "ZL") {
		partParam$nCat = rep(0, nPart)
	}
	
	for (i in 1:nPart) {
		
		pp = list()
		
		for (pn in names(partParam)) {
			pp[[pn]] = partParam[[pn]][i]
		}

		if (partParam$nCat[i] > 0) {
			pp$catMu = partParam$catMu[i,1:partParam$nCat[i]]
		}
		
		ppList[[i]] = pp
		
	}
	
	combinedParam = list()
	for (cond in 1:nCond) {
		combinedParam[[cond]] = ppList
	}
	
	probParams = getProbParams(NULL)
	sdParams = getSdParams(NULL)
	
	#apply condition param
	for (i in 1:nPart) {
		
		for (cond in 1:nCond) {
			
			for (pp in probParams) {
				lg = logit(ppList[[i]][[pp]]) #Transform back to latent space
				ce = conditionEffects[[pp]][cond]
				
				combinedParam[[cond]][[i]][[pp]] = logitInverse(lg + ce)
			}
			
			for (sp in sdParams) {
				part = ppList[[i]][[sp]]
				ce = conditionEffects[[sp]][cond]
				
				combinedParam[[cond]][[i]][[sp]] = pmax(part + ce, minSd)
			}
		}
	}
	
	if (dataType == "circular") {
		responseRange = c(0, 360)
	}
	
	if (is.null(studyRange)) {
		studyRange = responseRange
	}
	
	
	#Simulate the data
	simulatedData = NULL
	
	for (pn in 1:nPart) {
		for (cond in 1:nCond) {
			trials = trialsPerCondition[cond]
			
			param = combinedParam[[cond]][[pn]]
			
			samp = sampleDataFromModel(runif(trials, studyRange[1], studyRange[2]), param, modelVariant=modelVariant, 
																 dataType=dataType, responseRange=responseRange)
			
			data = data.frame(pnum=pn, cond=cond, study=samp$study, response=samp$response)
			
			
			simulatedData = rbind(simulatedData, data)
		}
	}
	
	list(simulatedData = simulatedData, combinedParam = combinedParam)
}

convertParameterListToDataFrame = function(paramList, maxCat) {
	
	nPart = length(paramList[[1]])
	
	df = NULL
	
	for (pnum in 1:nPart) {
		for (condInd in 1:length(paramList)) {
			
			pdata = paramList[[condInd]][[pnum]]
			
			for (n in names(pdata)) {
				
				if (length(pdata[[n]]) == 0) {
					pdata[[n]] = NA
				}
				if (n == "catSD" || n == "catSelectivity") {
					pdata[[n]] = pdata[[n]][1]
				}
				
				temp = data.frame(pnum = pnum, cond = condInd, param = n, cat = NA, value = pdata[[n]], stringsAsFactors = FALSE)
				if (n == "catMu") {
					temp$cat = 1:nrow(temp) - 1
				}
				df = rbind(df, temp)
			}
		}
	}
	
	if (length(unique(df$cond)) > 1) {
		for (n in unique(df$param)) {
			thisParam = df[ df$param == n, ]
			
			conds = unique(thisParam$cond)
			firstCond = conds[1]
			otherConds = conds[2:length(conds)]
			
			allEqual = TRUE
			for (cond in otherConds) {
				eq = thisParam$value[ thisParam$cond == firstCond ] == thisParam$value[ thisParam$cond == cond ]
				allEqual = allEqual && all(eq)
			}
			
			if (!is.na(allEqual) && allEqual) {
				df = df[ !(df$param == n & df$cond != firstCond), ]
			}
		}
	}
	
	df = df[ order(df$param, df$cond, df$pnum), ]
	
	df = df[ !is.na(df$value), ]
	
	df
}

compareTrueAndRecovered = function(results, trueParam) {
	recPart = participantPosteriorSummary(results)
	
	trueParam$pnum = as.character(trueParam$pnum)
	trueParam$cond = as.character(trueParam$cond)
	trueParam$param = as.character(trueParam$param)
	
	trueParam = trueParam[ order(trueParam$pnum), ]
	
	trueParam$param[ trueParam$param == "nCat" ] = "catActive"
	
	trueParam = trueParam[ trueParam$param != "catMu", ]
	
	df = NULL
	
	commonParameters = intersect( unique(trueParam$param), unique(recPart$param) )
	
	for (param in commonParameters) {
		
		conds = unique(trueParam$cond[ trueParam$param == param ])
		
		for (cond in conds) {
			rec = recPart$mean[ recPart$param == param & recPart$cond == cond ]
			true = trueParam$value[ trueParam$param == param & trueParam$cond == cond ]
			
			corr = cor(rec, true)
			slope = as.numeric(coef(lm(rec ~ true))[2])
			dif = mean(rec - true)
			percentDif = mean(rec - true) / mean((rec + true) / 2) * 100
			
			temp = data.frame(param = param, cond = cond, 
												true = mean(true), rec = mean(rec),
												cor = corr, slope = slope, 
												dif = dif, percentDif = percentDif, 
												stringsAsFactors = FALSE)
			
			df = rbind(df, temp)
		}
		
	}
	
	df
}
