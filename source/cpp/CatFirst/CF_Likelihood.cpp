#include "CF_Likelihood.h"

#include "CCM_Weights.h"
#include "MF_ModelUtil.h" // uniformDensity

#include "R_Interface.h" // readModelConfig


// [[Rcpp::export]]
Rcpp::NumericMatrix CCM_CF_likelihoodWrapper(Rcpp::List modelConfigList, Rcpp::DataFrame dataDF, Rcpp::NumericMatrix paramMat) {
	using namespace Rcpp;
	using namespace CatCont;

	// Calculate likelihood for each parameter iteration, for each row of the data. 
	// It's most efficient to process chunks of data per part/cond pair.

	CatFirst::CF_CompositeConfig cfConfig;
	cfConfig.setFromList(modelConfigList);

	//CatCont::ModelConfiguration modelConfig;
	//modelConfig.setFromList(modelConfigList);

	CatCont::DataCollection data;
	data.setup(dataDF, cfConfig.modelConfig.dataType, false);

	CatCont::CatFirst::LikelihoodCalculator likeCalc;
	likeCalc.setup(cfConfig);

	CharacterVector paramNamesCV = colnames(paramMat);
	std::vector<std::string> paramNames(paramNamesCV.begin(), paramNamesCV.end());

	ParamContainer paramCont;
	paramCont.setup(paramNames);

	CatCont::CatFirst::ParamExtractor paramExtractor;
	paramExtractor.setup(cfConfig, data, paramCont);

	// 1 likelihood for each observation and each iteration.
	NumericMatrix rval(dataDF.nrow(), paramMat.nrow());

	for (size_t iter = 0; iter < (size_t)paramMat.nrow(); iter++) {
		NumericVector iterParamNV = paramMat(iter, _);
		std::vector<double> iterParam(iterParamNV.begin(), iterParamNV.end());

		paramCont.setValues(iterParam);

		for (size_t i = 0; i < data.participants.size(); i++) {
			for (size_t j = 0; j < data.conditionNames.size(); j++) {
				const ConditionData& condData = data.participants[i].condData[j];
				if (condData.study.size() == 0) {
					continue;
				}

				CatCont::CatFirst::MPMP_complete completePar = paramExtractor.getCompleteMPMP(paramCont, i, j);

				//CatCont::CatFirst::MPMP_complete completePar;
				//completePar.setFromContainer(modelConfig, paramExtractor, paramCont, i, j);

				// DEBUG
				CatCont::message("Parameter values: \n" + completePar.printValues());

				vector<double> likes = likeCalc.likelihood(condData, completePar);
				
				for (size_t obs = 0; obs < condData.originalRow.size(); obs++) {
					rval(condData.originalRow[obs], iter) = likes[obs];
				}
				
			}
		}

	}

	return rval;
}



namespace CatCont {
namespace CatFirst {

bool LikelihoodCalculator::setup(CF_CompositeConfig config) {
	this->_config = config;

	if (_config.modelConfig.dataType == CatCont::DataType::Circular) {
		conditionalConfigureVMLut(_config.modelConfig.sdRanges.maxPrecision, VON_MISES_STEP_SIZE, false); // TODO: runConfig.verbose?
	}

	this->_setupLikelihoods();

	return true;
}

std::vector<double> LikelihoodCalculator::likelihood(const CatCont::ConditionData& data, const MPMP_complete& par) const {
	return this->_like_BW(data, par);
}


// If binding

void LikelihoodCalculator::_setupLikelihoods(void) {

	using namespace std::placeholders;

	_mainLikeImpl = nullptr;
	_withinCenterPredImpl = nullptr;

	if (_config.modelConfig.dataType == DataType::Circular) {

		// This was the first attempt, which fails because there are two dVonMises functions with different args.
		//_mainLikeImpl = std::bind(&CatCont::VonMisesLUT::dVonMises, &CatCont::vmLut, _1, _2, _3);

		// 1. One solution is a cast to specify the type signature of the function.
		//_mainLikeImpl = std::bind(static_cast<double(*)(CatCont::VonMisesLUT*,double,double,double)>(&CatCont::VonMisesLUT::dVonMises), &CatCont::vmLut, _1, _2, _3);
		
		// 2. Or use a lambda.
		_mainLikeImpl = [](double x, double mu, double kappa) -> double {
			return CatCont::vmLut.dVonMises(x, mu, kappa);
		};

		// There is a similar problem with circularMean, which has multiple versions.
		// Use a cast to disambiguate the type of overloaded function circularMean.
		_withinCenterPredImpl = static_cast<double(*)(double, double, double)>(&CatCont::Circular::circularMean);

	}
	else if (_config.modelConfig.dataType == DataType::Linear) {
		double linLower = _config.modelConfig.linearConfiguration.response.lower;
		double linUpper = _config.modelConfig.linearConfiguration.response.upper;

		_mainLikeImpl = std::bind(CatCont::Linear::dtnorm_noBoundCheck, _1, _2, _3, linLower, linUpper);

		_withinCenterPredImpl = &LikelihoodCalculator::_withinCenterPrediction_linear;
	}

	if (!_mainLikeImpl || !_withinCenterPredImpl) {
		// Error
	}

}




// Is it faster to bind than to branch?
double LikelihoodCalculator::_mainLikelihood(double response, double center, double variability) const {
	// Option 1 : Bind
	return this->_mainLikeImpl(response, center, variability);

	// Option 2: Branch
	/*
	if (_config.modelConfig.dataType == DataType::Circular) {
		return CatCont::vmLut.dVonMises(response, center, variability);
	}
	else if (_config.modelConfig.dataType == DataType::Linear) {
		const double& linLower = _config.modelConfig.linearConfiguration.response.lower;
		const double& linUpper = _config.modelConfig.linearConfiguration.response.upper;

		return CatCont::Linear::dtnorm_noBoundCheck(response, center, variability, linLower, linUpper);
	}

	return -1;
	*/
}


double LikelihoodCalculator::_withinCenterPrediction(double pContWithin, double studyLoc, double catLoc) const {

	// Option 1: Bind
	return this->_withinCenterPredImpl(pContWithin, studyLoc, catLoc);

	// Option 2: Branch
	/*
	double rval = -1;
	if (_config.modelConfig.dataType == DataType::Circular) {
		rval = CatCont::Circular::circularMean(pContWithin, studyLoc, catLoc);
	}
	else if (_config.modelConfig.dataType == DataType::Linear) {
		rval = pContWithin * studyLoc + (1 - pContWithin) * catLoc;
	}

	return rval;
	*/
}


// Move this to Linear?
double LikelihoodCalculator::_withinCenterPrediction_linear(double pContWithin, double contLoc, double catLoc) {
	return pContWithin * contLoc + (1 - pContWithin) * catLoc;
}

double LikelihoodCalculator::_within_combineVariability(double pContWithin, double contVariability, double catVariability) const {
	double rval = -1;
	if (_config.modelConfig.dataType == DataType::Circular) {
		rval = Circular::combineKappas(pContWithin, contVariability, catVariability);
	}
	else if (_config.modelConfig.dataType == DataType::Linear) {
		rval = Linear::combineSDs(pContWithin, contVariability, catVariability);
	}

	return rval;
}



vector<double> LikelihoodCalculator::_like_BW(const ConditionData& data, const MPMP_complete& par) const {
	
	// Parameter constraints for reduced models are currently done during parameter extraction. Should they be done here instead?

	// maybe?
	bool mv_BW = _config.modelConfig.modelVariant == ModelVariant::BetweenAndWithin;
	bool mv_BI = _config.modelConfig.modelVariant == ModelVariant::BetweenItem;
	bool mv_WI = _config.modelConfig.modelVariant == ModelVariant::WithinItem;
	//bool mv_ZL = _config.modelConfig.modelVariant == ModelVariant::ZL; // This variable is not used

	/*
	bool calcWithinMem = (modelConfig.modelVariant == ModelVariant::BetweenAndWithin) || 
		(modelConfig.modelVariant == ModelVariant::WithinItem);
	bool calcBetweenMem = (modelConfig.modelVariant == ModelVariant::BetweenAndWithin) || 
		(modelConfig.modelVariant == ModelVariant::BetweenItem) || 
		(modelConfig.modelVariant == ModelVariant::ZL);
	*/


	const double unifDens = CatCont::MemFirst::uniformDensity(_config.modelConfig);

	//const CategoryParam& catPar = par.categorization.cats; // alias
	//const size_t nCatActive = catPar.nCatActive();


	//const vector<double> activeCatMu = par.categorization.cats.getActiveCatMu();
	//const vector<size_t> activeCatType = par.categorization.cats.get
	const ActiveCategoryParam activeCatPar = par.categorization.cats.getActiveParam();
	const size_t nCatActive = activeCatPar.catMu.size();

	// Precalculate within variability because they do not depend on the data.
	// Only used for models with a within component.
	vector<double> withinVariability_ak(nCatActive);
	if (mv_BW || mv_WI) { 
		for (size_t ak = 0; ak < nCatActive; ak++) {

			size_t catTypeIndex = activeCatPar.catType.at(ak);
			const MPMP_memory& memPar = par.memory.at(catTypeIndex);

			if (_config.modelConfig.withinItem_contSD_is_withinSD) {
				// In this special case, withinKappa is based on only contSD, so contSD is withinSD
				withinVariability_ak[ak] = memPar.contSD; // For circular, contSD has already been converted into precision. (Bad variable naming)
			}
			else {
				// Normallly, contSD and catSD are combined based on pContWithin.
				withinVariability_ak[ak] = this->_within_combineVariability(memPar.pContWithin, memPar.contSD, memPar.catSD);
			}
		}
	}


	// Calculate the probability (weight) that each category will be used when making a catGuess.
	// This depends on using the bizarro pCatUsedToGuess parameter.
	vector<double> catGuess_catWeight(nCatActive);
	double catGuess_catWeight_sum = 0;
	for (size_t ak = 0; ak < nCatActive; ak++) {
		catGuess_catWeight[ak] = par.guessing.pCatUsedToGuess[ak] / nCatActive;
		catGuess_catWeight_sum += catGuess_catWeight[ak];
	}
	// Normalize catGuess_catWeight
	for (size_t ak = 0; ak < nCatActive; ak++) {
		catGuess_catWeight[ak] /= catGuess_catWeight_sum;
	}
	//vector<double> catGuess_catWeight(nCatActive, 1.0 / nCatActive); // DEBUG
	

	// Prepare category weights vector (which is recalculated for each study value)
	WeightsCalculator weightCalc(_config.modelConfig);
	weightCalc.config.zapsmallCutoff = 1e-6;


	// TODO: Finish implementing VP
	/*
	// Sample variable precisions for each observation per category
	vector<vector<double>> contSD_vpSamples(nCatActive); // indices: k, obs
	
	for (size_t ak = 0; ak < nCatActive; ak++) {
		const MPMP_memory& memPar = par.memory[ak];
		if (_config.vpConfig.usingVP) {
			contSD_vpSamples[ak] = _config.vpConfig.sampleVariableContSD(memPar.contSD, data.study.size());
		} else {
			// If not VP, just duplicate the value
			contSD_vpSamples[ak].resize(data.study.size(), memPar.contSD);
		}
	}
	*/
	// TODO: Anywhere that contSD would be used, use contSD_vpSamples[ak][obs] instead.
	// Including within mix?


	vector<double> likelihoods(data.study.size());

	for (size_t obs = 0; obs < data.study.size(); obs++) {

		// Special case here for 0 active cats (which includes ZL variant).
		if (nCatActive == 0) {
			const MPMP_memory& memPar = par.memory.front();

			double memoryDensity = _mainLikelihood(data.response[obs], data.study[obs], memPar.contSD);

			likelihoods[obs] = memPar.pMem * memoryDensity + (1 - memPar.pMem) * unifDens;

			continue;
		}
		// Below here, there is always at least one active category (not ZL variant).

		// Calculate guessing component of the model, which does not depend on memory.
		// This is the density for the guessing node, including both categorical and uniform guessing.
		double guessDens_weightedSum = 0;
		for (size_t ak = 0; ak < nCatActive; ak++) {



			double catGuessDens = _mainLikelihood(data.response[obs], activeCatPar.catMu[ak], par.guessing.catSD[ak]);

			double guessDens = par.guessing.pCatGuess * catGuessDens + (1 - par.guessing.pCatGuess) * unifDens;

			double weightedGuessDens = catGuess_catWeight[ak] * guessDens;

			guessDens_weightedSum += weightedGuessDens;
		}
		// No special case for 0 active cats or ZL because that is handled at the top of the loop.
		//if (mv_ZL || nCatActive == 0) {
		//	guessDens_weightedSum = unifDens;
		//}

		// Category weights apply to both the between and within components of the model
		weightCalc.calcWeights(data.study[obs], activeCatPar.catMu, par.categorization.catSelectivity);

		double fullTreeLikelihood = 0;

		for (size_t ak = 0; ak < nCatActive; ak++) {

			// If the weight for this category is 0, this branch can be skipped.
			if (weightCalc.weights[ak] == 0) {
				continue;
			}

			size_t catTypeIndex = activeCatPar.catType.at(ak);
			const MPMP_memory& memPar = par.memory.at(catTypeIndex);

			// Option 1: Calculate densities, then mix in a branching form

			/*
			// Continuous memory: For all variants but WI
			double contMemDens = 0;
			if (!mv_WI) {
				contMemDens = _mainLikelihood(data.response[obs], data.study[obs], memPar.contSD);
			}

			// The per-category between and within densities do not need to be weighted by the study weights (weightCalc) before being combined in the study tree.
			// The whole study tree is weighted.
			double betweenCatMemDens = 0;
			if (mv_BW || mv_BI) {
				betweenCatMemDens = _mainLikelihood(data.response[obs], activeCatPar.catMu[ak], memPar.catSD);
			}

			double withinMemDens = 0;
			if (mv_BW || mv_WI) {
				double predictedCenter = this->_withinCenterPrediction(memPar.pContWithin, data.study[obs], activeCatPar.catMu[ak]);

				withinMemDens = this->_mainLikelihood(data.response[obs], predictedCenter, withinVariability_ak[ak]);
			}

			double betweenDensity = memPar.pContBetween * contMemDens + (1 - memPar.pContBetween) * betweenCatMemDens;
			double memoryDensity = memPar.pBetween * betweenDensity + (1 - memPar.pBetween) * withinMemDens;

			double studyTreeDens = memPar.pMem * memoryDensity + (1 - memPar.pMem) * guessDens_weightedSum;

			fullTreeLikelihood += weightCalc.weights[ak] * studyTreeDens;
			*/


			// Option 2: Calculate prob times dens (pd*), then sum

			// Continuous memory. For all MV but WI.
			double pdContMem = 0;
			if (!mv_WI) {
				double pContMem = memPar.pMem * memPar.pBetween * memPar.pContBetween;
				double dContMem = _mainLikelihood(data.response[obs], data.study[obs], memPar.contSD);
				pdContMem = pContMem * dContMem;
			}

			// Between categorical memory. For BW and BI.
			double pdBetweenCatMem = 0;
			if (mv_BW || mv_BI) {
				double pBetweenCatMem = memPar.pMem * memPar.pBetween * (1 - memPar.pContBetween);
				double dBetweenCatMem = _mainLikelihood(data.response[obs], activeCatPar.catMu[ak], memPar.catSD);
				pdBetweenCatMem = pBetweenCatMem * dBetweenCatMem;
			}

			// Within-item mixture of cat and cont memory. For BW and WI.
			double pdWithinMem = 0;
			if (mv_BW || mv_WI) {
				double pWithinMem = memPar.pMem * (1 - memPar.pBetween);

				double predictedCenter = this->_withinCenterPrediction(memPar.pContWithin, data.study[obs], activeCatPar.catMu[ak]);

				double dWithinMem = this->_mainLikelihood(data.response[obs], predictedCenter, withinVariability_ak[ak]);

				pdWithinMem = pWithinMem * dWithinMem;
			}

			double pGuess = (1 - memPar.pMem);
			double pdGuess = pGuess * guessDens_weightedSum;

			double studyTreeDens = pdContMem + pdBetweenCatMem + pdWithinMem + pdGuess;

			fullTreeLikelihood += weightCalc.weights[ak] * studyTreeDens;
		}

		likelihoods[obs] = fullTreeLikelihood;

	}

	// Rescale likelihoods as though calculations had been done using degrees rather than radians.
	for (size_t i = 0; i < likelihoods.size(); i++) {
		likelihoods[i] *= (PI / 180.0);
	}

	return likelihoods;
}

} // namespace CatFirst
} // namespace CatCont