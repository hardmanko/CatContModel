#include "MF_Likelihood.h"

namespace CatCont {
namespace MemFirst {

namespace Circular {

	vector<double> zlLikelihood(const zlParameters& par, const ConditionData& data) {

		const double unifDens = 1 / (2 * PI);
		const double guessDens = (1 - par.pMem) * unifDens;

		vector<double> likelihoods(data.study.size());

		for (unsigned int i = 0; i < data.study.size(); i++) {
			double memDens = par.pMem * vmLut.dVonMises(data.response[i], data.study[i], par.contSD);

			likelihoods[i] = memDens + guessDens;

			likelihoods[i] *= (PI / 180.0); // Scale likelihood
		}
		return likelihoods;
	}

	vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config) {

		if (config.modelVariant == ModelVariant::ZL) {

			zlParameters zlPar;
			zlPar.pMem = par.pMem;
			zlPar.contSD = par.contSD;

			return zlLikelihood(zlPar, data);
		}

		bool calculateWithinComponent = (config.modelVariant == ModelVariant::BetweenAndWithin) || (config.modelVariant == ModelVariant::WithinItem);
		bool calculateBetweenComponent = (config.modelVariant == ModelVariant::BetweenAndWithin) || (config.modelVariant == ModelVariant::BetweenItem);

		double pBetween = par.pBetween;
		if (config.modelVariant == ModelVariant::BetweenItem) {
			pBetween = 1;
		}
		else if (config.modelVariant == ModelVariant::WithinItem) {
			pBetween = 0;
		}

		const double unifDens = 1 / (2 * PI);



		// Precalculate this because it doesn't depend on the data
		// This is bad variable naming. These are not SDs, they have been converted to kappas.
		//double combinedWithinKappa = Circular::combineKappas(par.pContWithin, par.contSD, par.cat.SD);

		// Precalculate withinKappa because it does not depend on the data.
		double withinKappa = 0;
		if (config.withinItem_contSD_is_withinSD) {
			// In this special case, withinKappa is based on only contSD, so contSD is withinSD
			withinKappa = par.contSD; // contSD has already been converted into precision. (Bad variable naming)

		}
		else {
			// Normallly, contSD and catSD are combined based on pContWithin.
			withinKappa = CatCont::Circular::combineKappas(par.pContWithin, par.contSD, par.cat.SD);
		}


		// When calculating the within component predicted response center, weight the 
		// categorical and continuous components of the response by pContWithin
		vector<double> cmLoc(2);
		vector<double> cmWeights(2);
		cmWeights[0] = par.pContWithin;
		cmWeights[1] = 1 - par.pContWithin;


		// Prepare category weights vector (which is recalculated for each study value)
		vector<double> weights(par.cat.mu.size());

		vector<double> likelihoods(data.study.size());

		for (unsigned int i = 0; i < data.study.size(); i++) {

			// Category weights apply to both the between and within components of the model
			CatCont::Circular::categoryWeights(data.study[i], par, weights.data());

			// Within component
			double withinDensity = 0;

			if (calculateWithinComponent) {

				if (par.cat.mu.size() == 0) {

					withinDensity = vmLut.dVonMises(data.response[i], data.study[i], par.contSD);

				}
				else {

					cmLoc[0] = data.study[i]; //The continuous component of the predicted response center

					for (unsigned int k = 0; k < weights.size(); k++) {

						cmLoc[1] = par.cat.mu[k]; //The categorical component of the center

						double predictedCenter = CatCont::Circular::circularMean(cmLoc, cmWeights);

						double thisDens = vmLut.dVonMises(data.response[i], predictedCenter, withinKappa);

						withinDensity += weights[k] * thisDens;
					}

				}
			}


			double betweenCatDensity = 0;
			double catGuessDens = 0;

			// For both between memory and guessing. This is always calculated for guessing even if between is not used.
			for (unsigned int k = 0; k < weights.size(); k++) {
				double dens = vmLut.dVonMises(data.response[i], par.cat.mu[k], par.cat.SD);

				betweenCatDensity += weights[k] * dens;
				catGuessDens += dens;
			}

			// Reweight guesses by 1/K
			if (par.cat.mu.size() > 0) {
				catGuessDens /= par.cat.mu.size();
			}

			double guessingDensity = par.pCatGuess * catGuessDens + (1 - par.pCatGuess) * unifDens;

			//Between memory component
			double betweenDensity = 0;
			if (calculateBetweenComponent) {
				double contDensity = vmLut.dVonMises(data.response[i], data.study[i], par.contSD);
				betweenDensity = (par.pContBetween * contDensity) + ((1 - par.pContBetween) * betweenCatDensity);
			}

			// Combine components
			double memoryDensity = pBetween * betweenDensity + (1 - pBetween) * withinDensity;

			double likelihood = par.pMem * memoryDensity + (1 - par.pMem) * guessingDensity;

			//The PI / 180 scales the likelihood so that it would be correct if degrees had been used instead of radians
			likelihoods[i] = likelihood * (PI / 180.0);

		}

		return likelihoods;

	}

	double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config) {
		vector<double> likelihoods = betweenAndWithinLikelihood(par, data, config);

		double llSum = 0;
		for (unsigned int i = 0; i < likelihoods.size(); i++) {
			llSum += log(likelihoods[i]);
		}

		return llSum;
	}

} // namespace Circular

namespace Linear {


	vector<double> zlLikelihood(const zlParameters& par, const ConditionData& data, const ModelConfiguration& config) {

		const CatCont::Linear::LinearConfiguration& lc = config.linearConfiguration;

		const double unifDens = 1 / (lc.response.upper - lc.response.lower);
		const double guessDens = (1 - par.pMem) * unifDens;

		vector<double> likelihoods(data.study.size());

		for (unsigned int i = 0; i < data.study.size(); i++) {

			double memDens = par.pMem * CatCont::Linear::dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

			likelihoods[i] = memDens + guessDens;
		}

		return likelihoods;
	}

	vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config) {

		const CatCont::Linear::LinearConfiguration& lc = config.linearConfiguration;

		if (config.modelVariant == ModelVariant::ZL) {
			zlParameters zlPar;
			zlPar.pMem = par.pMem;
			zlPar.contSD = par.contSD;

			return zlLikelihood(zlPar, data, config);
		}

		bool calculateWithinComponent = (config.modelVariant == ModelVariant::BetweenAndWithin) || (config.modelVariant == ModelVariant::WithinItem);
		bool calculateBetweenComponent = (config.modelVariant == ModelVariant::BetweenAndWithin) || (config.modelVariant == ModelVariant::BetweenItem);

		double pBetween = par.pBetween;
		if (config.modelVariant == ModelVariant::BetweenItem) {
			pBetween = 1;
		}
		else if (config.modelVariant == ModelVariant::WithinItem) {
			pBetween = 0;
		}

		// The width of the uniform depends on the response range
		const double unifDens = 1 / (lc.response.upper - lc.response.lower);

		// Precalculate withinSD because it does not depend on the data.
		double withinSD = 0;
		if (config.withinItem_contSD_is_withinSD) {
			// In this special case, withinSD is contSD (so contSD should be interpreted as withinSD).
			withinSD = par.contSD;
		}
		else {
			// Normallly, contSD and catSD are combined based on pContWithin.
			withinSD = CatCont::Linear::combineSDs(par.pContWithin, par.contSD, par.cat.SD);
		}

		// Prepare category weights vector (which is recalculated for each study value)
		vector<double> weights(par.cat.mu.size());
		//WeightsCalculator weightCalc(config); // TODO

		vector<double> likelihoods(data.study.size());

		for (unsigned int i = 0; i < data.study.size(); i++) {

			//Category weights apply to both the between and within components of the model
			CatCont::Linear::categoryWeights(data.study[i], par, lc, weights.data());

			//Within component
			double withinDensity = 0;

			if (calculateWithinComponent) {

				if (par.cat.mu.size() == 0) {
					
					withinDensity = CatCont::Linear::dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

				}
				else {

					for (unsigned int k = 0; k < weights.size(); k++) {

						double predictedCenter = par.pContWithin * data.study[i] + (1 - par.pContWithin) * par.cat.mu[k];

						double thisDens = CatCont::Linear::dtnorm_noBoundCheck(data.response[i], predictedCenter, withinSD, lc.response.lower, lc.response.upper);

						withinDensity += weights[k] * thisDens;
					}
				}
			}


			double betweenCatDens = 0;
			double catGuessDens = 0;

			// For both between memory and guessing. This is always calculated for guessing even if between memory is not used.
			for (unsigned int k = 0; k < weights.size(); k++) {
				double dens = CatCont::Linear::dtnorm_noBoundCheck(data.response[i], par.cat.mu[k], par.cat.SD, lc.response.lower, lc.response.upper);
				betweenCatDens += weights[k] * dens;
				catGuessDens += dens;
			}

			// Reweight guesses by 1/K
			if (par.cat.mu.size() > 0) {
				catGuessDens /= par.cat.mu.size();
			}

			double guessingDensity = par.pCatGuess * catGuessDens + (1 - par.pCatGuess) * unifDens;

			//Between memory component
			double betweenDensity = 0;
			if (calculateBetweenComponent) {
				double contDensity = CatCont::Linear::dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

				betweenDensity = (par.pContBetween * contDensity) + ((1 - par.pContBetween) * betweenCatDens);
			}

			//Combine components
			double memoryDensity = pBetween * betweenDensity + (1 - pBetween) * withinDensity;

			double likelihood = par.pMem * memoryDensity + (1 - par.pMem) * guessingDensity;

			likelihoods[i] = likelihood;

		}

		return likelihoods;

	}

	double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config) {
		vector<double> likelihoods = betweenAndWithinLikelihood(par, data, config);

		double llSum = 0;
		for (unsigned int i = 0; i < likelihoods.size(); i++) {
			llSum += log(likelihoods[i]);
		}

		return llSum;
	}

} // namespace Linear

} // namespace MemFirst




/*
TODO:

For Circular likelihood, convert parameters inside the function.
MemFirstModelRunner::combineParameters should be changed to not convert circular parameters.



LikelihoodCalculator::LikelihoodCalculator(ModelConfiguration* modCfg) :
	_modCfg(modCfg)
{

	ModelVariant modelVariant = _modCfg->modelVariant;
	DataType dataType = _modCfg->dataType;


	_calculatedUniformDensity = uniformDensity(*_modCfg);

		
	if (modelVariant == ModelVariant::ZL) {

		if (dataType == DataType::Linear) {

			const Linear::LinearConfiguration& lc = _modCfg->linearConfiguration;

			_zlMemDens = std::bind(Linear::dtnorm_noBoundCheck, 
				placeholders::_1, placeholders::_2, placeholders::_3, 
				lc.response.lower, lc.response.upper);

			//_zlMemDens = [this](double response, double study, double contSD) -> double {};

		} else if (dataType == DataType::Circular) {

			_zlMemDens = [this](double response, double study, double contSD) -> double {

				return vmLut.dVonMises(response, study, contSD);
			};

		}
	}
	else if (modelVariant == ModelVariant::BetweenItem) {



	}



}

void LikelihoodCalculator::calcLikelihoods(const ConditionData& data, const CombinedParameters& par) {

	this->likelihoods.resize(data.size());

	if (_modCfg->modelVariant == ModelVariant::ZL) {
		_zlLikelihood(data, par);
	} else {
		_bwLikelihood(data, par);
	}

}

double LikelihoodCalculator::sumLikelihoods(bool log) {

	double sum = 0;

	if (log) {
		for (size_t i = 0; i < this->likelihoods.size(); i++) {
			sum += log(this->likelihoods[i]);
		}
	} else {
		for (size_t i = 0; i < this->likelihoods.size(); i++) {
			sum += this->likelihoods[i];
		}
	}

	return sum;
}


void LikelihoodCalculator::_zlLikelihood(const ConditionData& data, const zlParameters& par) {

	//const double unifDens = 1 / (lc.response.upper - lc.response.lower);

	const double guessDens = (1 - par.pMem) * _calculatedUniformDensity;

	//this->likelihoods.resize(data.study.size());

	for (size_t i = 0; i < data.study.size(); i++) {

		double memDens = par.pMem * _zlMemDens(data.response[i], data.study[i], par.contSD);

		likelihoods[i] = memDens + guessDens;
	}

}

void LikelihoodCalculator::_bwLikelihood(const ConditionData& data, const CombinedParameters& par) {

	WeightsCalculator weightsCalc(_modCfg);
		



}


// ----------------------------------------------------------------------------------------
// Linear
//
void LikelihoodCalculator::_bwLikelihood_linear(const ConditionData& data, const CombinedParameters& par) {
	

	bool calculateWithinComponent = (config.modelVariant == ModelVariant::BetweenAndWithin) || (config.modelVariant == ModelVariant::WithinItem);
	bool calculateBetweenComponent = (config.modelVariant == ModelVariant::BetweenAndWithin) || (config.modelVariant == ModelVariant::BetweenItem);

	double pBetween = par.pBetween;
	if (config.modelVariant == ModelVariant::BetweenItem) {
		pBetween = 1;
	} else if (config.modelVariant == ModelVariant::WithinItem) {
		pBetween = 0;
	}

	//The width of the uniform depends on the response range
	//const double unifDens = 1 / (lc.response.upper - lc.response.lower);



	//vector<double> weights(par.cat.mu.size());
	WeightsCalculator weightCalc(this->_modCfg);
	size_t nCat = par.cat.nCat();

	//Precalculate this because it doesn't depend on the data
	double combinedWithinSd = Linear::combineSDs(par.pContWithin, par.contSD, par.cat.SD);

	//this->likelihoods.resize(data.study.size());

	const LinearConfiguration& lc = config.linearConfiguration;

	for (size_t i = 0; i < data.study.size(); i++) {

		//Category weights apply to both the between and within components of the model
		//categoryWeights(data.study[i], par, lc, weights.data()); // Old calculation
		weightCalc.calcWeights(data.study[i], par.cat);


		// TODO: Given the calculated weights, calculate new pMem and pContBetween based on weight sum
		//if (_modCfg->weightsDistribution == WeightsDistribution::PlatSpline) {
			//double weightLoss = 1 - weightCalc.sumWeights();
		//}

		//Within component
		double withinDensity = 0;

		if (calculateWithinComponent) {

			if (nCat == 0) {

				// If there are no categories, act as though it were fully continuous. It's not clear that this is right.
				withinDensity = dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

			} else {

				for (size_t j = 0; j < nCat; j++) {

					double predictedCenter = par.pContWithin * data.study[i] + (1 - par.pContWithin) * par.cat.mu[j];

					double thisDens = dtnorm_noBoundCheck(data.response[i], predictedCenter, combinedWithinSd, lc.response.lower, lc.response.upper);

					withinDensity += weightCalc.weights[j] * thisDens;
				}
			}
		}


		double betweenCatDens = 0;
		double catGuessDens = 0;

		//For both between and guessing. This is always calculated for guessing even if between isn't calculated.
		for (size_t j = 0; j < nCat; j++) {
			double dens = dtnorm_noBoundCheck(data.response[i], par.cat.mu[j], par.cat.SD, lc.response.lower, lc.response.upper);
			betweenCatDens += weightCalc.weights[j] * dens;
			catGuessDens += dens;
		}

		//Between component
		double betweenDensity = 0;
		if (calculateBetweenComponent) {
			double contDensity = dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

			betweenDensity = (par.pContBetween * contDensity) + ((1 - par.pContBetween) * betweenCatDens);
		}


		//Guessing component
		if (par.cat.mu.size() > 0) {
			catGuessDens /= par.cat.mu.size();
		}

		double guessingDensity = par.pCatGuess * catGuessDens + (1 - par.pCatGuess) * _calculatedUniformDensity;


		//Combine components
		double memoryDensity = pBetween * betweenDensity + (1 - pBetween) * withinDensity;

		double likelihood = par.pMem * memoryDensity + (1 - par.pMem) * guessingDensity;

		this->likelihoods[i] = likelihood;

	}

	//return likelihoods;
}

// ----------------------------------------------------------------------------------------
// Circular
//
void LikelihoodCalculator::_bwLikelihood_circular(const ConditionData& data, const CombinedParameters& par) {

	bool calculateWithinComponent = (modelVariant == ModelVariant::BetweenAndWithin) || (modelVariant == ModelVariant::WithinItem);
	bool calculateBetweenComponent = (modelVariant == ModelVariant::BetweenAndWithin) || (modelVariant == ModelVariant::BetweenItem);

	double pBetween = par.pBetween;
	if (modelVariant == ModelVariant::BetweenItem) {
		pBetween = 1;
	} else if (modelVariant == ModelVariant::WithinItem) {
		pBetween = 0;
	}

	vector<double> weights(par.cat.mu.size());

	//Precalculate this because it doesn't depend on the data
	// This is bad variable naming. These are not SDs, they have been converted to kappas.
	double combinedWithinKappa = Circular::combineKappas(par.pContWithin, par.contSD, par.cat.SD);


	//When calculating the within component predicted response center, weight the 
	//categorical and continuous components of the response by pContWithin
	vector<double> cmLoc(2);
	vector<double> cmWeights(2);
	cmWeights[0] = par.pContWithin;
	cmWeights[1] = 1 - par.pContWithin;


	vector<double> likelihoods(data.study.size());

	for (unsigned int i = 0; i < data.study.size(); i++) {

		//Category weights apply to both the between and with components of the model
		categoryWeights(data.study[i], par, weights.data());

		//Within component
		double withinDensity = 0;

		if (calculateWithinComponent) {

			if (par.cat.mu.size() == 0) {

				withinDensity = vmLut.dVonMises(data.response[i], data.study[i], par.contSD);

			} else {

				cmLoc[0] = data.study[i]; //The continuous component of the predicted response center

				for (unsigned int j = 0; j < weights.size(); j++) {

					cmLoc[1] = par.cat.mu[j]; //The categorical component of the center

					double predictedCenter = circularMean(cmLoc, cmWeights);

					double thisDens = vmLut.dVonMises(data.response[i], predictedCenter, combinedWithinKappa);

					withinDensity += weights[j] * thisDens;
				}

			}
		}


		double betweenCatDensity = 0;
		double catGuessDens = 0;

		//For both between and guessing. This is always calculated for guessing even if between isn't calculated.
		for (unsigned int j = 0; j < weights.size(); j++) {

			double dens = vmLut.dVonMises(data.response[i], par.cat.mu[j], par.cat.SD);

			betweenCatDensity += weights[j] * dens;
			catGuessDens += dens;
		}

		//Between component
		double betweenDensity = 0;

		if (calculateBetweenComponent) {
			double contDensity = vmLut.dVonMises(data.response[i], data.study[i], par.contSD);
			betweenDensity = (par.pContBetween * contDensity) + ((1 - par.pContBetween) * betweenCatDensity);
		}


		//Guessing component (sort of). catGuessDens is calculated above.
		if (par.cat.mu.size() > 0) {
			catGuessDens /= par.cat.mu.size();
		}

		double guessingDensity = par.pCatGuess * catGuessDens + (1 - par.pCatGuess) * _calculatedUniformDensity;


		//Combine components
		double memoryDensity = pBetween * betweenDensity + (1 - pBetween) * withinDensity;

		double likelihood = par.pMem * memoryDensity + (1 - par.pMem) * guessingDensity;

		//likelihoods[i] = likelihood * (PI / 180.0);

	}

	return likelihoods;

}

double LikelihoodCalculator::_bwLikelihood_within(const ConditionData& data, const CombinedParameters& par) {

}

double LikelihoodCalculator::_bwLikelihood_between(const ConditionData& data, const CombinedParameters& par) {

}

*/


} // namespace CatCont