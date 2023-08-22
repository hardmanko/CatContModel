#include "MF_ModelUtil.h"

#include "CCM_Circular.h"
#include "CCM_Linear.h"

namespace CatCont {
namespace MemFirst {

double uniformDensity(const ModelConfiguration& config) {

	if (config.dataType == DataType::Linear) {
		const Linear::LinearConfiguration& lc = config.linearConfiguration;

		return 1.0 / (lc.response.upper - lc.response.lower);
	}
	else if (config.dataType == DataType::Circular) {
		return 1.0 / (2.0 * PI);

	}

	// error
	return -1;
}

double catMuDeviateFunction(double currentCatMu, double candidateSD, bool clampTo360) {
	double newCatMu = CatCont::normalDeviate(currentCatMu, candidateSD);
	if (clampTo360) {
		newCatMu = CatCont::Circular::clampAngle(newCatMu, false, true); // 360, degrees
	}
	return newCatMu;
}

double catActiveDeviateFunction(double currentCatActive) {
	if (currentCatActive == 1.0) {
		return 0.0;
	}
	return 1.0;
}



bool CategoryParamPriorCalculator::setup(const ModelConfiguration& modelConfig, double distancePriorSD) {
	_modelConfig = modelConfig;

	//_catMuPriorData.sd = priors["catMuPriorSD"];
	_distancePriorData.sd = distancePriorSD;
	if (_modelConfig.dataType == DataType::Circular) {

		_distancePriorData.kappa = Circular::sdDeg_to_precRad(_distancePriorData.sd);
		_distancePriorData.maxLikelihood = vmLut.dVonMises(0, 0, _distancePriorData.kappa);

	}
	else if (_modelConfig.dataType == DataType::Linear) {

		_distancePriorData.kappa = std::numeric_limits<double>::infinity(); // kappa isn't used anywhere for linear data
		_distancePriorData.maxLikelihood = Linear::normalPDF(0, 0, _distancePriorData.sd); // NOT truncated.

	}

	return true;
}

//This is the f function from HVR17
double CategoryParamPriorCalculator::scaledDensity(const vector<double>& catMus, const vector<unsigned int>& catActives, unsigned int catIndex, unsigned int steps) const {

	if (_modelConfig.dataType == DataType::Linear) {
		// The range of catMu defines a uniform distribution that is multiplied by the rest of the prior.
		// If mu[catIndex] is out of range, it has 0 density.
		// This is where the catMu range is applied to the catMu parameters.
		if (catMus[catIndex] < _modelConfig.linearConfiguration.catMu.lower || catMus[catIndex] > _modelConfig.linearConfiguration.catMu.upper) {
			return 0;
		}
	}

	// If this category is inactive, the density is just the height of a uniform distribution.
	if (catActives[catIndex] == 0) {
		return uniformDensity(this->_modelConfig);
	}

	double unscaled = unscaledDensity(catMus, catActives, catIndex);
	double scaleFactor = estimateScaleFactor(catMus, catActives, catIndex, steps);

	return unscaled * scaleFactor;
}

//This is a combination of the g and h functions from HVR17
double CategoryParamPriorCalculator::unscaledDensity(const vector<double>& catMus, const vector<unsigned int>& catActives, unsigned int catIndex) const {

	double density = 1;

	for (unsigned int k = 0; k < catMus.size(); k++) {

		if (k != catIndex) {
			if (catActives[catIndex] == 1 && catActives[k] == 1) {

				double like;
				if (_modelConfig.dataType == DataType::Circular) {
					like = vmLut.dVonMises(catMus[catIndex], catMus[k], this->_distancePriorData.kappa);
				}
				else {
					like = Linear::normalPDF(catMus[catIndex], catMus[k], _distancePriorData.sd); //NOT truncated
				}

				double ratio = 1 - like / _distancePriorData.maxLikelihood;

				density *= (ratio * ratio); //square the ratio to make flatter bottoms on the notches
			}
		}
	}

	return density; //return non-log density

}

// This is an approximation of the S function from HVR17
// catMus are copied not referenced because they are modified
double CategoryParamPriorCalculator::estimateScaleFactor(vector<double> catMus, const vector<unsigned int>& catActives, unsigned int catIndex, unsigned int steps) const {
	// default to circular
	double stepSize = 2 * PI / steps;
	double startPoint = 0;

	if (_modelConfig.dataType == DataType::Linear) {
		stepSize = (_modelConfig.linearConfiguration.catMu.upper - _modelConfig.linearConfiguration.catMu.lower) / steps;
		startPoint = _modelConfig.linearConfiguration.catMu.lower;
	}

	double densSum = 0;

	for (unsigned int i = 0; i < steps; i++) {
		catMus[catIndex] = startPoint + i * stepSize;

		densSum += unscaledDensity(catMus, catActives, catIndex);
	}

	return 1.0 / (densSum * stepSize);
}


double CategoryParamPriorCalculator::distancePrior_catMu(const std::vector<double>& catMus, const std::vector<unsigned int>& catActives, size_t thisCatIndex, bool logDensity) const {
	double density = this->scaledDensity(catMus, catActives, thisCatIndex, _modelConfig.catMuPriorApproximationPrecision);

	if (logDensity) {
		density = std::log(density);
	}

	return density;
}

double CategoryParamPriorCalculator::distancePrior_catActive(const std::vector<double>& catMus, const std::vector<unsigned int>& catActives, size_t thisCatIndex, bool logDensity) const {

	vector<unsigned int> activeCatActives = catActives;
	vector<unsigned int> inactiveCatActives = catActives;

	activeCatActives[thisCatIndex] = 1;
	inactiveCatActives[thisCatIndex] = 0;

	double activeDens = this->scaledDensity(catMus, activeCatActives, thisCatIndex, _modelConfig.catMuPriorApproximationPrecision);
	double inactiveDens = this->scaledDensity(catMus, inactiveCatActives, thisCatIndex, _modelConfig.catMuPriorApproximationPrecision);

	double numDens = (catActives[thisCatIndex] == 1) ? activeDens : inactiveDens;
	double denDens = activeDens + inactiveDens;

	double density = numDens / denDens;
	
	if (logDensity) {
		density = std::log(density);
	}

	return density;
}


/*
double CategoryParamPriorCalculator::distancePrior_catMu(const ParameterList& param, const std::string& pnum, unsigned int catIndex) const {
	//string catIStrBase = "[" + data.participants.at(pIndex).pnum + ",";
	string catIStrBase = "[" + pnum + ",";

	vector<double> mus(_modelConfig.maxCategories);
	vector<unsigned int> catActives(_modelConfig.maxCategories);

	for (unsigned int k = 0; k < _modelConfig.maxCategories; k++) {
		string catIStr = catIStrBase + CatCont::catIndexString(k) + "]";

		mus[k] = param.at("catMu" + catIStr);
		catActives[k] = (unsigned int)param.at("catActive" + catIStr);
	}


	double density = this->scaledDensity(mus, catActives, catIndex, _modelConfig.catMuPriorApproximationPrecision);

	return std::log(density);
}

double CategoryParamPriorCalculator::distancePrior_catActive(const ParameterList& param, const std::string& pnum, unsigned int catIndex) const {
	//string catIStrBase = "[" + data.participants.at(pIndex).pnum + ",";
	string catIStrBase = "[" + pnum + ",";

	vector<double> mus(_modelConfig.maxCategories);
	vector<unsigned int> catActives(_modelConfig.maxCategories);

	for (unsigned int k = 0; k < _modelConfig.maxCategories; k++) {
		string catIStr = catIStrBase + CatCont::catIndexString(k) + "]";

		mus[k] = param.at("catMu" + catIStr);
		catActives[k] = (unsigned int)param.at("catActive" + catIStr);
	}

	vector<unsigned int> activeCatActives = catActives;
	vector<unsigned int> inactiveCatActives = catActives;

	activeCatActives[catIndex] = 1;
	inactiveCatActives[catIndex] = 0;

	double activeDens = this->scaledDensity(mus, activeCatActives, catIndex, _modelConfig.catMuPriorApproximationPrecision);
	double inactiveDens = this->scaledDensity(mus, inactiveCatActives, catIndex, _modelConfig.catMuPriorApproximationPrecision);

	double numDens = (catActives[catIndex] == 1) ? activeDens : inactiveDens;
	double denDens = activeDens + inactiveDens;

	return std::log(numDens / denDens);
}
*/




} // namespace MemFirst
} // namespace CatCont