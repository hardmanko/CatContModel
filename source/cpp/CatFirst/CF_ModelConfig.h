#pragma once

#include "CCM_Main.h"

#include "CCM_ModelConfig.h"

namespace CatCont {
namespace CatFirst {

// Some of this stuff should move to other structs.
struct CatFirstConfig {
	//CatClassMap ccm; // maybe

	map<string, vector<string>> condEffects; // TODO: Move to main config?
	//map<string, vector<string>> catEffects; // TODO: Remove
	map<string, vector<string>> catTypes;

	size_t getCatTypeCornerstoneIndex(void) const;

	//std::set<string> sharedParticipantParam; // set is easily searchable. TODO: move to main config?

	bool setFromList(Rcpp::List cl);
};

struct VariablePrecisionConfig {

	bool useVP = false;

	bool useAlpha = true; // if false, use beta
	double paramValue = std::numeric_limits<double>::signaling_NaN(); // alpha or beta
	//double beta = std::numeric_limits<double>::quiet_NaN();

	bool setFromList(Rcpp::List cl) {

		int uvpi = cl["useVP"];
		this->useVP = uvpi == 1;

		int uai = cl["useAlpha"];
		this->useAlpha = uai == 1;

		this->paramValue = cl["paramValue"];

		return true;
	}

	std::vector<double> sampleVariableContSD(double contSD, size_t nObs) const {
		std::vector<double> rval(nObs, contSD);
		if (!useVP) {
			return rval;
		}
		//mean = alpha / beta
		// contSD == mean
		double alpha;
		double beta;
		if (useAlpha) {
			alpha = this->paramValue;
			beta = alpha / contSD;
		}
		else {
			beta = this->paramValue;
			alpha = beta * contSD;
		}
		double scale = 1 / beta;
		for (size_t i = 0; i < rval.size(); i++) {
			rval[i] = R::rgamma(alpha, scale); // shape and scale
		}
		return rval;
	}

};

// maybe?
struct CF_CompositeConfig {

	bool setFromList(Rcpp::List cfgList);

	ModelConfiguration modelConfig; // main config

	CatFirstConfig cfConfig; // currently specialized things

	VariablePrecisionConfig vpConfig;

};

} // namespace CatFirst
} // namespace CatCont