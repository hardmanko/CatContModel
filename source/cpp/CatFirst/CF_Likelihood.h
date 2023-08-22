#pragma once

#include "CCM_Main.h"
#include "CCM_Data.h"
#include "CCM_Util.h"

#include "CCM_Circular.h"
#include "CCM_Linear.h"
#include "CCM_Weights.h"

#include "CF_Parameters.h"
#include "CF_ModelConfig.h"


namespace CatCont {
namespace CatFirst {

class LikelihoodCalculator {
public:

	bool setup(CF_CompositeConfig config);

	std::vector<double> likelihood(const CatCont::ConditionData& data, const MPMP_complete& par) const;

private:

	CF_CompositeConfig _config;
	
	// If binding
	std::function<double(double, double, double)> _mainLikeImpl;
	std::function<double(double, double, double)> _withinCenterPredImpl;

	
	//double _uniformDensity(void) const;

	void _setupLikelihoods(void);



	// Is it faster to bind than to branch?
	double _mainLikelihood(double response, double center, double variability) const;


	double _withinCenterPrediction(double pContWithin, double studyLoc, double catLoc) const;
	double _within_combineVariability(double pContWithin, double contVariability, double catVariability) const;

	// Move this to Linear?
	static double _withinCenterPrediction_linear(double pContWithin, double contLoc, double catLoc);

	
	vector<double> _like_BW(const ConditionData& data, const MPMP_complete& par) const;

};



} // namespace CatFirst
} // namespace CatCont