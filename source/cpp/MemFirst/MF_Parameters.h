#pragma once

#include "CCM_Main.h"
#include "CCM_ModelConfig.h"
#include "CCM_Util.h"
#include "CCM_Circular.h"

#include "GibbsParameters.h"

//#define USING_CAT_SPLINE

namespace CatCont {
namespace MemFirst {

struct zlParameters {
	double pMem;
	double contSD;
};

#ifdef USING_CAT_SPLINE
struct PlatSplineParameters {
	vector<double> platHW; // half width

	vector<double> splineHW; // half width (or standard deviation)

	vector<double> height;

	size_t nCat(void) const {
		if (!isValid()) {
			return 0;
		}

		return platHW.size();
	}

	bool isValid(void) const {
		return platHW.size() == splineHW.size() && platHW.size() == height.size();
	}
};
#endif


struct CategoryParameters {

	size_t nCat(void) const {
		return mu.size();
	}

	// Parameters related to categorization (the weights function)

	std::vector<double> mu; // catMu: Location of the category
	double selectivity; // catSel: How categories transition into one another

#ifdef USING_CAT_SPLINE
	vector<double> platHW; // PlatSpline plateau half width
	//double beta; // PlatSpline lambda multiplier. Not really a category parameter, but closely related to categorization.
	//PlatSplineParameters psPar;
#endif

// Parameters related to memory.

	double SD; // catSD: Memory precision for category parameters
};



struct Parameters {

	double pMem;
	double pBetween;
	double pContBetween;
	double pContWithin;

	double pCatGuess;

	double contSD;

	CategoryParameters cat;
};

struct ConditionCategoryParameters {
	double selectivity;
	double SD;

	//double beta; // lambda
};

struct ConditionParameters {
	double pMem;
	double pBetween;
	double pContBetween;
	double pContWithin;

	double pCatGuess;

	double contSD;

	ConditionCategoryParameters cat;
};

//typedef Parameters ParticipantParameters;
//typedef Parameters CombinedParameters;

struct ParticipantParameters : Parameters {};
struct CombinedParameters : Parameters {};

// Parameter extractors
ParticipantParameters getParticipantParameters(const ParameterList& param, std::string pnum, unsigned int maxCategories);
ConditionParameters getConditionParameters(const ParameterList& param, std::string condName);

ParticipantParameters getParticipantParameters(const ParamContainer& param, std::string pnum, unsigned int maxCategories);
ConditionParameters getConditionParameters(const ParamContainer& param, std::string condName);

CombinedParameters combineParameters(const ParticipantParameters& part, const ConditionParameters& cond, const SDRanges& sdRanges, DataType dataType);


} // namespace MemFirst
} // namespace CatCont