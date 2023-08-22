#pragma once

#include "CCM_Main.h"

#include "CCM_ModelConfig.h"
#include "CCM_Util.h"

namespace CatCont {
namespace MemFirst {

	// TODO: Does not belong here, but it's better than the last place it was in.
	double uniformDensity(const ModelConfiguration& config);

	double catMuDeviateFunction(double currentCatMu, double candidateSD, bool clampTo360);

	// "Samples" a deviate for the MH step for catActive.
	double catActiveDeviateFunction(double currentCatActive);

	class CategoryParamPriorCalculator {
	public:

		// External interface
		bool setup(const ModelConfiguration& modelConfig, double distancePriorSD);

		double distancePrior_catMu(const std::vector<double>& catMus, const std::vector<unsigned int>& catActives, size_t thisCatIndex, bool logDensity) const;
		double distancePrior_catActive(const std::vector<double>& catMus, const std::vector<unsigned int>& catActives, size_t thisCatIndex, bool logDensity) const;


		// Ummm, these do parameter extraction
		//double distancePrior_catMu(const ParameterList& param, const std::string& pnum, unsigned int catIndex) const;
		//double distancePrior_catActive(const ParameterList& param, const std::string& pnum, unsigned int catIndex) const;

		// Internal helper functions
		double scaledDensity(const vector<double>& mus, const vector<unsigned int>& catActives, unsigned int catIndex, unsigned int steps) const;
		double unscaledDensity(const vector<double>& mus, const vector<unsigned int>& catActives, unsigned int catIndex) const;
		double estimateScaleFactor(vector<double> mus, const vector<unsigned int>& catActives, unsigned int catIndex, unsigned int steps) const;



	private:

		ModelConfiguration _modelConfig; // copied

		struct {
			double sd; //standard deviation in degrees/units
			double kappa; //precision in radians
			double maxLikelihood; //non-log
		} _distancePriorData;

	};

} // namespace MemFirst
} // namespace CatCont