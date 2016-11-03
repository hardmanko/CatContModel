#pragma once

#include "CCM_Util.h"

//#include "GibbsSampler.h"
//#include "UtilityFunctions.h"

namespace CatCont {

	namespace Linear {

		struct LinearConfiguration {

			ModelVariant modelVariant; //? This is duplicate information, but the likelihood function needs it...

			struct {
				double lower;
				double upper;
			} response;
			//study??

			//just make these twice the response ones to sort of constrain the mus?
			struct {
				double lower;
				double upper;
			} catMu;

		};

		//Minial version?
		struct LinearLikelihoodConfig {
			ModelVariant modelVariant;

			double lower;
			double upper;
		};

		double normalPDF(double x, double mu, double sd);
		double normalCDF(double x, double mu, double sd);
		double dtnorm_denominator(double mu, double sd, double lower, double upper);
		double dtnorm(double x, double mu, double sd, double lower, double upper, bool log = false);
		double dtnorm_noBoundCheck(double x, double mu, double sd, double lower, double upper);

		vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const LinearConfiguration& config);
		double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const LinearConfiguration& config);
	}
}