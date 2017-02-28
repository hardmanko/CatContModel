#pragma once

#include "CCM_Util.h"

namespace CatCont {

	namespace Linear {

		struct LinearConfiguration {

			ModelVariant modelVariant; //? This is duplicate information, but the likelihood function needs it...

			struct {
				double lower;
				double upper;
			} response;

			//Study ranges are not needed

			struct {
				double lower;
				double upper;
			} catMu;

		};

		double normalPDF(double x, double mu, double sd);
		double normalCDF(double x, double mu, double sd);
		double dtnorm_denominator(double mu, double sd, double lower, double upper);
		double dtnorm(double x, double mu, double sd, double lower, double upper, bool log_ = false);
		vector<double> dtnorm(const vector<double>& x, double mu, double sd, double lower, double upper, bool log_);
		double dtnorm_noBoundCheck(double x, double mu, double sd, double lower, double upper);

		vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const LinearConfiguration& config);
		double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const LinearConfiguration& config);
	}
}