#pragma once

#include "CCM_Util.h"

namespace CatCont {

	namespace Linear {

		double normalPDF(double x, double mu, double sd);
		double normalCDF(double x, double mu, double sd);
		double dtnorm_denominator(double mu, double sd, double lower, double upper);
		double dtnorm(double x, double mu, double sd, double lower, double upper, bool log_ = false);
		vector<double> dtnorm(const vector<double>& x, double mu, double sd, double lower, double upper, bool log_);
		double dtnorm_noBoundCheck(double x, double mu, double sd, double lower, double upper);

		vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& mCfg);
		double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& mCfg);

		// OLD
		void categoryWeights(double study, const CombinedParameters& par, const LinearConfiguration& lc, double* OUT_weights);
	}
}