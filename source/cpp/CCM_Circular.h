#pragma once

#include "CCM_Util.h"

namespace CatCont {
	namespace Circular {

		double degreesToRadians(double degrees);
		double radiansToDegrees(double radians);
		vector<double> degreesToRadians(const vector<double>& degrees);
		vector<double> radiansToDegrees(const vector<double>& radians);

		double sdDeg_to_precRad(double sdDeg);
		double precRad_to_sdDeg(double precRad);

		double combineKappas(double contKappa, double catKappa, double pCont);

		double circularMean(const vector<double>& radians, const vector<double>& weights);
		double circularMean(const vector<double>& radians);
		double circularMean(double rad1, double rad2);

		double circularAbsoluteDistance(double x, double y, bool degrees);


		vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant);
		double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant);
		double betweenAndWithinNll(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant); //whatever
		vector<double> zlLikelihood(const zlParameters& par, const ConditionData& data);

		void categoryWeights(double study, const CombinedParameters& par, double* OUT_weights);

	}
}