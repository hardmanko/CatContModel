#pragma once

#include "CCM_Util.h"
#include "CCM_DistributionLUTs.h"

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

	// Degrees or radians
	double clampAngle(double x, bool pm180, bool degrees);
	double circularDistance(double x, double y, bool absDist, bool degrees);

	// Radians only
	double clampAngle180(double x);
	double clampAngle360(double x);
	vector<double> clampAngle180(const vector<double>& xs);
	vector<double> clampAngle360(const vector<double>& xs);
	double circularSignedDistance(double x, double y);
	double circularAbsoluteDistance(double x, double y);

	


	vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant);
	double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant);
	double betweenAndWithinNll(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant); //whatever
	vector<double> zlLikelihood(const zlParameters& par, const ConditionData& data);

	// OLD
	void categoryWeights(double study, const CombinedParameters& par, double* OUT_weights);

} // namespace Circular
} // namespace CatCont