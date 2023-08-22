#pragma once

#include "CCM_Main.h"

#include "CCM_DistributionLUTs.h"

namespace CatCont {
namespace Circular {

	double degreesToRadians(double degrees);
	double radiansToDegrees(double radians);
	vector<double> degreesToRadians(const vector<double>& degrees);
	vector<double> radiansToDegrees(const vector<double>& radians);

	double sdDeg_to_precRad(double sdDeg);
	double precRad_to_sdDeg(double precRad);

	// For within-item
	double combineKappas(double pContWithin, double contKappa, double catKappa);

	double circularMean(const vector<double>& radians, const vector<double>& weights);
	double circularMean(const vector<double>& radians);
	double circularMean(double rad1, double rad2);
	double circularMean(double p1, double rad1, double rad2);

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

} // namespace Circular
} // namespace CatCont