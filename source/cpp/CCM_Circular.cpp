#include "CCM_Circular.h"

namespace CatCont {
namespace Circular {

double degreesToRadians(double degrees) {
	return degrees * PI / 180.0;
}
double radiansToDegrees(double radians) {
	return radians * 180.0 / PI;
}

vector<double> degreesToRadians(const vector<double>& degrees) {
	vector<double> radians(degrees.size());
	for (unsigned int i = 0; i < degrees.size(); i++) {
		radians[i] = degrees[i] * PI / 180.0;
	}
	return radians;
}

vector<double> radiansToDegrees(const vector<double>& radians) {
	vector<double> degrees(radians.size());
	for (unsigned int i = 0; i < radians.size(); i++) {
		degrees[i] = radians[i] * 180.0 / PI;
	}
	return degrees;
}

double precRad_to_sdDeg(double precRad) {
	double sdRad = 1.0 / sqrt(precRad); //precision to sd
	return sdRad * 180.0 / PI; //rad to deg
}

double sdDeg_to_precRad(double sdDeg) {
	double sdRad = sdDeg * PI / 180.0; //degrees to rad
	return 1.0 / (sdRad * sdRad); //sd to precision
}



double circularMean(double rad1, double rad2) {
	vector<double> rads = { rad1, rad2 };
	//rads[0] = rad1;
	//rads[1] = rad2;
	return circularMean(rads);

	// or
	//circularMean(0.5, rad1, rad2);
}

double circularMean(const vector<double>& radians) {
	vector<double> weights(radians.size(), 1.0 / radians.size());
	return circularMean(radians, weights);
}

// weights should sum to 1
double circularMean(const vector<double>& radians, const vector<double>& weights) {
	double cosx = 0;
	double sinx = 0;
	for (size_t i = 0; i < radians.size(); i++) {
		cosx += cos(radians[i]) * weights[i];
		sinx += sin(radians[i]) * weights[i];
	}
	double atn = atan2(sinx, cosx);
	if (atn < 0) {
		atn += 2 * PI;
	}
	return atn;
}

// Special case for what is actually used in the likelihood
// p1 is probability for rad1, (1 - p1) is probability for rad2
double circularMean(double p1, double rad1, double rad2) {

	double cosx = p1 * cos(rad1) + (1 - p1) * cos(rad2);
	double sinx = p1 * sin(rad1) + (1 - p1) * sin(rad2);

	double atn = atan2(sinx, cosx);
	if (atn < 0) { // while?
		atn += 2 * PI;
	}
	return atn;

}


// 0 <= rval < 2pi
double clampAngle360(double x) {

	x = std::fmod(x, 2*PI);

	if (x < 0) {
		x += 2*PI;
	}

	return x;
}

// -pi <= rval < pi
double clampAngle180(double x) {

	x = std::fmod(x, 2 * PI);

	if (x < -PI) {
		x += 2 * PI;
	} else if (x >= PI) {
		x -= 2 * PI;
	}

	return x;
}

vector<double> clampAngle180(const vector<double>& xs) {
	vector<double> rval(xs.size());
	for (size_t i = 0; i < xs.size(); i++) {
		rval[i] = clampAngle180(xs[i]);
	}
	return rval;
}

vector<double> clampAngle360(const vector<double>& xs) {
	vector<double> rval(xs.size());
	for (size_t i = 0; i < xs.size(); i++) {
		rval[i] = clampAngle360(xs[i]);
	}
	return rval;
}

double clampAngle(double x, bool pm180, bool degrees) {

	if (degrees) {
		x = degreesToRadians(x);
	}

	double rval = pm180 ? clampAngle180(x) : clampAngle360(x);

	if (degrees) {
		rval = radiansToDegrees(rval);
	}

	return rval;

	/*

	double fcirc = degrees ? 360 : 2 * PI;
	double hcirc = degrees ? 180 : PI;

	x = std::fmod(x, fcirc);

	if (pm180) {
		// -180 <= x <= 180
		if (x <= -hcirc) {
			x += fcirc;
		} else if (x > hcirc) {
			x -= fcirc;
		}
	} else {
		// 0 <= x <= 360
		if (x < 0) {
			x += fcirc;
		}
	}

	return x;
	*/
}


// assume radians
double circularSignedDistance(double x, double y) {
	return clampAngle180(y - x);
}

double circularAbsoluteDistance(double x, double y) {
	return std::abs(circularSignedDistance(x, y));
}

double circularDistance(double x, double y, bool absDist, bool degrees) {
	if (degrees) {
		x = degreesToRadians(x);
		y = degreesToRadians(y);
	}

	double d = circularSignedDistance(x, y);
	if (absDist) {
		d = std::abs(d);
	}

	if (degrees) {
		d = radiansToDegrees(d);
	}

	return d;
}





// See Equation 25 in the Appendix. This is a rewritten form. sigma^2 = 1 / kappa
double combineKappas(double pContWithin, double contKappa, double catKappa) {

	// This is based on the idea that the study and category locations both have some amount of
	// noise added to them, then are averaged. Thus, the average has reduced variance.
	double combinedVariance = (pow(pContWithin, 2) / contKappa) + (pow((1 - pContWithin), 2) / catKappa);

	double newKappa = 1 / combinedVariance; //convert variance back to precision

	// Make sure the newKappa is in range. This is invisible to the rest of the system!
	// However, it is also very rare that this is needed, at least for my data.
	newKappa = std::min(newKappa, vmLut.getMaxKappa());

	return newKappa;
}



} // namespace Circular
} // namespace CatCont