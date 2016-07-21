#include "VonMisesLut.h"


#ifndef PI
#define PI       3.14159265358979323846
#endif

//besselFun is a pointer to a function that is the partial bessel function of the first kind with k = 0 and a scaled exponent
void VonMisesLut::setup(double rangeUpper_, double stepSize_, double (*besselFun)(double)) {
	rangeUpper = rangeUpper_;
	stepSize = stepSize_;

	_lut.clear();

	//Avoid off by one errors/numerical imprecision when the input is at or near rangeUpper by adding steps
	double realUpperRange = rangeUpper_ + (20 * stepSize_);

	double x = 0;
	unsigned int step = 0;

	while (x <= realUpperRange) {

		double val = besselFun(x);

		_lut.push_back(val);

		//x += stepSize; //Does this maintain precision? The answer is no, it does not, but it's negligible.
		++step;
		x = step * stepSize; //This is better
	}

	//TODO: The LUT is biased because the function is monotonically decreasing.
	//Adjust by taking the midpoint between two LUT values, calculate the bessel there, then adjust the LUT
	//endpoints to something like half the distance between the midpoint lerp and midpoint bessel.
}

//x must be non-negative and must be less than rangeUpper/maxValue(). It is up to the user to verify this.
double VonMisesLut::besselLerp(double x) {
	x /= stepSize; //so this helps with precision vs doing just fmod(x, stepSize)...

	int i = (int)x;
	double rem = fmod(x, 1);

	return _lut[i] + rem * (_lut[i + 1] - _lut[i]);
}

//Calculates likelihood, not log likelihood
//x and mu are in radians, kappa is precision
double VonMisesLut::dVonMises(double x, double mu, double kappa) {

	double bessel = this->besselLerp(kappa);
	
	double num = pow(exp(cos(x - mu) - 1), kappa);

	return num / (2 * PI * bessel);
}

double VonMisesLut::dVonMises(double x, double mu, double kappa, bool log) {

	double rval = this->dVonMises(x, mu, kappa);

	if (log) {
		rval = std::log(rval);
	}

	return rval;
}

double VonMisesLut::maxValue(void) {
	return rangeUpper;
}