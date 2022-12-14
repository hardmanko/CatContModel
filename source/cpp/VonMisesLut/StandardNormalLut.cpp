#include "StandardNormalLut.h"

void StandardNormalLUT::setup(const StandardNormalLUT::Config& config) {
	_config = config;

	_config.sdRange = std::min(_config.sdRange, 39.0);

	if (_config.stepSize > 0.1) {
		// warn?
	}

	size_t length = config.sdRange / config.stepSize;
	length += 10; // padding

	_dnormLUT.resize(length);
	_pnormLUT.resize(length);

	for (size_t i = 0; i < _dnormLUT.size(); i++) {

		double x = config.stepSize * i;

		_dnormLUT[i] = config.dnormFun(x);
		_pnormLUT[i] = config.pnormFun(x);
	}
}


const StandardNormalLUT::Config& StandardNormalLUT::getConfig(void) {
	return _config;
}

  ///////////
 // dnorm //
///////////


double StandardNormalLUT::dnorm(double x) {
	// The normal is symmetrical, so the sign of x doesn't matter
	x = std::abs(x);

	// If x is out of range, assume zero density
	if (x > _config.sdRange) {
		return 0.0;
	}

	x /= stepSize; // This helps with precision vs doing just fmod(x, stepSize).

	double rem = fmod(x, 1);
	size_t i = (size_t)x;

	double lutI = _dnormLUT[i];

	return lutI + rem * (_dnormLUT[i + 1] - lutI);
}

double StandardNormalLUT::dnorm(double x, double mu, double sigma) {
	return this->dnorm((x - mu) / sigma);
}


double StandardNormalLUT::dnorm(double x, double mu, double sigma, bool log) {

	double rval = this->dnorm(x, mu, sigma);

	if (log) {
		rval = std::log(rval);
	}

	return rval;
}

  ///////////
 // pnorm //
///////////

double StandardNormalLUT::pnorm(double x) {
	
	bool neg = x < 0;

	// Once the sign of x is noted, strip the sign to calculate index
	x = std::abs(x);

	// If x is out of range, return 0 if x < 0 and 1 if x > 0
	if (x > _config.sdRange) {
		return neg ? 0.0 : 1.0;
	}

	x /= stepSize; // This helps with precision vs doing just fmod(x, stepSize).

	double rem = fmod(x, 1);
	size_t i = (size_t)x;

	double lutI = _pnormLUT[i];

	double lerpedValue = lutI + rem * (_pnormLUT[i + 1] - lutI);

	// If x < 0, then the integral should be flipped.
	if (neg) {
		lerpedValue = 1.0 - lerpedValue;
	}

	return lerpedValue;
}

double StandardNormalLUT::pnorm(double x, double mu, double sigma) {
	return this->pnorm((x - mu) / sigma);
}

