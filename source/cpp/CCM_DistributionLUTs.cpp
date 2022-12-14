#include "CCM_DistributionLUTs.h"

#ifndef PI
#define PI       3.14159265358979323846
#endif

namespace CatCont {

// Static instances
VonMisesLUT vmLut;
//NormalLUT normalLUT;



/////////////////
// VonMisesLUT //
/////////////////

void VonMisesLUT::setup(std::function<double(double)> besselFun) {
	Config cfg;
	cfg.besselFun = besselFun;
	cfg.useLUT = false;
	this->setup(cfg);
}

void VonMisesLUT::setup(const Config& cfg) {

	Config prevConfig = _config;

	_config = cfg;

	//if (_config.stepSize > 0.1) {
		// warn?
	//}

	if (!_config.useLUT) {
		_besselLUT.clear();

		_config.stepSize = 0;
		_config.maxKappa = std::numeric_limits<double>::max();

		_internalBessel = _config.besselFun;
		return;
	}

	// Using the LUT

	// Set the internal bessel function to be the lerp function
	_internalBessel = [this](double kappa) -> double {
		return this->besselLerp(kappa);
	};

	if (_config.skipRecomputationIfAble) {

		bool hasLut = _besselLUT.size() > 0;

		bool sufficientMaxKappa = prevConfig.maxKappa >= _config.maxKappa;
		bool sufficientPrecision = prevConfig.stepSize <= _config.stepSize;

		//bool rangeCorrect = abs(_config.maxKappa - cfg.maxKappa) < 0.00001;
		//bool stepSizeCorrect = abs(_config.stepSize - cfg.stepSize) < 0.00001;

		bool skipComputation = hasLut && sufficientMaxKappa && sufficientPrecision;
		if (skipComputation) {
			return;
		}
	}

	// Recompute the LUT
	size_t newLength = cfg.maxKappa / cfg.stepSize;
	newLength += 10; // padding

	_besselLUT.resize(newLength);

	for (size_t i = 0; i < newLength; i++) {

		double x = _config.stepSize * i;

		_besselLUT[i] = _config.besselFun(x);
	}

	// TODO: The LUT is biased because the bessel function is an exponential-type scoop shape.
	// This means that the curvature of the function always undershoots the straight line linear interpolation.
	// Adjust by taking the midpoint between two LUT values, calculate the bessel there, then adjust the LUT
	// endpoints to something like half the distance between the midpoint lerp and midpoint bessel.
}

const VonMisesLUT::Config& VonMisesLUT::getConfig(void) const {
	return _config;
}

bool VonMisesLUT::ready(double maxKappa) {
	if (_config.useLUT) {
		bool hasLut = _besselLUT.size() > 0;
		bool maxKappaBigEnough = _config.maxKappa >= maxKappa;

		return hasLut && maxKappaBigEnough;
	} else {

		bool hasBesselFun = (bool)_config.besselFun;

		return hasBesselFun;
	}

	return false;
}

// kappa must be non-negative and must be less than maxValue(). It is up to the user to verify this.
double VonMisesLUT::besselLerp(double kappa) {
	kappa /= _config.stepSize; // This helps with precision vs doing just fmod(x, stepSize).

	size_t i = (size_t)kappa;
	double rem = fmod(kappa, 1.0);

	return _besselLUT[i] + rem * (_besselLUT[i + 1] - _besselLUT[i]);
}

// Calculates likelihood, not log likelihood
// x and mu are in radians, kappa is precision
double VonMisesLUT::dVonMises(double x, double mu, double kappa) {

	double bessel = _internalBessel(kappa);

	double num = pow(exp(cos(x - mu) - 1), kappa);

	return num / (2 * PI * bessel);
}

double VonMisesLUT::dVonMises(double x, double mu, double kappa, bool log) {

	double rval = this->dVonMises(x, mu, kappa);

	if (log) {
		rval = std::log(rval);
	}

	return rval;
}

double VonMisesLUT::getMaxKappa(void) const {
	return _config.maxKappa;
}

///////////////////////
// NormalLUT //
///////////////////////

void NormalLUT::setup(const NormalLUT::Config& cfg) {
	_config = cfg;

	_config.sdRange = std::min(_config.sdRange, 39.0);

	//if (_config.stepSize > 0.1) {
		// warn?
	//}

	size_t length = _config.sdRange / _config.stepSize;
	length += 10; // padding

	_dnormLUT.resize(length);
	_pnormLUT.resize(length);

	for (size_t i = 0; i < length; i++) {

		double x = _config.stepSize * i;

		_dnormLUT[i] = _config.dnormFun(x);
		_pnormLUT[i] = _config.pnormFun(x);
	}
}


const NormalLUT::Config& NormalLUT::getConfig(void) {
	return _config;
}

///////////
// dnorm //
///////////


double NormalLUT::dnorm(double x) {
	// The normal is symmetrical, so the sign of x doesn't matter
	x = std::abs(x);

	// If x is out of range, assume zero density
	if (x > _config.sdRange) {
		return 0.0;
	}

	x /= _config.stepSize; // This helps with precision vs doing just fmod(x, stepSize).

	double rem = fmod(x, 1);
	size_t i = (size_t)x;

	double lutI = _dnormLUT[i];

	return lutI + rem * (_dnormLUT[i + 1] - lutI);
}

double NormalLUT::dnorm(double x, double mu, double sigma) {
	return this->dnorm((x - mu) / sigma);
}


double NormalLUT::dnorm(double x, double mu, double sigma, bool log) {

	double rval = this->dnorm(x, mu, sigma);

	if (log) {
		rval = std::log(rval);
	}

	return rval;
}

///////////
// pnorm //
///////////

double NormalLUT::pnorm(double x) {

	bool neg = x < 0;

	// Once the sign of x is noted, strip the sign to calculate index
	x = std::abs(x);

	// If x is out of range, return 0 if x < 0 and 1 if x > 0
	if (x > _config.sdRange) {
		return neg ? 0.0 : 1.0;
	}

	x /= _config.stepSize; // This helps with precision vs doing just fmod(x, stepSize).

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

double NormalLUT::pnorm(double x, double mu, double sigma) {
	return this->pnorm((x - mu) / sigma);
}




/////////////////
// LogisticLUT //
/////////////////

/*

void LogisticLUT::setup(const LogisticLUT::Config& cfg) {
	_config = cfg;

	_config.sdRange = std::min(_config.sdRange, 17.0);

	//if (_config.stepSize > 0.1) {
		// warn?
	//}

	size_t length = _config.sdRange / _config.stepSize;
	length += 10; // padding

	_densityLUT.resize(length);
	_cumulativeLUT.resize(length);

	for (size_t i = 0; i < length; i++) {

		double x = _config.stepSize * i;

		_densityLUT[i] = _config.dlogisFun(x);
		_cumulativeLUT[i] = _config.plogisFun(x);
	}
}


const LogisticLUT::Config& LogisticLUT::getConfig(void) {
	return _config;
}


double LogisticLUT::dlogis(double x) {
	// The logistic is symmetrical, so the sign of x doesn't matter
	x = std::abs(x);

	// If x is out of range, assume zero density
	if (x > _config.sdRange) {
		return 0.0;
	}

	x /= _config.stepSize; // This helps with precision vs doing just fmod(x, stepSize).

	double rem = fmod(x, 1);
	size_t i = (size_t)x;

	double lutI = _densityLUT[i];

	return lutI + rem * (_densityLUT[i + 1] - lutI);
}



double LogisticLUT::dlogis(double x, bool log) {

	double rval = this->dnorm(x);

	if (log) {
		rval = std::log(rval);
	}

	return rval;
}


double LogisticLUT::plogis(double x) {

	bool neg = x < 0;

	// Once the sign of x is noted, strip the sign to calculate index
	x = std::abs(x);

	// If x is out of range, return 0 if x < 0 and 1 if x > 0
	if (x > _config.sdRange) {
		return neg ? 0.0 : 1.0;
	}

	x /= _config.stepSize; // This helps with precision vs doing just fmod(x, stepSize).

	double rem = fmod(x, 1);
	size_t i = (size_t)x;

	double lutI = _cumulativeLUT[i];

	double lerpedValue = lutI + rem * (_cumulativeLUT[i + 1] - lutI);

	// If x < 0, then the integral should be flipped.
	if (neg) {
		lerpedValue = 1.0 - lerpedValue;
	}

	return lerpedValue;
}

*/


} // namespace CatCont
