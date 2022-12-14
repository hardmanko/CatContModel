#pragma once


#include <limits>
#include <vector>
#include <cmath>
#include <functional>


namespace CatCont {

class VonMisesLUT {
public:

	struct Config {

		// A bessel i function of the 0th order with the exponent scaled:
		// e.g. R::bessel_i(x, 0, 2); // where 2 == true
		std::function<double(double)> besselFun;

		bool useLUT = true;
		bool skipRecomputationIfAble = true;

		double maxKappa;
		double stepSize;
		
	};

	void setup(std::function<double(double)> besselFun);
	void setup(const Config& cfg);
	const Config& getConfig(void) const;

	bool ready(double maxKappa = -1.0);

	// Radians and precision
	double dVonMises(double x, double mu, double kappa);
	double dVonMises(double x, double mu, double kappa, bool log);

	double besselLerp(double kappa);

	double getMaxKappa(void) const;
	

private:
	Config _config;

	std::vector<double> _besselLUT;
	std::function<double(double)> _internalBessel;
};



/*!
Calculates dnorm and pnorm by using internal standard normal lookup tables.
The LUTs are one-sided to cut memory use in half.

sdRange and stepSize control the range and precision of the LUTs.
Maybe sdRange = 39 and stepSize = 0.0001.
*/
class NormalLUT {
public:

	struct Config {

		// LUT range in standard deviation units.
		// Maximum value of 39, where dnorm(39) is approximately 0 in 64-bit floating point precision.
		// dnorm(x) == 0, x < -sdRange || x > sdRange
		// pnorm(x) == 0, x < -sdRange 
		// pnorm(x) == 1, x > sdRange.
		double sdRange = 39;

		// Steps are taken until sdRange is reached.
		// Recommend values between 0.001 and 0.00001 (10^-3 to 10^-5).
		double stepSize = 1e-4;

		// Standard normal density
		std::function<double(double)> dnormFun;

		// Standard normal cumulative distribution
		std::function<double(double)> pnormFun;
	};

	void setup(const Config& config);
	const Config& getConfig(void);

	double dnorm(double x);
	double dnorm(double x, double mu, double sigma);
	double dnorm(double x, double mu, double sigma, bool log);

	double pnorm(double x);
	double pnorm(double x, double mu, double sigma);

	// TODO: Implement
	//double estimateDnormBias(std::vector<double> testXs);
	//double estimatePnormBias(std::vector<double> testXs);

private:

	Config _config;

	std::vector<double> _dnormLUT;
	std::vector<double> _pnormLUT;
};

/*
class LogisticLUT {
public:

	struct Config {

		// LUT range in standard deviation units.
		// Maximum value of 17, where dnorm(17) is approximately 0 in 64-bit floating point precision.
		// dnorm(x) == 0, x < -sdRange || x > sdRange
		// pnorm(x) == 0, x < -sdRange 
		// pnorm(x) == 1, x > sdRange.
		double sdRange = 17;

		// Steps are taken until sdRange is reached.
		// Recommend values between 0.001 and 0.00001 (10^-3 to 10^-5).
		double stepSize = 1e-4;

		// Logistic density function
		std::function<double(double)> dlogisFun;

		// Logistic cumulative distribution
		std::function<double(double)> plogisFun;
	};

	void setup(const Config& config);
	const Config& getConfig(void);

	double dlogis(double x);
	double dlogis(double x, bool log);

	double plogis(double x);

	// TODO: Implement
	//double estimateDensityBias(std::vector<double> testXs);
	//double estimateCumulativeBias(std::vector<double> testXs);

private:

	Config _config;

	std::vector<double> _densityLUT;
	std::vector<double> _cumulativeLUT;
};
*/

extern VonMisesLUT vmLut;
//extern NormalLUT normalLUT;
//extern LogisticLUT logitLUT;

} // namespace CatCont