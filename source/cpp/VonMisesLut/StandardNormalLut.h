#pragma once

#include <vector>
#include <cmath>
#include <functional>


/*! 
Calculates dnorm and pnorm by using internal standard normal lookup tables.
The LUTs are one-sided to cut memory use in half.

sdRange and stepSize control the range and precision of the LUTs.
Maybe sdRange = 39 and stepSize = 0.0001.
*/
class StandardNormalLUT {
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

		std::function<double(double)> dnormFun;
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
