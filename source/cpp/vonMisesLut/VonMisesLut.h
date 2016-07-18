#pragma once

#include <vector>
#include <cmath>


class VonMisesLut {
public:
	void setup(double rangeUpper_, double stepSize_, double(*besselFun)(double));

	double rangeUpper;
	double stepSize;

	double dVonMises(double x, double mu, double kappa);
	double dVonMises(double x, double mu, double kappa, bool log);

	double maxValue(void);

	double besselLerp(double x);

private:
	std::vector<double> _lut;

	unsigned int _maxI;
};

