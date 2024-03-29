#include "UtilityFunctions.h"

#ifdef COMPILING_WITH_RCPP
using namespace R;
#endif

double inverseGammaLogLikelihood(double x, double alpha, double beta) {
	return alpha * std::log(beta) - lgammafn(alpha) - (alpha + 1) * std::log(x) - (beta / x);
}

// Sample from condition posterior distribution of variance (var_est) using an inverse-gamma conjugate prior.
// y are observations. y ~ Normal(mu_est, var_est)
// mu_est is estimated mu that pairs with var_est.
// a0 and b0 are prior shape and rate: var_est ~ InverseGamma(a0, b0).
double normal_varPostSample(const std::vector<double>& y, double mu_est, double a0, double b0) {

	double SSE = 0;
	for (size_t i = 0; i < y.size(); i++) {
		double dif = y[i] - mu_est;
		SSE += dif * dif;
	}

	double a = a0 + y.size() / 2.0;
	double b = b0 + SSE / 2;

	double r = rgamma(a, 1 / b); // rgamma takes shape and scale, so 1/b is scale, b being rate
	return 1 / r;
}

// Sample from condition posterior distribution of  location (mu_est) using a normal conjugate prior.
// y are observations. y ~ Normal(mu_est, var_est).
// var_est is the estimated variance that pairs with mu_est.
// mu0 and var0 are prior mean and variance: mu_est ~ Normal(mu0, var0).
double normal_muPostSample(const std::vector<double>& y, double var_est, double mu0, double var0) {
#ifdef COMPILING_WITH_CX
	double mean = Util::mean(y);
#else
	double sum = 0;
	for (size_t i = 0; i < y.size(); i++) {
		sum += y[i];
	}
	double mean = sum / y.size();
#endif

	return normal_muPostSample(mean, y.size(), var_est, mu0, var0);
}

// mean and N are the sample mean and count of observations.
double normal_muPostSample(double mean, double N, double var_est, double mu0, double var0) {
	double a = N / var_est + 1 / var0;
	double b = mean * N / var_est + mu0 / var0;

	double sd1 = sqrt(1 / a);
	double mean1 = b / a;

	return rnorm(mean1, sd1);
}


// goes from [0,1] to (-inf,inf)
// i.e., this is qlogis
double logitTransform(double unitInterval) {
	return std::log(unitInterval / (1 - unitInterval));
}

// goes from (-inf,inf) to [0,1]
// i.e. this is plogis
double logitInverse(double fullScale) {
	double ex = exp(fullScale);
	return ex / (1 + ex);
}

// x are counts??, p are probs.
//x and p must be the same length. sum(p) = 1
double multinomialProportionalLogLikelihood(const std::vector<unsigned>& x, const std::vector<double>& p) {
	double d = 0;
	for (unsigned int i = 0; i < p.size(); i++) {
		d += x[i] * std::log(p[i]);
	}
	return d;
}


double binomialProportionalLogLikelihood(double successes, double trials, double p) {
	return successes * std::log(p) + (trials - successes) * std::log(1 - p);
}



#ifdef COMPILING_WITH_CX

void outputDataVector(std::string filename, std::string header, std::vector<double> values) {
	std::ostringstream ss;
	ss << std::scientific << std::setprecision(20);
	ss << header << std::endl;
	for (unsigned int i = 0; i < values.size(); i++) {
		ss << values[i] << std::endl;
	}
	Util::writeToFile(filename, ss.str(), false);
}

#endif