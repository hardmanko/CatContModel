#include "UtilityFunctions.h"

#ifdef COMPILING_WITH_RCPP
using namespace R;
#endif

double inverseGammaLogLikelihood(double x, double alpha, double beta) {
	return alpha * std::log(beta) - lgammafn(alpha) - (alpha + 1) * std::log(x) - (beta / x);
}

double normal_varPostSample(const std::vector<double>& y, double mu, double a0, double b0) {

	double SSE = 0;
	for (unsigned int i = 0; i < y.size(); i++) {
		double dif = y[i] - mu;
		SSE += dif * dif;
	}

	double a = a0 + y.size() / 2.0;
	double b = b0 + SSE / 2;

	double r = rgamma(a, 1 / b); //this takes shape and scale, so 1/b, b being rate
	return 1 / r;
}

double normal_muPostSample(double mean, double N, double data_var, double mu0, double var0) {
	double a = N / data_var + 1 / var0;
	double b = mean * N / data_var + mu0 / var0;

	double sd1 = sqrt(1 / a);
	double mean1 = b / a;

	return rnorm(mean1, sd1);

	//return mean1 + sd1 * K4B_Bayes::_standardNormal(K4B_Bayes::_generator);
}

double normal_muPostSample(const std::vector<double>& y, double data_var, double mu0, double var0) {
#ifdef COMPILING_WITH_CX
	double mean = Util::mean(y);
#else
	double sum = 0;
	for (unsigned int i = 0; i < y.size(); i++) {
		sum += y[i];
	}
	double mean = sum / y.size();
#endif

	return normal_muPostSample(mean, y.size(), data_var, mu0, var0);
}

//goes from [0,1] to (-inf,inf)
double logitTransform(double unitInterval) {
	return std::log(unitInterval / (1 - unitInterval));
}

//goes from (-inf,inf) to [0,1]
double logitInverse(double fullScale) {
	double ex = exp(fullScale);
	return ex / (1 + ex);
}

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



