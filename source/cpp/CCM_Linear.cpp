#include "CCM_Linear.h"

#ifdef COMPILING_WITH_RCPP
#define R_TRUE 1
#define R_FALSE 0
#endif

namespace CatCont {
namespace Linear {

double normalPDF(double x, double mu, double sd) {
#ifdef COMPILING_WITH_RCPP
	return R::dnorm(x, mu, sd, R_FALSE); // log = false
#else
	return dnorm(x, mu, sd, Rboolean::FALSE); // log = false
#endif
}

double normalCDF(double x, double mu, double sd) {
#ifdef COMPILING_WITH_RCPP
	return R::pnorm(x, mu, sd, R_TRUE, R_FALSE); // lowerTail = TRUE, log = FALSE
#else
	return pnorm(x, mu, sd, Rboolean::TRUE, Rboolean::FALSE); // lowerTail = TRUE, log = FALSE
#endif
}


//This does NOT test for infinity (which pnorm does anyway)
double dtnorm_denominator(double mu, double sd, double lower, double upper) {

	double a = normalCDF(lower, mu, sd);
	double b = normalCDF(upper, mu, sd);

	return b - a;
}


//Does not enforce lower < upper, but why would that ever happen?
double dtnorm(double x, double mu, double sd, double lower, double upper, bool log_) {

	// This bounds check is not needed as the data bounds always contain the data.
	if (x < lower || x > upper) {
		return 0;
	}

	double denominator = dtnorm_denominator(mu, sd, lower, upper);

	double d = normalPDF(x, mu, sd) / denominator;
	if (log_) {
		d = std::log(d);
	}
	return d;
}

vector<double> dtnorm(const vector<double>& x, double mu, double sd, double lower, double upper, bool log_) {

	double denominator = dtnorm_denominator(mu, sd, lower, upper);

	vector<double> dens(x.size());

	for (unsigned int i = 0; i < x.size(); i++) {

		if (x[i] < lower || x[i] > upper) {
			dens[i] = 0;
		}

		double d = normalPDF(x[i], mu, sd) / denominator;
		if (log_) {
			d = std::log(d);
		}
		dens[i] = d;

	}

	return dens;
}


//This does NOT check that x is in the interval. It assumes that the data range is wider than the
//observed data, which is enforced in the R interface.
double dtnorm_noBoundCheck(double x, double mu, double sd, double lower, double upper) {

	double denominator = dtnorm_denominator(mu, sd, lower, upper);

	return normalPDF(x, mu, sd) / denominator;
}


//This limits how much the density can be scaled by only allowing the prob within interval to be so small.
double dtnorm_limitProbWithin(double x, double mu, double sd, double lower, double upper) {
	double pWithinInterval = dtnorm_denominator(mu, sd, lower, upper);

	// No less than a 0.5 scale value.
	pWithinInterval = std::max(pWithinInterval, 0.5); 
	// 0.5 is an important point. if x == upper, then if p < 0.5, the density is higher when mu > x than when mu == x.
	// For a standard truncated normal, the dens increases as mu goes beyond x, which is wrong for the purposes of CatContModel.

	double dens = normalPDF(x, mu, sd);

	return dens / pWithinInterval;
}


double combineSDs(double pContWithin, double contSD, double catSD) {
	double combinedVar = pow(pContWithin * contSD, 2) + pow((1 - pContWithin) * catSD, 2);
	return sqrt(combinedVar);
}


} // namespace Linear
} // namespace CatCont