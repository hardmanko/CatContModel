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


//OUT_weights must be an array of the correct size (par.cat.mu.size()). 
//This code is nasty for the sake of efficiency, avoiding lots of little allocations of weights vectors.
//This function is called  in the likelihood function for every observation, so it gets a lot of use.
void categoryWeights(double study, const CombinedParameters& par, const LinearConfiguration& lc, double* OUT_weights) {

	unsigned int n = par.cat.mu.size();

	double densSum = 0;
	for (unsigned int i = 0; i < n; i++) {
		//This distribution is not truncated: Category assignment is independent of the study/response space.
		double d = normalPDF(study, par.cat.mu[i], par.cat.selectivity);
		densSum += d;
		OUT_weights[i] = d;
	}

	if (densSum < 1e-250) {

		//If densSum is tiny, give equal weights
		for (unsigned int i = 0; i < n; i++) {
			OUT_weights[i] = 1.0 / n;
		}
	} else {
		for (unsigned int i = 0; i < n; i++) {
			OUT_weights[i] /= densSum;
		}
	}
}

double combineSds(double cont, double cat, double pContWithin) {
	double combinedVar = pow(pContWithin * cont, 2) + pow((1 - pContWithin) * cat, 2);
	return sqrt(combinedVar);
}

vector<double> zlLikelihood(const zlParameters& par, const ConditionData& data, const ModelConfiguration& config) {

	const LinearConfiguration& lc = config.linearConfiguration;

	const double unifDens = 1 / (lc.response.upper - lc.response.lower);
	const double guessDens = (1 - par.pMem) * unifDens;

	vector<double> likelihoods(data.study.size());

	for (unsigned int i = 0; i < data.study.size(); i++) {

		double memDens = par.pMem * dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

		likelihoods[i] = memDens + guessDens;
	}

	return likelihoods;
}

vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config) {

	const LinearConfiguration& lc = config.linearConfiguration;

	if (config.modelVariant == ModelVariant::ZL) {
		zlParameters zlPar;
		zlPar.pMem = par.pMem;
		zlPar.contSD = par.contSD;

		return zlLikelihood(zlPar, data, config);
	}

	bool calculateWithinComponent = (config.modelVariant == ModelVariant::BetweenAndWithin) || (config.modelVariant == ModelVariant::WithinItem);
	bool calculateBetweenComponent = (config.modelVariant == ModelVariant::BetweenAndWithin) || (config.modelVariant == ModelVariant::BetweenItem);

	double pBetween = par.pBetween;
	if (config.modelVariant == ModelVariant::BetweenItem) {
		pBetween = 1;
	} else if (config.modelVariant == ModelVariant::WithinItem) {
		pBetween = 0;
	}

	//The width of the uniform depends on the response range
	const double unifDens = 1 / (lc.response.upper - lc.response.lower);


	vector<double> weights(par.cat.mu.size());

	//Precalculate this because it doesn't depend on the data
	double combinedWithinSd = combineSds(par.contSD, par.cat.SD, par.pContWithin);

	vector<double> likelihoods(data.study.size());

	for (unsigned int i = 0; i < data.study.size(); i++) {

		//Category weights apply to both the between and within components of the model
		categoryWeights(data.study[i], par, lc, weights.data());

		//Within component
		double withinDensity = 0;

		if (calculateWithinComponent) {

			if (par.cat.mu.size() == 0) {

				withinDensity = dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

			} else {

				for (unsigned int j = 0; j < weights.size(); j++) {

					double predictedCenter = par.pContWithin * data.study[i] + (1 - par.pContWithin) * par.cat.mu[j];

					double thisDens = dtnorm_noBoundCheck(data.response[i], predictedCenter, combinedWithinSd, lc.response.lower, lc.response.upper);

					withinDensity += weights[j] * thisDens;
				}
			}
		}

				
		double betweenCatDens = 0;
		double catGuessDens = 0;

		//For both between and guessing. This is always calculated for guessing even if between isn't calculated.
		for (unsigned int j = 0; j < weights.size(); j++) {
			double dens = dtnorm_noBoundCheck(data.response[i], par.cat.mu[j], par.cat.SD, lc.response.lower, lc.response.upper);
			betweenCatDens += weights[j] * dens;
			catGuessDens += dens;
		}

		//Between component
		double betweenDensity = 0;
		if (calculateBetweenComponent) {
			double contDensity = dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

			betweenDensity = (par.pContBetween * contDensity) + ((1 - par.pContBetween) * betweenCatDens);
		}


		//Guessing component
		if (par.cat.mu.size() > 0) {
			catGuessDens /= par.cat.mu.size();
		}

		double guessingDensity = par.pCatGuess * catGuessDens + (1 - par.pCatGuess) * unifDens;


		//Combine components
		double memoryDensity = pBetween * betweenDensity + (1 - pBetween) * withinDensity;

		double likelihood = par.pMem * memoryDensity + (1 - par.pMem) * guessingDensity;

		likelihoods[i] = likelihood;

	}

	return likelihoods;

}

double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config) {
	vector<double> likelihoods = betweenAndWithinLikelihood(par, data, config);

	double llSum = 0;
	for (unsigned int i = 0; i < likelihoods.size(); i++) {
		llSum += log(likelihoods[i]);
	}

	return llSum;
}

} // namespace Linear
} // namespace CatCont