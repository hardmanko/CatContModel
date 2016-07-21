#include "CCM_Linear.h"

#ifdef COMPILING_WITH_RCPP
#define R_TRUE 1
#define R_FALSE 0
#endif

namespace CatCont {
	namespace Linear {

		/*
		double dtnorm_denominator(double mu, double sd, double lower, double upper) {
			double a;
			double b;

			if (lower == -std::numeric_limits<double>::infinity()) {
				a = 0;
			} else {
				a = R::pnorm(lower, mu, sd, R_TRUE, R_FALSE); //q, mu, sd, lower tail=T, log.p=F
			}

			if (upper == std::numeric_limits<double>::infinity()) {
				b = 1;
			} else {
				b = R::pnorm(upper, mu, sd, R_TRUE, R_FALSE);
			}

			return b - a;
		}
		*/

		//This does NOT test for infinity (which pnorm probably does anyway)
		double dtnorm_denominator(double mu, double sd, double lower, double upper) {

			double a = R::pnorm(lower, mu, sd, R_TRUE, R_FALSE); //q, mu, sd, lower tail=T, log.p=F
			double b = R::pnorm(upper, mu, sd, R_TRUE, R_FALSE);

			return b - a;
		}

		double dnorm(double x, double mu, double sd) {
			return R::dnorm(x, mu, sd, R_FALSE);
		}

		//Does not enforce lower < upper, but why would that ever happen?
		double dtnorm(double x, double mu, double sd, double lower, double upper, bool log_) {

			//This may not be needed as the data bounds always contain the data...
			if (x < lower || x > upper) {
				return 0;
			}

			double denominator = dtnorm_denominator(mu, sd, lower, upper);

			double d = dnorm(x, mu, sd) / denominator;
			if (log_) {
				d = std::log(d);
			}
			return d;
		}

		double dtnorm(double x, double mu, double sd, double lower, double upper) {

			//This may not be needed as the data bounds always contain the data...
			if (x < lower || x > upper) {
				return 0;
			}

			double denominator = dtnorm_denominator(mu, sd, lower, upper);

			return R::dnorm(x, mu, sd, R_FALSE) / denominator;
		}

		//This does NOT check that x is in the interval. It assumes that the data range is wider than the
		//observed data, which is enforced in R.
		double dtnorm_noBoundCheck(double x, double mu, double sd, double lower, double upper) {

			double denominator = dtnorm_denominator(mu, sd, lower, upper);

			return R::dnorm(x, mu, sd, R_FALSE) / denominator;
		}


		vector<double> dtnorm(const vector<double>& x, double mu, double sd, double lower, double upper, bool log_) {

			double denominator = dtnorm_denominator(mu, sd, lower, upper);

			vector<double> dens(x.size());

			for (unsigned int i = 0; i < x.size(); i++) {

				if (x[i] < lower || x[i] > upper) {
					dens[i] = 0;
				}
				
				double d = dnorm(x[i], mu, sd) / denominator;
				if (log_) {
					d = std::log(d);
				}
				dens[i] = d;

			}

			return dens;
		}

		double dtnorm_unifEdge(double x, double mu, double sd, double lower, double upper, double edgeWidth) {

			double probInInterval = dtnorm_denominator(mu, sd, lower, upper);
			double excessProb = 1 - probInInterval;

			double normalDens = R::dnorm(x, mu, sd, R_FALSE);
			double excessUnifDens = 0;
			double linearCenter = (upper + lower) / 2;

			if (mu > linearCenter) {
				if (x >= upper - edgeWidth) {
					excessUnifDens = excessProb / edgeWidth;
				}
			} else {
				if (x <= lower + edgeWidth) {
					excessUnifDens = excessProb / edgeWidth;
				}
			}

			double dens = normalDens + excessUnifDens;

			return dens;
		}

		//This isn't really important
		double dtnorm_tnormEdge(double x, double mu, double sd, double lower, double upper, double edgeSD) {

			if (x < lower || x > upper) {
				return 0;
			}

			double probInInterval = dtnorm_denominator(mu, sd, lower, upper);

			double probOutsideInterval = 1 - probInInterval;

			double normalDens = R::dnorm(x, mu, sd, R_FALSE);
			double linearCenter = (upper + lower) / 2;
			double edgeNormalDens;

			if (mu > linearCenter) {
				edgeNormalDens = R::dnorm(x, upper, edgeSD, R_FALSE);
			} else {
				edgeNormalDens = R::dnorm(x, lower, edgeSD, R_FALSE);
			}

			//You know that half of the edge normal's mass is outside of the interval, so multiply by 2 to scale it up to have a mass of 1.
			//Then multiply by the amount of mass outside of the interval to scale it back down.
			edgeNormalDens *= 2 * probOutsideInterval;

			double dens = normalDens + edgeNormalDens;

			return dens;
		}

		//If mu and sd put the distribution away from the bounds so far that the proportion within the bounds is tiny,
		//return 0 density.
		//NOTE: Don't use this, it's conceptually difficult. Should pWithinInterval < 0.01?
		double dtnorm_clampToZero(double x, double mu, double sd, double lower, double upper) {

			double pWithinInterval = dtnorm_denominator(mu, sd, lower, upper);
			if (pWithinInterval < 1e-100) {
				return 0;
			}

			double dens = R::dnorm(x, mu, sd, R_FALSE);

			return dens / pWithinInterval;
		}

		//This limits how much the density can be scaled by only allowing the prob within interval to be so small.
		double dtnorm_limitProbWithin(double x, double mu, double sd, double lower, double upper) {
			double pWithinInterval = dtnorm_denominator(mu, sd, lower, upper);
			pWithinInterval = std::max(pWithinInterval, 0.5); //No less that a 0.5 scale value.
			//0.5 is an important point. if x == upper, then if p < 0.5, the density is higher when mu > x than when mu == x.
			//If p == 0 (i.e. a standard truncated normal), then dens increases as mu goes beyond x. This is stupid.

			double dens = R::dnorm(x, mu, sd, R_FALSE);

			return dens / pWithinInterval;
		}


		//OUT_weights must be an array of the correct size (par.cat.mu.size()). 
		//This code is nasty for the sake of efficiency, avoiding lots of little allocations of weights vectors.
		//This function is called  in the likelihood function for every observation, so it gets a lot of use.
		void categoryWeights(double study, const CombinedParameters& par, const LinearConfiguration& lc, double* OUT_weights) {

			unsigned int n = par.cat.mu.size();

			double densSum = 0;
			for (unsigned int i = 0; i < n; i++) {
				//There is an argument that this distribution is not truncated: Category assignment is perceptual, independent of the study/response space.
				double d = R::dnorm(study, par.cat.mu[i], par.cat.selectivity, R_FALSE);
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

		vector<double> zlLikelihood(const zlParameters& par, const ConditionData& data, const LinearConfiguration& lc) {

			const double unifDens = 1 / (lc.response.upper - lc.response.lower);
			const double guessDens = (1 - par.pMem) * unifDens;

			vector<double> likelihoods(data.study.size());

			for (unsigned int i = 0; i < data.study.size(); i++) {

				double memDens = par.pMem * dtnorm_noBoundCheck(data.response[i], data.study[i], par.contSD, lc.response.lower, lc.response.upper);

				likelihoods[i] = memDens + guessDens;
			}

			return likelihoods;
		}

		vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const LinearConfiguration& lc) {

			if (lc.modelVariant == ModelVariant::ZL) {
				zlParameters zlPar;
				zlPar.pMem = par.pMem;
				zlPar.contSD = par.contSD;

				return zlLikelihood(zlPar, data, lc);
			}

			bool calculateWithinComponent = (lc.modelVariant == ModelVariant::BetweenAndWithin) || (lc.modelVariant == ModelVariant::WithinItem);
			bool calculateBetweenComponent = (lc.modelVariant == ModelVariant::BetweenAndWithin) || (lc.modelVariant == ModelVariant::BetweenItem);

			double pBetween = par.pBetween;
			if (lc.modelVariant == ModelVariant::BetweenItem) {
				pBetween = 1;
			} else if (lc.modelVariant == ModelVariant::WithinItem) {
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

		double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const LinearConfiguration& config) {
			vector<double> likelihoods = betweenAndWithinLikelihood(par, data, config);

			double llSum = 0;
			for (unsigned int i = 0; i < likelihoods.size(); i++) {
				llSum += log(likelihoods[i]);
			}

			return llSum;
		}

	} // namespace Linear
} // namespace CatCont