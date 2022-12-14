#pragma once

#include <functional>

#include "CCM_Util.h"
#include "CCM_Linear.h"
#include "CCM_Circular.h"

namespace CatCont {

#ifdef USING_PLAT_SPLINE

	double zeroDerivativeCubicSplineDensity(double z);
	double zeroDerivativeCubicSplineDensity2(double distFromPlatEdge, double splineHW);

	double dPlatSpline(double absDist, double platHW, double splineHW);

	// This is a very different function
	double dPlatSplineFull(double x, double mu, double platHW, double splineHW, bool linear = false, bool degrees = true);

	std::vector<double> dPlatSplineWeights(double study, CategoryParameters catPar, const ModelConfiguration& modCfg);

#endif

	//double calcPlatSplineLambda(const std::vector<double>& weights, LambdaVariant variant);

	//double dPlatSpline(double study, const CategoryParameters& catPar, const ModelConfiguration& modCfg);

/*
	struct PlatSplineWeightsCalculator {

		ModelConfiguration* modCfg;

		std::function<double(double, double)> distFun; //depends on DataType

		void setup(ModelConfiguration* modCfg);

		bool calcWeights(double study, const CategoryParameters& catPar);


		double dPlatSpline(double absDist, double platHW, double splineHW);
		// This is a very different function
		double dPlatSpline(double x, double mu, double platW, double splineHW, bool linear = false, bool degrees = true);

		double zeroDerivativeCubicSplineDensity(double z);

	};
*/


    struct WeightsCalculator {

		WeightsCalculator(ModelConfiguration* modCfg);

		ModelConfiguration* modCfg; // or ref/smart pointer?

		// Calculated weights
		size_t weightCount;
		vector<double> weights;

		// returns calculated weights in weightCount and weights members
		void calcWeights(double study, const CategoryParameters& catPar);

		// uses calculated weights
		double sumWeights(void) const;

#ifdef USING_PLAT_SPLINE
		double calcLambda(void) const; 
		//double calcLambda(double study, const CategoryParameters& catPar);

		// I dunno if this sould be member function. It calculates multiple weights.
		double calcLambdaIntegral(const vector<double>& studys, const CategoryParameters& catPar);

#endif

		vector<double> copyFilledWeights(void) const;
		void reset(void);

	private:

		void _calcWeights_platSpline(double study, const CategoryParameters& catPar);
		void _calcWeights_linear(double study, const CategoryParameters& catPar);
		void _calcWeights_circular(double study, const CategoryParameters& catPar);

		void _rescaleWeights(double densSum);
		
    };




	/*
	struct WeightsConfig {

		std::function<double(double, double, double)> densityFunction; // x, mu, catSel

	};

	
	struct WeightsManager {

		void setup(size_t maxCat, DataType dataType) {

			weights.resize(maxCat);

			if (dataType == DataType::Circular) {

				densityFunction = [](double study, double mu, double catSel) -> double {
					return vmLut.dVonMises(study, mu, catSel); // assume catSel is precision at this point
				};

			} else if (dataType == DataType::Linear) {

				//This distribution is not truncated: Category assignment is independent of the study/response space.
				densityFunction = Linear::normalPDF;
			}

		}

		//WeightsConfig config;
		std::function<double(double, double, double)> densityFunction; // x, mu, catSel

		// The density function for the PlatSpline method takes the whole vector of mus

		std::vector<double> weights;

		void calcWeights(double study, const CombinedParameters& par, const Linear::LinearConfiguration& lc) {

			unsigned int n = par.cat.mu.size();

			double densSum = 0;
			for (unsigned int i = 0; i < n; i++) {

				double d = this->densityFunction(study, par.cat.mu[i], par.cat.selectivity);
				densSum += d;
				this->weights[i] = d;
			}

			if (densSum < 1e-250 && !_usingLambda) {

				//If densSum is tiny, give equal weights but only if not using lambda
				for (unsigned int i = 0; i < n; i++) {
					this->weights[i] = 1.0 / n;
				}

			} else {
				for (unsigned int i = 0; i < n; i++) {
					this->weights[i] /= densSum;
				}
			}

			//return n;
		}


		vector<double> lambdaStudy; // The study values that are "integrated" over.

		void calcLambda(const CombinedParameters& par) {

		}

		// per participant?
		double getLambda(void) {
			if (!_usingLambda) {
				return 0;
			}

			return _lastLambda;
		}

	private:

		bool _usingLambda;
		double _lastLambda;


	};
	*/

}