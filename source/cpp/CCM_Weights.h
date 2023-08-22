#pragma once

#include <functional>

#include "CCM_Main.h"
#include "CCM_Linear.h"
#include "CCM_Circular.h"

#include "MF_Parameters.h"

//#define USING_CAT_SPLINE

std::vector<double> CCM_CPP_CalculateWeights(std::string dataType, std::string weightType, double study, std::vector<double> activeCatMu, double catSelectivity);

namespace CatCont {

	enum class WeightsDistribution : int {
		Default, // HVR17
		Nearest, // no parameters
		NearestInRange // temp using catSel as range, but should use a plateau parameter
#ifdef USING_CAT_SPLINE
		, 
		PlatSpline // multiple parameters. unimplemented
#endif
	};
	WeightsDistribution weightsDistributionFromString(string weightsDistStr);

#ifdef USING_CAT_SPLINE
	enum class LambdaVariant : int {
		None, // lambda is always 0
		CatWeightSum, // 0 to 1, increases with cat weight
		InverseCatWeightSum, // 0 to 1, decreases with cat weight
		Minus1to1 // -1 to 1, increases with cat weight
	};
	LambdaVariant lambdaVariantFromString(string lambdaVariantStr);
#endif



	struct WeightsConfig {
		DataType dataType;
		size_t maxCategories;

		WeightsDistribution weightDist = WeightsDistribution::Default;

		double zeroSumCutoff = 1e-250; // If the weight sum is below this, weights are set to equal values.
		double zapsmallCutoff = 0; // If a (rescaled) weight is below this, it is set to 0. Maybe 1e-6?
	};

    struct WeightsCalculator {

		WeightsCalculator(const ModelConfiguration& modCfg);
		WeightsCalculator(const WeightsConfig& wConfig);

		WeightsConfig config;

		//shared_ptr<ModelConfiguration> modCfg; // or ref/smart pointer?
		//DataType dataType;
		//double zeroSumCutoff = 1e-250; // If the weight sum is below this, weights are set to equal values.
		//double zapsmallCutoff = 0; // If a (rescaled) weight is below this, it is set to 0. Maybe 1e-6?

		// Calculated weights
		
		vector<double> weights; // The weights vector is sized on construction to be large enough for any number of weights.
		size_t weightCount = 0; // Weight count tracks how many weights are stored.

		// returns calculated weights in weightCount and weights members
		void calcWeights(double study, const vector<double>& catMu, double catSelectivity);

		vector<double> copyFilledWeights(void) const;
		void clearWeights(void);

		// uses calculated weights
		double sumWeights(void) const;

#ifdef USING_PLAT_SPLINE
		double calcLambda(void) const; 
		//double calcLambda(double study, const CategoryParameters& catPar);

		// I dunno if this sould be member function. It calculates multiple weights.
		double calcLambdaIntegral(const vector<double>& studys, const CategoryParameters& catPar);

		void _calcWeights_platSpline(double study, const CategoryParameters& catPar);

#endif



	private:
		
		void _calcWeights_default_HVR17(double study, const vector<double>& catMu, double catSelectivity);
		//void _calcWeights_linear(double study, const vector<double>& catMu, double catSelectivity);
		//void _calcWeights_circular(double study, const vector<double>& catMu, double catSelectivity);
		void _rescaleWeights(double densSum);
		void _zapSmallWeights(void);

		void _calcWeights_nearest(double study, const vector<double>& catMu);
		void _calcWeights_nearestInRange(double study, const vector<double>& catMu, double catRange);
		std::vector<double> _studyMuDistance(double study, const vector<double>& catMu) const;

		
    };


	// These functions are the old way of doing things.
	namespace Circular {
		void categoryWeights(double study, const MemFirst::CombinedParameters& par, double* OUT_weights);
	}
	namespace Linear {
		void categoryWeights(double study, const MemFirst::CombinedParameters& par, const LinearConfiguration& lc, double* OUT_weights);
	}

#ifdef USING_PLAT_SPLINE

	double zeroDerivativeCubicSplineDensity(double z);
	double zeroDerivativeCubicSplineDensity2(double distFromPlatEdge, double splineHW);

	double dPlatSpline(double absDist, double platHW, double splineHW);

	// This is a very different function
	double dPlatSplineFull(double x, double mu, double platHW, double splineHW, bool linear = false, bool degrees = true);

	std::vector<double> dPlatSplineWeights(double study, CategoryParameters catPar, const ModelConfiguration& modCfg);

	double calcPlatSplineLambda(const std::vector<double>& weights, LambdaVariant variant);

	double dPlatSpline(double study, const CategoryParameters& catPar, const ModelConfiguration& modCfg);

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

#endif

}