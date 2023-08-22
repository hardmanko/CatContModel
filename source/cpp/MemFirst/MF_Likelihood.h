#pragma once

#include <functional>

#include "CCM_Main.h"
#include "CCM_Data.h"
#include "CCM_DistributionLUTs.h"

#include "CCM_Weights.h"
#include "CCM_Linear.h"
#include "CCM_Circular.h"

#include "MF_Parameters.h"
//#include "MF_ModelUtil.h" // for uniform


namespace CatCont {

	namespace MemFirst {

		namespace Circular {
			vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config);
			double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config);
			vector<double> zlLikelihood(const zlParameters& par, const ConditionData& data);
		}

		namespace Linear {
			vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config);
			double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, const ModelConfiguration& config);
		}

	}










	///////////////////////////////
	// I dunno about below
	/*

	struct LikelihoodCalculator {

		LikelihoodCalculator(CatCont::ModelConfiguration* modCfg);

		// functions for calculating various parts of the likelihood
		function<double(double, double, double)> _mainDensity; // either VM or normal. normal takes 5 args, so must be bound in.
		function<double(double, double, double)> _withinMean; // linear or circular

		vector<double> likelihoods;

		void calcLikelihoods(const ConditionData& data, const CombinedParameters& par);

		double sumLikelihoods(bool log);

		// ?
		void setErrorCallback(function<void(string)> cbFun);
		string getLastError(void);

	private:

		ModelConfiguration* _modCfg; // shared_ptr?

		// Constant w.r.t config
		double _calculatedUniformDensity;


		function<double(double, double, double)> _zlMemDens;
		void _zlLikelihood(const ConditionData& data, const zlParameters& par);

		// Different ideas
		void _bwLikelihood(const ConditionData& data, const CombinedParameters& par);

		void _bwLikelihood_linear(const ConditionData& data, const CombinedParameters& par);
		void _bwLikelihood_circular(const ConditionData& data, const CombinedParameters& par);

		double _bwLikelihood_within(const ConditionData& data, const CombinedParameters& par);
		double _bwLikelihood_between(const ConditionData& data, const CombinedParameters& par);
		function<double(double, double, double)> _bwMemDens;



		string _lastError;

	};

	*/

} // namespace CatCont