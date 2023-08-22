// Things related to model-general configuration

#pragma once

#include "CCM_Main.h"

#include <set>

using namespace std;

namespace CatCont {


	enum class ModelVariant : int {
		BetweenAndWithin, //between and within
		BetweenItem, //between item variant
		WithinItem, //within item variant, 
		ZL //zhang and luck model
	};
	ModelVariant modelVariantFromString(string modelVariantStr);

	enum class DataType : int {
		Circular,
		Linear
	};
	DataType dataTypeFromString(string dataTypeStr);




	///////////////////////////////////////////////////////////////////////////////////////////////
	// Stuff related to linear config

	namespace Linear {
		struct LinearConfiguration {

			bool setFromList(Rcpp::List cfgList);

			struct {
				double lower;
				double upper;
			} response;
			// Study ranges are not needed, only response

			struct {
				double lower;
				double upper;
			} catMu;

		};
	}

	struct SDRanges {
		double minSd;
		double maxSd;

		double minPrecision;
		double maxPrecision;
	};



	///////////////////////////////////////////////////////////////////////////////////////////////
	// Model and Run configurations

	struct RunConfig {

		bool setFromList(Rcpp::List cfgList);

		unsigned int iterations;
		unsigned int iterationsPerStatusUpdate = 20;

		bool verbose = true;
		bool profileParameterTypes = false; // profile parameter timing?
	};

	struct PrivateConfig {
		bool useVonMisesLookupTable = true; // This doesn't appear to be used
	};

	struct ModelConfiguration {

		bool setFromList(Rcpp::List cfgList);

		DataType dataType;
		ModelVariant modelVariant;

		// Parameters that are shared by all participants.
		std::set<string> sharedParameters; // set is easily searchable.

		// Category config
		unsigned int maxCategories;
		unsigned int catMuPriorApproximationPrecision;

		SDRanges sdRanges;

		// Conditions
		string cornerstoneConditionName;

		map<string, vector<string> > conditionEffects;
		vector<string> paramWithConditionEffects;

		// TODO: These functions seem like model-specific rather than model-general and maybe should be moved.
		vector<string> getParamWithAndWithoutConditionEffects(void) const;
		vector<string> getParamWithConditionEffects(void) const;
		vector<string> getParamWithoutConditionEffects(void) const;
		vector<string> getParamWithHierachicalPriors(void) const;

		Linear::LinearConfiguration linearConfiguration;

		struct {
			map<string, double> mhTunings;
			map<string, double> priors;
			map<string, double> startingValues;
			map<string, double> constantValues;
			map<string, string> equalityConstraints; // TODO: Keep or drop?
		} overrides;

		bool calculateParticipantLikelihoods = false;

		////////////////////////////////////
		// Things below here are experimental and are only used for testing. 
		// They should not be used by normal users of the package.
		PrivateConfig privateConfig;

		// Experimental category sharing stuff
		//bool catMuShared = false;
		//bool catActiveShared = false; // catActive can only be shared if catMu is shared. Possibilities: neither shared, catMu shared, both shared.
		bool catActiveHierPrior = false; // can only be true if catMuShared && !catActiveShared. The only parameter needed for the prior is something to control the range of probabilities.
		bool catActiveDistancePrior = true; // joint prior on catMu and catActive. This is the Hardman, Vergauwe, & Ricker 2017 default.

		// For within-item model
		bool withinItem_contSD_is_withinSD = false; // If true, contSD is freely estimated separately from catSD, not combined with catSD. (i.e. treat contSD as withinSD.)

#ifdef USING_CAT_SPLINE
		WeightsDistribution weightsDistribution;
		LambdaVariant lambdaVariant;
#endif
	};

} // namespace CatCont
