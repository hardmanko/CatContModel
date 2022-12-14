#pragma once

#include "Compilation.h"

#if COMPILING_WITH_CX
#include "CX.h"
#include "ofxRMath.h"
#endif

//#define USING_CAT_SPLINE


#include <vector>
#include <map>
#include <string>


using namespace std;

#ifndef PI
	#define PI       3.14159265358979323846
#endif

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

	enum class WeightsDistribution : int {
		Default,
		PlatSpline
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


	double normalLL(double x, double mu, double var);
	double cauchyLL(double x, double loc, double scale);
	double normalDeviate(double x, double sd);
	double uniformDeviate(double low, double high);

	double clamp(double x, double minimum, double maximum);

	// Unused
	string extractIndex(string s);
	string extractBaseParameterName(string s);

	string joinIndexedParam(string name, vector<string> index);
	vector<string> splitIndexedParam(string joined);
	vector<string> splitString(string str, string delim);

	void logMessage(string module, string message, bool endLine = true);

	struct ConditionData {
		string condition; // redundant with Data::conditionNames

		vector<double> study;
		vector<double> response;
	};

	struct ParticipantData {
		string pnum;
		vector<ConditionData> condData;
	};

	vector<ParticipantData> copyParticipantData(vector<string> pnumsCol, vector<string> condsCol, 
		vector<double> studyCol, vector<double> responseCol, DataType dataType, bool verbose);

#if COMPILING_WITH_CX
	vector<ParticipantData> getParticipantData(string filename, DataType dataType, bool verbose);
#endif
#if COMPILING_WITH_RCPP
	vector<ParticipantData> getParticipantData(Rcpp::DataFrame df, DataType dataType, bool verbose);
#endif


	struct Data {
		Data(void) {
			studyRange.lower = numeric_limits<double>::max();
			studyRange.upper = numeric_limits<double>::min();

			responseRange.lower = numeric_limits<double>::max();
			responseRange.upper = numeric_limits<double>::min();
		}

		struct {
			double lower;
			double upper;
		} studyRange;

		struct {
			double lower;
			double upper;
		} responseRange;

		vector<ParticipantData> participants;
		vector<string> conditionNames;
	};




	struct zlParameters {
		double pMem;
		double contSD;
	};

#ifdef USING_CAT_SPLINE
	struct PlatSplineParameters {
		vector<double> platHW; // half width

		vector<double> splineHW; // half width (or standard deviation)

		vector<double> height;

		size_t nCat(void) const {
			if (!isValid()) {
				return 0;
			}

			return platHW.size();
		}

		bool isValid(void) const {
			return platHW.size() == splineHW.size() && platHW.size() == height.size();
		}
	};
#endif
	

	struct CategoryParameters {

		size_t nCat(void) const {
			return mu.size();
		}

		// Parameters related to categorization (the weights function)

		vector<double> mu; // catMu: Location of the category
		double selectivity; // catSel: How categories transition into one another
		
#ifdef USING_CAT_SPLINE
		vector<double> platHW; // PlatSpline plateau half width
		//double beta; // PlatSpline lambda multiplier. Not really a category parameter, but closely related to categorization.
		//PlatSplineParameters psPar;
#endif

		// Parameters related to memory.

		double SD; // catSD: Memory precision for category parameters
	};



	struct Parameters {

		double pMem;
		double pBetween;
		double pContBetween;
		double pContWithin;

		double pCatGuess;

		double contSD;

		CategoryParameters cat;
	};

	struct ConditionCategoryParameters {
		double selectivity;
		double SD;

		//double beta; // lambda
	};

	struct ConditionParameters {
		double pMem;
		double pBetween;
		double pContBetween;
		double pContWithin;

		double pCatGuess;

		double contSD;

		ConditionCategoryParameters cat;
	};

	//typedef Parameters ParticipantParameters;
	//typedef Parameters CombinedParameters;

	struct ParticipantParameters : Parameters {};
	struct CombinedParameters : Parameters {};

	namespace Linear {
		struct LinearConfiguration {

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

	struct PrivateConfig {
		bool useVonMisesLookupTable = true; // This doesn't appear to be used
	};

	struct RunConfig {
		unsigned int iterations;
		unsigned int iterationsPerStatusUpdate;

		bool verbose = true;
		bool profileParameterTypes = false; // profile parameter timing?
	};

	struct ModelConfiguration {

		//unsigned int iterations;
		//unsigned int iterationsPerStatusUpdate;

		DataType dataType;
		ModelVariant modelVariant;

		// Category config
		unsigned int maxCategories;
		unsigned int catMuPriorApproximationPrecision;

		// Experimental category sharing stuff
		bool catMuShared = false;
		bool catActiveShared = false; // catActive can only be shared if catMu is shared. Possibilities: neither shared, catMu shared, both shared.
		bool catActiveHierPrior = false; // can only be true if catMuShared && !catActiveShared. No parameters are needed for the prior.
		bool catActiveDistancePrior = true; // joint prior on catMu and catActive.

		string cornerstoneConditionName;
		unsigned int cornerstoneConditionIndex;
		// sums-to-zero option? (might be only practical for ZL)

		SDRanges ranges;

		bool calculateParticipantLikelihoods;

		map<string, vector<string> > conditionEffects;
		vector<string> paramWithConditionEffects;

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
			map<string, string> equalityConstraints;
		} overrides;

		// Should not be used by normal users of the package
		PrivateConfig privateConfig;

#ifdef USING_CAT_SPLINE
		WeightsDistribution weightsDistribution;
		LambdaVariant lambdaVariant;
#endif
	};



	template <typename T>
	string _convertToString(T val) {
		ostringstream oss;
		oss << val;
		return oss.str();
	}

	template <typename T>
	T _convertFromString(string s) {
		T data;
		stringstream ss;
		ss << s;
		ss >> data;
		return data;
	}

	template<typename T>
	vector<T> uniqueElements(vector<T> v) {
		std::sort(v.begin(), v.end());
		auto lastIndex = std::unique(v.begin(), v.end());
		v.erase(lastIndex, v.end());
		return v;
	}


#ifdef COMPILING_WITH_CX
	void useModel_Bayesian(void); // Freestanding interface
#endif

}