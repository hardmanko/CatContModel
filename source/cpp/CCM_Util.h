#pragma once

#include "Compilation.h"

#if COMPILING_WITH_CX
#include "CX.h"
#include "ofxRMath.h"
#endif

#include "VonMisesLut.h"

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

	enum class DataType : int {
		Circular,
		Linear
	};

	DataType dataTypeFromString(string dataTypeString);
	ModelVariant modelVariantFromString(string modelVariantStr);
	

	extern VonMisesLut vmLut;


	double normalLL(double x, double mu, double var);
	double cauchyLL(double x, double loc, double scale);
	double normalDeviate(double x, double sd);
	double uniformDeviate(double low, double high);

	double cubicSplineDensity(double x, double scale);

	double clamp(double x, double minimum, double maximum);

	string extractIndex(string s);
	string extractBaseParameterName(string s);

	void logMessage(string module, string message, bool endLine = true);

	struct ConditionData {
		string condition;

		vector<double> study;
		vector<double> response;
	};

	struct ParticipantData {
		string pnum;
		vector<ConditionData> condData;
	};

	vector<ParticipantData> copyParticipantData(vector<string> pnumsCol, vector<string> condsCol, vector<double> studyCol, vector<double> responseCol, DataType dataType, bool verbose);
#if COMPILING_WITH_CX
	vector<ParticipantData> getParticipantData(string filename, DataType dataType, bool verbose);
#endif
#if COMPILING_WITH_RCPP
	vector<ParticipantData> getParticipantData(Rcpp::DataFrame df, DataType dataType, bool verbose);
#endif


	struct Data {
		Data(void) {
			studyRange.lower = std::numeric_limits<double>::max();
			studyRange.upper = std::numeric_limits<double>::min();

			responseRange.lower = std::numeric_limits<double>::max();
			responseRange.upper = std::numeric_limits<double>::min();
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

	/*
	struct CategoryParameters {
		vector<double> mu;
		double selectivity;
		double kappa;
	};

	struct BaseParameters {
		double pMem;
		double contKappa;
		double pCatGuess;

		CategoryParameters cat;
	};

	struct bwParameters : public BaseParameters {
		double pBetween;

		double pContWithin;
		double pContBetween;
	};
	*/

	struct zlParameters {
		double pMem;
		double contSD;
	};

	struct Parameters {

		double pMem;
		double pBetween;
		double pContBetween;
		double pContWithin;

		double pCatGuess;

		double contSD;

		struct {
			vector<double> mu;
			double selectivity;
			double SD;
		} cat;
	};

	struct ConditionParameters {
		double pMem;
		double pBetween;
		double pContBetween;
		double pContWithin;

		double pCatGuess;

		double contSD;

		struct {
			double selectivity;
			double SD;
		} cat;
	};

	typedef Parameters ParticipantParameters;
	typedef Parameters CombinedParameters;


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