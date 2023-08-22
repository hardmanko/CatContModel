#pragma once

#include "CCM_Main.h"

#include "CCM_ModelConfig.h"

#include "UtilityFunctions.h"

#if COMPILING_WITH_CX
#include "CX.h"
#include "ofxRMath.h"
#endif




#include <vector>
#include <map>
#include <string>

using namespace std;

namespace CatCont {

	// Should these be in Main?
	void stop(string message, string context = "");
	void message(string message, string context = "", bool endLine = true);


	double paramTransform_probability(double latentP); // = logitInverse
	double paramTransform_sd(double sd, const SDRanges& sdRanges, DataType dataType);

	
	
	double normalLL(double x, double mu, double var);
	double cauchyLL(double x, double loc, double scale);
	double normalDeviate(double x, double sd);
	double uniformDeviate(double low, double high);
	double gammaDeviate(double shape, double rate);

	// TODO: Move to DistributionLUTs? Or keep those files pure?
	double curriedBessel_i(double kappa);
	void conditionalConfigureVMLut(double maxValue, double stepSize, bool message = true);

	double clamp(double x, double minimum, double maximum);

	string catIndexString(size_t index);

	string extractIndex(string s);
	string extractBaseParameterName(string s);

	string joinIndexedParam(string name, vector<string> index);
	vector<string> splitIndexedParam(string joined);
	vector<string> splitString(string str, string delim);

	template <typename T>
	string toString(T val) {
		ostringstream oss;
		oss << val;
		return oss.str();
	}

	template <typename T>
	T fromString(string s) {
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


	// TODO: Move elsewhere
#ifdef COMPILING_WITH_CX
	void useModel_Bayesian(void); // Freestanding interface
#endif

}