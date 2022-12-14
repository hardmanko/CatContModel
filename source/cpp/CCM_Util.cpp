#include "CCM_Util.h"

#include "CCM_Linear.h"
#include "CCM_Circular.h"

#ifdef COMPILING_WITH_RCPP
#define R_TRUE 1
#define R_FALSE 0
#endif


namespace CatCont {

ModelVariant modelVariantFromString(string modelVariantStr) {
	ModelVariant rval = ModelVariant::BetweenItem;

	if (modelVariantStr == "betweenAndWithin") {
		rval = ModelVariant::BetweenAndWithin;
	} else if (modelVariantStr == "betweenItem") {
		rval = ModelVariant::BetweenItem;
	} else if (modelVariantStr == "withinItem") {
		rval = ModelVariant::WithinItem;
	} else if (modelVariantStr == "ZL") {
		rval = ModelVariant::ZL;
	}

	return rval;
}

DataType dataTypeFromString(string dataTypeString) {
	DataType dataType = DataType::Circular;
	if (dataTypeString == "linear") {
		dataType = DataType::Linear;
	} else if (dataTypeString == "circular") {
		dataType = DataType::Circular;
	}
	return dataType;
}

WeightsDistribution weightsDistributionFromString(string weightsDistStr) {
	WeightsDistribution wd = WeightsDistribution::Default;
	if (weightsDistStr == "PlatSpline") {
		wd = WeightsDistribution::PlatSpline;
	}
	return wd;
}

#ifdef USING_CAT_SPLINE
LambdaVariant lambdaVariantFromString(string lambdaVariantStr) {

	LambdaVariant lambdaVariant = LambdaVariant::None;
	if (lambdaVariantStr == "CatWeightSum") {
		lambdaVariant = LambdaVariant::CatWeightSum;
	} else if (lambdaVariantStr == "InverseCatWeightSum") {
		lambdaVariant = LambdaVariant::InverseCatWeightSum;
	} else if (lambdaVariantStr == "Minus1to1") {
		lambdaVariant = LambdaVariant::Minus1to1;
	}
	return lambdaVariant;

}
#endif

#if COMPILING_WITH_RCPP

double normalLL(double x, double mu, double var) {
	return R::dnorm(x, mu, sqrt(var), R_TRUE); //log = TRUE
}

double cauchyLL(double x, double loc, double scale) {
	return R::dcauchy(x, loc, scale, R_TRUE); //log = TRUE
}

double normalDeviate(double x, double sd) {
	return R::rnorm(x, sd);
}

double uniformDeviate(double low, double high) {
	return R::runif(low, high);
}

#else
double normalLL(double x, double mu, double var) {
	return dnorm(x, mu, sqrt(var), Rboolean::TRUE);
}

double cauchyLL(double x, double loc, double scale) {
	return dcauchy(x, loc, scale, Rboolean::TRUE);
}

double normalDeviate(double x, double sd) {
	return rnorm(x, sd);
}

double uniformDeviate(double low, double high) {
	return runif(low, high);
}
#endif


double clamp(double x, double minimum, double maximum) {
	return std::min(std::max(x, minimum), maximum);
}

string extractIndex(string s) {
	size_t start = s.find('[');
	size_t end = s.find(']');
	if (start == string::npos || end == string::npos) {
		return "";
	}
	return s.substr(start + 1, end - start - 1);
}

string extractBaseParameterName(string s) {
	unsigned int baseEnd = s.find_first_of("_.["); //weak heuristic: pMem_cond[1], pMem.mu, pMem[2]
	return s.substr(0, baseEnd);
}

string joinIndexedParam(string name, vector<string> index) {
	string rval = name + "[";
	for (size_t i = 0; i < index.size(); i++) {
		rval += index[i];
		if (i < index.size() - 1) {
			rval += ",";
		}
	}
	rval += "]";

	return rval;
}

// First index is the parameter name
vector<string> splitIndexedParam(string joined) {

	vector<string> rval;
	rval.push_back(extractBaseParameterName(joined));

	string indexStr = extractIndex(joined);
	vector<string> indexParts = splitString(indexStr, ",");
	rval.insert(rval.end(), indexParts.begin(), indexParts.end());

	return rval;
}

vector<string> splitString(string str, string delim) {
	vector<string> rval;

	size_t offset = 0;
	while (offset < str.size()) {
		size_t delimStart = str.find(delim, offset);
		size_t start = offset;
		size_t end = delimStart - 1;
		rval.push_back(str.substr(start, end));
		offset = delimStart + delim.size();
	}

	if (rval.size() == 0) {
		rval.push_back(str);
	}

	return rval;
}

void logMessage(string module, string message, bool endLine) {
#if COMPILING_WITH_CX
	Log.notice(module) << message;
	Log.flush();
#elif COMPILING_WITH_RCPP
	if (module != "") {
		Rcpp::Rcout << "<" << module << "> " << message;
	} else {
		Rcpp::Rcout << message;
	}
	if (endLine) {
		Rcpp::Rcout << std::endl;
	}
#else
	if (module != "") {
		std::cout << "<" << module << "> " << message;
	} else {
		std::cout << message;
	}
	if (endLine) {
		std::cout << std::endl;
	}
#endif
}


vector<ParticipantData> copyParticipantData(vector<string> pnumsCol, vector<string> condsCol, vector<double> studyCol, vector<double> responseCol, DataType dataType, bool verbose) {

	vector<string> uniquePnums = uniqueElements(pnumsCol);
	vector<string> uniqueConditions = uniqueElements(condsCol);

	map<string, unsigned int> trialsPerCondition;
	for (unsigned int i = 0; i < uniqueConditions.size(); i++) {
		trialsPerCondition[uniqueConditions[i]] = 0;
	}

	vector<ParticipantData> data;

	for (unsigned int p = 0; p < uniquePnums.size(); p++) {

		ParticipantData thisPart;
		thisPart.pnum = uniquePnums[p];

		//Find which rows of the data correspond to the given pnum
		vector<unsigned int> pnumRows;
		for (unsigned int row = 0; row < pnumsCol.size(); row++) {
			if (pnumsCol[row] == uniquePnums[p]) {
				pnumRows.push_back(row);
			}
		}

		for (unsigned int c = 0; c < uniqueConditions.size(); c++) {

			//Copy from the secondary data frame to primitive data types in a ConditionData struct
			ConditionData condData;

			condData.condition = uniqueConditions[c];

			for (unsigned int i = 0; i < pnumRows.size(); i++) {

				unsigned int thisRow = pnumRows[i];

				//If the condition on the current row is equal to the cth condition, store the data
				if (condsCol[thisRow] == uniqueConditions[c]) {
					double study = studyCol[thisRow];
					double response = responseCol[thisRow];

					if (dataType == DataType::Circular) {
						study = Circular::degreesToRadians(study); //TODO: This should probably be done by the model...
						response = Circular::degreesToRadians(response);
					}
					condData.study.push_back(study);
					condData.response.push_back(response);
				}
			}

			thisPart.condData.push_back(condData);

			trialsPerCondition[condData.condition] += condData.study.size();

			/*
			if (verbose) {
				std::stringstream ss;
				ss << "Participant " << uniquePnums[p] << ", condition " << uniqueConditions[c] << ": " <<
					_convertToString(condData.study.size()) << " trials found.";
				logMessage("", ss.str());
			}
			*/
		}

		data.push_back(thisPart);

	}

	if (verbose) {
		std::stringstream ss;
		ss << "Participants found: " << uniquePnums.size() << endl;
		ss << "Total trials per condition: " << endl;
		for (unsigned int i = 0; i < uniqueConditions.size(); i++) {
			string cond = uniqueConditions[i];
			ss << cond << ": " << trialsPerCondition[cond] << endl;
		}
		logMessage("", ss.str());
	}

	return data;

}


//Participant level parameters with estimated priors on them. 
//Basically, every parameter except the catMu and catActive.
vector<string> ModelConfiguration::getParamWithHierachicalPriors(void) const {
	vector<string> paramWithHierarchicalPrior;

	paramWithHierarchicalPrior.push_back("pMem");
	paramWithHierarchicalPrior.push_back("pBetween");
	paramWithHierarchicalPrior.push_back("pContBetween");
	paramWithHierarchicalPrior.push_back("pContWithin");
	paramWithHierarchicalPrior.push_back("pCatGuess");

	paramWithHierarchicalPrior.push_back("contSD");
	paramWithHierarchicalPrior.push_back("catSelectivity");
	paramWithHierarchicalPrior.push_back("catSD");

	return paramWithHierarchicalPrior;
}

vector<string> ModelConfiguration::getParamWithAndWithoutConditionEffects(void) const {
	return this->getParamWithHierachicalPriors();
}

vector<string> ModelConfiguration::getParamWithConditionEffects(void) const {
	return this->paramWithConditionEffects;
}

vector<string> ModelConfiguration::getParamWithoutConditionEffects(void) const {

	vector<string> withAndWithout = this->getParamWithAndWithoutConditionEffects();
	vector<string> with = this->getParamWithConditionEffects();

	vector<string> without;
	for (unsigned int i = 0; i < withAndWithout.size(); i++) {

		bool notInWith = std::find(with.begin(), with.end(), withAndWithout[i]) == with.end();
		if (notInWith) {
			without.push_back(withAndWithout[i]);
		}

	}

	return without;
}



} // namespace CatCont
