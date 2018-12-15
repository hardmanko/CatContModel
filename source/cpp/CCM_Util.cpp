#include "CCM_Util.h"

#include "CCM_Circular.h"

#ifdef COMPILING_WITH_RCPP
#define R_TRUE 1
#define R_FALSE 0
#endif


namespace CatCont {

	VonMisesLut vmLut;


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


	double cubicSplineDensity(double x, double scale) {

		x /= scale;

		if (x < 0) {
			return 1;
		} else if (x > 2) {
			return 0;
		}

		double coef[4] = { 1.0, 0.0, -0.75, 0.25 };
		double dens = 0;
		for (unsigned int i = 0; i < 4; i++) {
			dens += pow(x, (double)i) * coef[i];
		}

		return dens;
	}

	/*
	double dPlateauSpline(double x, double mu, double plateauHalfWidth, double tailScale) {
		double d = circularAbsoluteDistance(x, mu); // TODO: This is wrong for linear
		return dPlateauSpline(d, dPlateauSpline, tailScale);
	}
	*/

	// d is distance from center of plateau.
	// plateauHalfWidth should be chosen depending on which side of the center you are (not this function's job).
	double dPlateauSpline(double d, double plateauHalfWidth, double tailScale) {
		double remD = max(d - plateauHalfWidth, 0.0); // You don't need to do this max
		return cubicSplineDensity(remD, tailScale);
	}


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



}