#include "CCM_Util.h"

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


	double clamp(double x, double minimum, double maximum) {
		return std::min(std::max(x, minimum), maximum);
	}

	void logMessage(string module, string message, bool endLine) {
#if COMPILING_WITH_CX
		Log.notice(module) << message;
		Log.flush();
#elif COMPILING_WITH_RCPP
		if (module != "") {
			Rcpp::Rcout << "<" << module << "> " << message << std::endl;
		} else {
			Rcpp::Rcout << message << std::endl;
		}
		if (endLine) {
			Rcpp::Rcout << std::endl;
		}
#else
		if (module != "") {
			std::cout << "<" << module << "> " << message << std::endl;
		} else {
			std::cout << message << std::endl;
		}
		if (endLine) {
			std::cout << std::endl;
		}
#endif
	}





}