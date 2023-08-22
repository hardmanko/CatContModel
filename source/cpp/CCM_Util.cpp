#include "CCM_Util.h"

#include "CCM_Linear.h"
#include "CCM_Circular.h"

#ifdef COMPILING_WITH_RCPP
#define R_TRUE 1
#define R_FALSE 0
#endif


namespace CatCont {


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

double gammaDeviate(double shape, double rate) {
	return R::rgamma(shape, 1 / rate);
}

double curriedBessel_i(double kappa) {
	// nu (order) = 0, exponent scaled = TRUE (2 = TRUE for some reason)
	return R::bessel_i(kappa, 0, 2);
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
double curriedBessel_i(double kappa) {
	// TODO
	CatCont::stop("curriedBessel_i not implemented.");
	return -1;
}
#endif

void conditionalConfigureVMLut(double maxValue, double stepSize, bool message) {

	if (CatCont::vmLut.readyToUseLUT(maxValue, stepSize)) {
		if (message) {
			CatCont::message("Von Mises look up table already set up.");
		}
		return; // DEBUG
	}

	if (message) {
		CatCont::message("Setting up Von Mises look up table.");
	}

	CatCont::VonMisesLUT::Config vmlConfig;
	vmlConfig.useLUT = true;
	vmlConfig.skipRecomputationIfAble = true; // kinda redundant
	vmlConfig.maxKappa = maxValue;
	vmlConfig.stepSize = stepSize;

	vmlConfig.besselFun = CatCont::curriedBessel_i;

	CatCont::vmLut.setup(vmlConfig);
}

double clamp(double x, double minimum, double maximum) {
	return std::min(std::max(x, minimum), maximum);
}


double paramTransform_probability(double latentP) {
	return logitInverse(latentP);
}

// For circular, applies ranges then transforms to precision (kappa)
double paramTransform_sd(double sd, const SDRanges& sdRanges, DataType dataType) {
	if (dataType == DataType::Circular) {
		// Transform from infinite space to bounded space
		sd = CatCont::clamp(sd, sdRanges.minSd, sdRanges.maxSd);

		// Convert to precision
		double kappa = Circular::sdDeg_to_precRad(sd);

		// Clamp again for extra safety. This is not needed.
		//double clampedKappa = clamp(kappa, ranges.minPrecision, ranges.maxPrecision);

		return kappa;
	}
	else if (dataType == DataType::Linear) {
		sd = CatCont::clamp(sd, sdRanges.minSd, sdRanges.maxSd);
		return sd;
	}

	return -1;
}

string catIndexString(size_t index) {
	return toString(index + 1);
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


std::vector<std::string> splitString(std::string str, std::string delim) {
	std::vector<std::string> rval;

	size_t offset = 0;
	if (delim == "") {
		offset = str.size();
	}

	while (offset < str.size()) {
		size_t delimStart = str.find(delim, offset);

		size_t start = offset;
		size_t end = delimStart;

		offset = end + delim.size();

		if (delimStart == std::string::npos) {
			end = str.size();
			offset = str.size();
		}

		rval.push_back(str.substr(start, end - start));
	}

	if (rval.size() == 0) {
		rval.push_back(str);
	}

	return rval;
}


void stop(string message, string context) {
#if COMPILING_WITH_RCPP
	if (context != "") {
		message = "<" + context + "> " + message;
	}
	Rcpp::stop(message);
#else
	CX::Instances::Log.error(context) << message;
	std::abort();
#endif
}


void message(string message, string context, bool endLine) {
#if COMPILING_WITH_CX
	Log.notice(context) << message;
	Log.flush();
#elif COMPILING_WITH_RCPP
	if (context != "") {
		Rcpp::Rcout << "<" << context << "> " << message;
	} else {
		Rcpp::Rcout << message;
	}
	if (endLine) {
		Rcpp::Rcout << std::endl;
	}
#else
	if (context != "") {
		std::cout << "<" << context << "> " << message;
	} else {
		std::cout << message;
	}
	if (endLine) {
		std::cout << std::endl;
	}
#endif
}


} // namespace CatCont
