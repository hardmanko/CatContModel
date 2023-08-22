#pragma once

// Included from Compilation.h
//#include "Rcpp.h"

#include "CCM_Main.h"

#include "CCM_Util.h"

namespace CatCont {
	//vector<ParticipantData> getParticipantData(Rcpp::DataFrame df, CatCont::DataType dataType, bool verbose)
}

vector<ParameterList> transposePosteriors(Rcpp::List posteriors);

//CatCont::ModelConfiguration readModelConfigFromList(Rcpp::List configList);

template <typename T>
std::map<std::string, T> convertListToMap(Rcpp::List list) {
	std::map<std::string, T> rval;

	if (list.size() == 0) {
		return rval;
	}

	Rcpp::CharacterVector rawNames = list.names();
	std::vector<string> names(rawNames.begin(), rawNames.end());
	for (size_t i = 0; i < names.size(); i++) {
		rval[names[i]] = Rcpp::as<T>(list[names[i]]);
	}

	return rval;
}

std::vector<std::string> convertCharVector(const Rcpp::CharacterVector& cv);
std::vector<double> convertNumericVector(const Rcpp::NumericVector& nv);