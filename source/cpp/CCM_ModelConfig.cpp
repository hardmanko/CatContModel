#include "CCM_ModelConfig.h"

#include "CCM_Util.h"
#include "CCM_Circular.h"

#include "R_Interface.h" // convertListToMap
#include "CCM_DistributionLUTs.h" // VON_MISES_MAX_SD

namespace CatCont {


ModelVariant modelVariantFromString(string modelVariantStr) {
	ModelVariant rval = ModelVariant::BetweenItem;

	if (modelVariantStr == "betweenAndWithin") {
		rval = ModelVariant::BetweenAndWithin;
	}
	else if (modelVariantStr == "betweenItem") {
		rval = ModelVariant::BetweenItem;
	}
	else if (modelVariantStr == "withinItem") {
		rval = ModelVariant::WithinItem;
	}
	else if (modelVariantStr == "ZL") {
		rval = ModelVariant::ZL;
	}

	return rval;
}

DataType dataTypeFromString(string dataTypeString) {
	DataType dataType = DataType::Circular;
	if (dataTypeString == "linear") {
		dataType = DataType::Linear;
	}
	else if (dataTypeString == "circular") {
		dataType = DataType::Circular;
	}
	return dataType;
}






// TODO: These functions seem more like part of an individual model than model-general.
//Participant level parameters with estimated priors on them. 
//Basically, every parameter except the catMu and catActive.
vector<string> ModelConfiguration::getParamWithHierachicalPriors(void) const {
	vector<string> paramWithHierarchicalPrior = { "pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess", "contSD", "catSD", "catSelectivity" };

	/*
	paramWithHierarchicalPrior.push_back("pMem");
	paramWithHierarchicalPrior.push_back("pBetween");
	paramWithHierarchicalPrior.push_back("pContBetween");
	paramWithHierarchicalPrior.push_back("pContWithin");
	paramWithHierarchicalPrior.push_back("pCatGuess");

	paramWithHierarchicalPrior.push_back("contSD");
	paramWithHierarchicalPrior.push_back("catSelectivity");
	paramWithHierarchicalPrior.push_back("catSD");
	*/

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




bool Linear::LinearConfiguration::setFromList(Rcpp::List cfgList) {
	*this = LinearConfiguration();

	Rcpp::NumericVector rr = cfgList["responseRange"];
	Rcpp::NumericVector cmr = cfgList["catMuRange"];

	this->response.lower = rr[0];
	this->response.upper = rr[1];

	this->catMu.lower = cmr[0];
	this->catMu.upper = cmr[1];

	return true;
}

bool RunConfig::setFromList(Rcpp::List cfgList) {

	*this = RunConfig(); // reset

	if (cfgList.containsElementNamed("iterations")) {
		this->iterations = cfgList["iterations"];
	}
	else {
		CatCont::message("Run config must contain \"iterations\"");
		return false;
	}

	if (cfgList.containsElementNamed("iterationsPerStatusUpdate")) {
		this->iterationsPerStatusUpdate = cfgList["iterationsPerStatusUpdate"];
	}

	if (cfgList.containsElementNamed("verbose")) {
		int verbose_int = cfgList["verbose"];
		this->verbose = (verbose_int == 1);
	}

	if (cfgList.containsElementNamed("profileParameterTypes")) {
		int profile_int = cfgList["profileParameterTypes"];
		this->profileParameterTypes = (profile_int == 1);
	}

	return true;
}

bool ModelConfiguration::setFromList(Rcpp::List cfgList) {

	*this = ModelConfiguration(); // reset

	string modelVariantStr = cfgList["modelVariant"];
	this->modelVariant = CatCont::modelVariantFromString(modelVariantStr);

	string dataTypeStr = cfgList["dataType"];
	this->dataType = CatCont::dataTypeFromString(dataTypeStr);

	string ccName = cfgList["cornerstoneConditionName"];
	this->cornerstoneConditionName = ccName;

	this->catMuPriorApproximationPrecision = cfgList["catMuPriorApproximationPrecision"];

	if (this->modelVariant == CatCont::ModelVariant::ZL) {
		this->maxCategories = 0;
	}
	else {
		this->maxCategories = cfgList["maxCategories"];
	}

	if (this->dataType == CatCont::DataType::Linear) {
		this->linearConfiguration.setFromList(cfgList);
	}

	if (cfgList.containsElementNamed("calculateParticipantLikelihoods")) {
		// Going straight to bool doesn't work right
		int calculateParticipantLikelihoods_temp = cfgList["calculateParticipantLikelihoods"];
		this->calculateParticipantLikelihoods = (calculateParticipantLikelihoods_temp == 1);
	}

	// Set ranges
	this->sdRanges.minSd = cfgList["minSD"];
	this->sdRanges.maxSd = VON_MISES_MAX_SD;
	this->sdRanges.minPrecision = CatCont::Circular::sdDeg_to_precRad(this->sdRanges.maxSd);
	this->sdRanges.maxPrecision = CatCont::Circular::sdDeg_to_precRad(this->sdRanges.minSd);


	// Condition Effects
	Rcpp::List conditionEffects = cfgList["conditionEffects"];
	Rcpp::CharacterVector rawNames = conditionEffects.names();
	for (int i = 0; i < rawNames.size(); i++) {
		string name = Rcpp::as<string>(rawNames[i]);
		Rcpp::CharacterVector ce = conditionEffects[name];

		vector<string> cev(ce.size());

		for (int j = 0; j < ce.size(); j++) {
			cev[j] = Rcpp::as<string>(ce[j]);
		}

		this->conditionEffects[name] = cev;
	}

	for (auto it = this->conditionEffects.begin(); it != this->conditionEffects.end(); it++) {
		const vector<string>& ce = it->second;

		bool isNone = ce.size() == 0 || (ce.size() == 1 && ce[0] == "none");

		if (!isNone) {
			this->paramWithConditionEffects.push_back(it->first);
		}
	}

	
	//Rcpp::CharacterVector paramWithConditionEffects = cfgList["parametersWithConditionEffects"];
	//for (unsigned int i = 0; i < (unsigned int)paramWithConditionEffects.size(); i++) {
	//	this->paramWithConditionEffects.push_back((string)paramWithConditionEffects[i]);
	//}
	

	// Overrides are now part of the config list
	if (cfgList.containsElementNamed("priorOverrides")) {
		this->overrides.priors = convertListToMap<double>(cfgList["priorOverrides"]);
	}
	if (cfgList.containsElementNamed("mhTuningOverrides")) {
		this->overrides.mhTunings = convertListToMap<double>(cfgList["mhTuningOverrides"]);
	}
	if (cfgList.containsElementNamed("startingParamValues")) {
		this->overrides.startingValues = convertListToMap<double>(cfgList["startingParamValues"]);
	}
	if (cfgList.containsElementNamed("constantParamValues")) {
		this->overrides.constantValues = convertListToMap<double>(cfgList["constantParamValues"]);
	}

	//if (cfgList.containsElementNamed("equalityConstraints")) {
	//	this->overrides.equalityConstraints = convertListToMap<double>(cfgList["equalityConstraints"]);
	//}

	// Copy shared participant parameters
	if (cfgList.containsElementNamed("sharedParameters")) {
		Rcpp::CharacterVector spp = cfgList["sharedParameters"];
		vector<string> sppStr(spp.begin(), spp.end());
		this->sharedParameters.clear();
		for (const string& str : sppStr) {
			this->sharedParameters.insert(str);
		}
	}
	else {
		// TODO: Re-add warning
		//CatCont::stop("Config list is missing \"sharedParameters\".", "ModelConfiguration::setFromList");
		//return false;
	}

	// privateConfig
	if (cfgList.containsElementNamed("privateConfig")) {
		Rcpp::List pcList = cfgList["privateConfig"];

		// TODO: Where does this go? It doesn't appear to be used anywhere (why not?)
		if (pcList.containsElementNamed("useVonMisesLookupTable")) {
			int vmlut_int = pcList["useVonMisesLookupTable"];
			this->privateConfig.useVonMisesLookupTable = (vmlut_int == 1);
		}

		// Temporarily, these settings are passed to c++ in privateConfig
		/*
		if (pcList.containsElementNamed("catMuShared")) {
			int cms_int = pcList["catMuShared"];
			this->catMuShared = (cms_int == 1);
		}
		if (pcList.containsElementNamed("catActiveShared")) {
			int cas_int = pcList["catActiveShared"];
			this->catActiveShared = (cas_int == 1);
		}
		*/
		if (pcList.containsElementNamed("catActiveHierPrior")) {
			int cahp_int = pcList["catActiveHierPrior"];
			this->catActiveHierPrior = (cahp_int == 1);
		}
		if (pcList.containsElementNamed("catActiveDistancePrior")) {
			int cadp_int = pcList["catActiveDistancePrior"];
			this->catActiveDistancePrior = (cadp_int == 1);
		}

		//
		if (pcList.containsElementNamed("withinItem_contSD_is_withinSD")) {
			int bool_int = pcList["withinItem_contSD_is_withinSD"];
			this->withinItem_contSD_is_withinSD = (bool_int == 1);
		}
	}

	#ifdef USING_CAT_SPLINE
	if (cfgList.containsElementNamed("weightsDistribution")) {
		string weightsDistributionStr = cfgList["weightsDistribution"];
		this->weightsDistribution = CatCont::weightsDistributionFromString(weightsDistributionStr);
	}
	else {
		this->weightsDistribution = CatCont::WeightsDistribution::Default;
	}
	this->lambdaVariant = CatCont::LambdaVariant::None; // TODO: Choose value.
	#endif

	return true;
}



} // namespace CatCont
