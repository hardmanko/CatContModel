#include "CF_ModelConfig.h"

#include "CCM_Util.h"

namespace CatCont {
namespace CatFirst {

bool CatFirstConfig::setFromList(Rcpp::List cfgList) {

	// Copy condition effects
	{
		if (!cfgList.containsElementNamed("condEffects")) {
			CatCont::stop("Config list is missing \"condEffects\".", "CatFirstConfig::setFromList");
			return false;
		}
		Rcpp::DataFrame condDF = Rcpp::as<Rcpp::DataFrame>(cfgList["condEffects"]);
		vector<string> requiredCE = { "cond", "pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess", "contSD", "catSD", "catSelectivity" };
		this->condEffects.clear();
		for (const string& ce : requiredCE) {
			if (!condDF.containsElementNamed(ce.data())) {
				CatCont::stop("\"condEffects\" data frame is missing column \"" + ce + "\"", "CatFirstConfig::setFromList");
			}
			Rcpp::CharacterVector elem = condDF[ce];
			vector<string> ceStr(elem.begin(), elem.end());
			this->condEffects[ce] = ceStr;
		}
	}

	// Copy category effects
	/*
	{
		if (!cfgList.containsElementNamed("catEffects")) {
			CatCont::stop("Config list is missing \"catEffects\".", "CatFirstConfig::setFromList");
			return false;
		}
		Rcpp::DataFrame catDF = cfgList["catEffects"];
		vector<string> requiredCatEff = { "pMem", "pBetween", "pContBetween", "pContWithin", "contSD", "catSD" };
		this->catEffects.clear();
		for (const string& ce : requiredCatEff) {
			if (!catDF.containsElementNamed(ce.data())) {
				CatCont::stop("\"catEffects\" data frame is missing parameter \"" + ce + "\"", "CatFirstConfig::setFromList");
			}
			Rcpp::CharacterVector elem = catDF[ce];
			vector<string> ceStr(elem.begin(), elem.end());
			this->catEffects[ce] = ceStr;
		}
	}
	*/

	// Copy category types
	{
		if (!cfgList.containsElementNamed("catTypes")) {
			CatCont::stop("Config list is missing \"catTypes\".", "CatFirstConfig::setFromList");
			return false;
		}
		Rcpp::DataFrame ctdf = Rcpp::as<Rcpp::DataFrame>(cfgList["catTypes"]);
		vector<string> requiredParam = { "type", "pMem", "pBetween", "pContBetween", "pContWithin", "contSD", "catSD" };
		this->catTypes.clear();
		for (const string& pn : requiredParam) {
			if (!ctdf.containsElementNamed(pn.data())) {
				CatCont::stop("\"catTypes\" data frame is missing column \"" + pn + "\"", "CatFirstConfig::setFromList");
			}
			Rcpp::CharacterVector elem = ctdf[pn];
			vector<string> ceStr(elem.begin(), elem.end());
			this->catTypes[pn] = ceStr;
		}
	}

	// Copy shared participant param
	/*
	{
		if (!cfgList.containsElementNamed("sharedParticipantParam")) {
			CatCont::stop("Config list is missing \"sharedParticipantParam\".", "CatFirstConfig::setFromList");
			return false;
		}

		Rcpp::CharacterVector spp = cfgList["sharedParticipantParam"];
		vector<string> sppStr(spp.begin(), spp.end());
		this->sharedParticipantParam.clear();
		for (const string& str : sppStr) {
			this->sharedParticipantParam.insert(str);
		}
	}
	*/

	return true;
}

size_t CatFirstConfig::getCatTypeCornerstoneIndex(void) const {

	size_t rval = std::numeric_limits<size_t>::max();
	if (this->catTypes.find("type") == this->catTypes.end()) {
		CatCont::stop("catTypes does not have a \"type\" column.", "CatFirstConfig::getCatTypeCornerstoneIndex");
	}

	const vector<string>& catTypeNames = this->catTypes.at("type");
	for (size_t i = 0; i < catTypeNames.size(); i++) {
		if (catTypeNames[i] == "CS") {
			rval = i;
			break;
		}
	}
	if (rval == std::numeric_limits<size_t>::max()) {
		CatCont::stop("No cornerstone type found. \"CS\" should be one element in the \"type\" column.", "CatFirstConfig::getCatTypeCornerstoneIndex");
	}

	return rval;
}

bool CF_CompositeConfig::setFromList(Rcpp::List cfgList) {
	bool mainSuccess = this->modelConfig.setFromList(cfgList);

	Rcpp::List privateCfg = cfgList["privateConfig"];

	bool cfSuccess = false;
	if (privateCfg.containsElementNamed("CatFirst")) {
		cfSuccess = this->cfConfig.setFromList(privateCfg["CatFirst"]);
	}
	else {
		CatCont::stop("privateConfig$CatFirst required in config list.");
	}

	bool vpSuccess = true; // not required
	if (privateCfg.containsElementNamed("VPConfig")) {
		vpSuccess = this->vpConfig.setFromList(privateCfg["VPConfig"]);
	}
			
	return mainSuccess && cfSuccess && vpSuccess;
}

} // namespace CatFirst
} // namespace CatCont