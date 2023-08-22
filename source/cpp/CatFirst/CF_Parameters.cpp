#include "CF_Parameters.h"

#include "R_Interface.h"


// [[Rcpp::export]]
Rcpp::List CCM_CF_testParameterStuff(Rcpp::List args) {

	using namespace CatCont;
	using namespace CatCont::CatFirst;

	Rcpp::List rval;

	StandardParam sp;
	Rcpp::List spValues = args["standardParam"];

	ParamMap paramList = convertListToMap<double>(spValues);

	sp.setFromList(paramList, "");

	rval["standardParamValues"] = sp.toMap();

	SDRanges sdRanges;
	sdRanges.minSd = 1;
	sdRanges.maxSd = 1000;
	sp.transformToManifest(sdRanges, DataType::Circular);

	rval["standardParamTransformed"] = sp.toMap();

	return rval;
}


// [[Rcpp::export]]
Rcpp::List CCM_CPP_CF_getMPMP(Rcpp::List modCfgList, Rcpp::List partList, Rcpp::List condEffList, Rcpp::DataFrame catTypeEffectsDF) {

	// partList has only active categories.
	// Single condition only.
	// catTypes???

	using namespace Rcpp;
	using namespace CatCont;

	ModelConfiguration modelConfig;
	modelConfig.setFromList(modCfgList);

	//
	CatFirst::StandardParam partStd;

	partStd.pMem = partList["pMem"];
	partStd.pBetween = partList["pBetween"];
	partStd.pContBetween = partList["pContBetween"];
	partStd.pContWithin = partList["pContWithin"];
	partStd.pCatGuess = partList["pCatGuess"];

	partStd.contSD = partList["contSD"];
	partStd.catSD = partList["catSD"];
	partStd.catSelectivity = partList["catSelectivity"];

	//
	CatFirst::StandardParam condStd;

	condStd.pMem = condEffList["pMem"];
	condStd.pBetween = condEffList["pBetween"];
	condStd.pContBetween = condEffList["pContBetween"];
	condStd.pContWithin = condEffList["pContWithin"];
	condStd.pCatGuess = condEffList["pCatGuess"];

	condStd.contSD = condEffList["contSD"];
	condStd.catSD = condEffList["catSD"];
	condStd.catSelectivity = condEffList["catSelectivity"];

	//
	CatFirst::CategoryParam partCat;

	vector<double> catMu = Rcpp::as<vector<double>>(partList["catMu"]);
	vector<double> catActive = Rcpp::as<vector<double>>(partList["catActive"]);
	vector<double> catType = Rcpp::as<vector<double>>(partList["catType"]);

	partCat.setFromVectors(catMu, catActive, catType);

	size_t nCatActive = partCat.nCatActive();



	// Cat Type Effects. Set per catType
	CatFirst::CatTypeEffects catEff;

	size_t nCatTypes = catTypeEffectsDF.nrows();
	//catEff.typeNames(nCatTypes);
	catEff.catTypeEffects.resize(nCatTypes);

	catEff.typeNames = Rcpp::as<vector<string>>(catTypeEffectsDF["type"]);

	vector<double> cte_pMem = Rcpp::as<vector<double>>(catTypeEffectsDF["pMem"]);
	//vector<double> cte_pMem = convertNumericVector(catTypeEffectsDF["pMem"]);
	vector<double> cte_pBetween = convertNumericVector(catTypeEffectsDF["pBetween"]);
	vector<double> cte_pContBetween = convertNumericVector(catTypeEffectsDF["pContBetween"]);
	vector<double> cte_pContWithin = convertNumericVector(catTypeEffectsDF["pContWithin"]);
	vector<double> cte_contSD = convertNumericVector(catTypeEffectsDF["contSD"]);
	vector<double> cte_catSD = convertNumericVector(catTypeEffectsDF["catSD"]);

	for (size_t t = 0; t < nCatTypes; t++) {
		CatFirst::SingleCategoryEffect& cak = catEff.catTypeEffects[t];

		cak.pMem = cte_pMem[t];
		cak.pBetween = cte_pBetween[t];
		cak.pContBetween = cte_pContBetween[t];
		cak.pContWithin = cte_pContWithin[t];
		cak.contSD = cte_contSD[t];
		cak.catSD = cte_catSD[t];
	}


	//
	CatFirst::SpecialParam special_empty;
	special_empty.pCatUsedToGuess.resize(nCatActive, 0);

	// Make and set the mpmp
	CatFirst::MPMP_complete mpmp;
	mpmp.setFromLatent(modelConfig, partCat, partStd, condStd, catEff, special_empty);

	// DEBUG Print to console
	//CatCont::message(mpmp.printValues());

	// Copy mpmp into the rval
	Rcpp::List rval;

	rval["categorization"] = List::create(
		Named("catSelectivity") = wrap(mpmp.categorization.catSelectivity),
		Named("catMu") = wrap(mpmp.categorization.cats.catMu),
		Named("catActive") = wrap(mpmp.categorization.cats.catActive),
		Named("catType") = wrap(mpmp.categorization.cats.catType)
	);

	rval["guessing"] = List::create(
		Named("pCatGuess") = mpmp.guessing.pCatGuess,
		Named("catSD") = wrap(mpmp.guessing.catSD),
		Named("pCatUsedToGuess") = wrap(mpmp.guessing.pCatUsedToGuess)
	);

	List memList;
	for (size_t t = 0; t < mpmp.memory.size(); t++) {
		string typeName = catEff.typeNames[t];
		const CatFirst::MPMP_memory& memPar = mpmp.memory[t];

		memList[typeName] = List::create(
			Named("pMem") = memPar.pMem,
			Named("pBetween") = memPar.pBetween,
			Named("pContBetween") = memPar.pContBetween,
			Named("pContWithin") = memPar.pContWithin,
			Named("contSD") = memPar.contSD,
			Named("catSD") = memPar.catSD
		);
	}
	rval["memory"] = memList;

	return rval;
}

namespace CatCont {
namespace CatFirst {

double extractParam_standard(const ParameterList& param, string paramBaseName, string pnum, string cond, int cat) {
	double partV = 0;
	double condV = 0;
	double catV = 0;
	if (pnum != "") {
		partV = param.at(paramBaseName + "_part[" + pnum + "]");
	}
	if (cond != "") {
		condV = param.at(paramBaseName + "_cond[" + cond + "]");
	}
	if (cat >= 0) {
		catV = param.at(paramBaseName + "_catClass[" + CatCont::catIndexString(cat) + "]");
	}
	return partV + condV + catV;
}

double extractParam_category(const ParameterList& param, string paramBaseName, string pnum, int cat) {
	string indexStr = paramBaseName + "_part[" + pnum + "," + CatCont::catIndexString(cat) + "]";
	return param.at(indexStr);
}

  /////////////////////////
 // SingleStandardParam //
/////////////////////////

void SingleStandardParam::setFromList(const ParameterList& param, string paramBaseName, string pnum, string cond, int cat) {
	double partV = 0;
	double condV = 0;
	double catV = 0;
	if (pnum != "") {
		partV = param.at(paramBaseName + "_part[" + pnum + "]");
	}
	if (cond != "") {
		condV = param.at(paramBaseName + "_cond[" + cond + "]");
	}
	if (cat >= 0) {
		catV = param.at(paramBaseName + "_catClass[" + CatCont::catIndexString(cat) + "]");
	}
	this->value = partV + condV + catV;
}

void SingleStandardParam::toManifestProb(void) {
	this->value = CatCont::paramTransform_probability(this->value);
}

void SingleStandardParam::toManifestSD(const SDRanges& sdRanges, DataType dataType) {
	this->value = CatCont::paramTransform_sd(this->value, sdRanges, dataType);
}

  ///////////////////
 // StandardParam //
///////////////////

// indexStr can be like _pnum[123] or _cond[C2] or _catClass[5]
bool StandardParam::setFromList(const ParameterList& param, string indexStr) {
	this->pMem = param.at("pMem" + indexStr);
	this->pBetween = param.at("pBetween" + indexStr);
	this->pContBetween = param.at("pContBetween" + indexStr);
	this->pContWithin = param.at("pContWithin" + indexStr);
	this->pCatGuess = param.at("pCatGuess" + indexStr);

	this->contSD = param.at("contSD" + indexStr);
	this->catSD = param.at("catSD" + indexStr);
	this->catSelectivity = param.at("catSelectivity" + indexStr);

	return true;
}

map<string, double> StandardParam::toMap(void) const {
	map<string, double> rval;

	rval["pMem"] = this->pMem;
	rval["pBetween"] = this->pBetween;
	rval["pContBetween"] = this->pContBetween;
	rval["pContWithin"] = this->pContWithin;
	rval["pCatGuess"] = this->pCatGuess;

	rval["contSD"] = this->contSD;
	rval["catSD"] = this->catSD;
	rval["catSelectivity"] = this->catSelectivity;

	return rval;
}

void StandardParam::transformToManifest(const SDRanges& sdRanges, DataType dataType) {
	this->pMem = CatCont::paramTransform_probability(this->pMem);
	this->pBetween = CatCont::paramTransform_probability(this->pBetween);
	this->pContBetween = CatCont::paramTransform_probability(this->pContBetween);
	this->pContWithin = CatCont::paramTransform_probability(this->pContWithin);
	this->pCatGuess = CatCont::paramTransform_probability(this->pCatGuess);

	this->contSD = CatCont::paramTransform_sd(this->contSD, sdRanges, dataType);
	this->catSD = CatCont::paramTransform_sd(this->catSD, sdRanges, dataType);
	this->catSelectivity = CatCont::paramTransform_sd(this->catSelectivity, sdRanges, dataType);
}

StandardParam StandardParam::operator+(const StandardParam& rhs) const {
	StandardParam rval = *this;

	rval.pMem += rhs.pMem;
	rval.pBetween += rhs.pBetween;
	rval.pContBetween += rhs.pContBetween;
	rval.pContWithin += rhs.pContWithin;
	rval.pCatGuess += rhs.pCatGuess;

	rval.contSD += rhs.contSD;
	rval.catSD += rhs.catSD;
	rval.catSelectivity += rhs.catSelectivity;

	return rval;
}

StandardParam StandardParam::operator+(const SingleCategoryEffect& rhs) const {
	StandardParam rval = *this;

	rval.pMem += rhs.pMem;
	rval.pBetween += rhs.pBetween;
	rval.pContBetween += rhs.pContBetween;
	rval.pContWithin += rhs.pContWithin;
	//rval.pCatGuess += rhs.pCatGuess;

	rval.contSD += rhs.contSD;
	rval.catSD += rhs.catSD;
	//rval.catSelectivity += rhs.catSelectivity;

	return rval;
}

///////////////////////////////////
/*
void SingleCategoryEffect::setFromList(const ParameterList& param, size_t k) {

	string indexStr = "_cat[" + CatCont::catIndexString(k) + "]";

	this->pMem = param.at("pMem" + indexStr);
	this->pBetween = param.at("pBetween" + indexStr);
	this->pContBetween = param.at("pContBetween" + indexStr);
	this->pContWithin = param.at("pContWithin" + indexStr);
	//this->pCatGuess = param.at("pCatGuess" + indexStr);

	this->contSD = param.at("contSD" + indexStr);
	this->catSD = param.at("catSD" + indexStr);
	//this->catSelectivity = param.at("catSelectivity" + indexStr);

	//this->pCatUsedToGuess = param.at("pCatUsedToGuess" + indexStr);

}

void CompleteCategoryEffects::setFromList(const ParameterList& param, const vector<size_t>& activeCatK) {

	this->paramPerAK.resize(activeCatK.size());

	for (size_t i = 0; i < activeCatK.size(); i++) {
		size_t ak = activeCatK[i];
		this->paramPerAK[i].setFromList(param, ak);
	}
}

bool SpecialParam::setFromList(const ParameterList& param, string pnum, size_t maxCategories) {

	string nameBase = "pCatUsedToGuess_part[" + pnum + ",";
	this->pCatUsedToGuess.clear();
	for (size_t k = 0; k < maxCategories; k++) {
		pCatUsedToGuess.push_back(param.at(nameBase + CatCont::catIndexString(k) + "]"));
	}

	catSD_guessShift = param.at("catSD_guessShift_part[" + pnum + "]");

	return true;
}

void SpecialParam::setFromList(const ParameterList& param, string pnum, const vector<size_t>& aks) {

	string nameBase = "pCatUsedToGuess_part[" + pnum + ",";
	this->pCatUsedToGuess.clear();
	for (size_t ak : aks) {
		//for (size_t k = 0; k < maxCategories; k++) {
		pCatUsedToGuess.push_back(param.at(nameBase + CatCont::catIndexString(ak) + "]"));
	}

	catSD_guessShift = param.at("catSD_guessShift_part[" + pnum + "]");
}
*/

  ///////////////////
 // CategoryParam //
///////////////////

vector<double> CategoryParam::getActiveCatMu(void) const {
	vector<double> acm;
	for (size_t k = 0; k < this->catActive.size(); k++) {
		if (this->catActive[k] == 1) {
			acm.push_back(this->catMu[k]);
		}
	}
	return acm;
}

vector<size_t> CategoryParam::getActiveCatK(void) const {
	vector<size_t> ack;
	for (size_t k = 0; k < this->catActive.size(); k++) {
		if (this->catActive[k] == 1) {
			ack.push_back(k);
		}
	}
	return ack;
}

vector<size_t> CategoryParam::getActiveCatType(void) const {
	vector<size_t> act;
	for (size_t k = 0; k < this->catActive.size(); k++) {
		if (this->catActive[k] == 1) {
			act.push_back(this->catType[k]);
		}
	}
	return act;
}

ActiveCategoryParam CategoryParam::getActiveParam(void) const {
	ActiveCategoryParam rval;

	for (size_t k = 0; k < this->catActive.size(); k++) {
		if (this->catActive[k] == 1) {
			rval.k.push_back(k);
			rval.catMu.push_back(this->catMu[k]);
			rval.catType.push_back(this->catType[k]);
		}
	}

	return rval;
}

int CategoryParam::nCatActive(int version) const {
	if (version == 1) {
		int sum = 0;
		for (size_t i = 0; i < this->catActive.size(); i++) {
			sum += this->catActive[i];
		}
		return sum;
	} 
	//else if (version == 2) {
	//	return this->activeCatMu.size();
	//}
	return 0;
}

void CategoryParam::clear(void) {
	this->catMu.clear();
	this->catActive.clear();
	this->catType.clear();

	//this->activeCatMu.clear();
	//this->activeCatK.clear();
}

/*/
bool CategoryParam::setFromList(const ParameterList& param, string pnum, size_t maxCategories) {

	this->clear();

	string catIStrBase = "_part[" + pnum + ",";

	for (size_t k = 0; k < maxCategories; k++) {
		string catIstr = catIStrBase + CatCont::catIndexString(k) + "]";

		double catMu = param.at("catMu" + catIstr);
		unsigned int catActive = param.at("catActive" + catIstr);

		this->catMu.push_back(catMu);
		this->catActive.push_back(catActive);

		if (catActive == 1) {
			this->activeCatMu.push_back(catMu);
			this->activeCatK.push_back(k);
		}
		//if not active, don't add to the active lists
	}

	return true;
}
*/

bool CategoryParam::setFromVectors(vector<double> catMu_, vector<double> catActive_, vector<double> catType_) {
	this->clear();

	if (catActive_.size() == 0) {
		catActive_.resize(catMu_.size(), 1);
	}

	if (catActive_.size() != catMu_.size() || catType_.size() != catMu_.size()) {
		return false;
	}

	for (size_t k = 0; k < catMu_.size(); k++) {

		this->catMu.push_back(catMu_[k]);
		this->catActive.push_back((unsigned int)catActive_[k]);
		this->catType.push_back((size_t)catType_[k]);

		/*
		if (catActiveInt == 1) {
			this->activeCatMu.push_back(catMu_[k]);
			this->activeCatK.push_back(k);
		}
		*/
		//if not active, don't add to the active lists
	}

	return true;
}

  //////////
 // MPMP //
//////////

void MPMP_categorization::setFromLatent(const ModelConfiguration& modelConfig, const CategoryParam& catPar, double catSelectivity_) {
	this->cats = catPar;

	this->cats.catMu = CatCont::Circular::degreesToRadians(this->cats.catMu);

	this->catSelectivity = CatCont::paramTransform_sd(catSelectivity_, modelConfig.sdRanges, modelConfig.dataType);
}

void MPMP_memory::setFromLatent(const ModelConfiguration& modelConfig, StandardParam standard) {

	standard.transformToManifest(modelConfig.sdRanges, modelConfig.dataType);

	this->pMem = standard.pMem;
	this->pBetween = standard.pBetween;
	this->pContBetween = standard.pContBetween;
	this->pContWithin = standard.pContWithin;

	this->contSD = standard.contSD;
	this->catSD = standard.catSD;
	//this->catSelectivity = standard.catSelectivity;

	// Model variant parameter constraints
	if (modelConfig.modelVariant == ModelVariant::BetweenItem) {
		this->pBetween = 1;
	}
	else if (modelConfig.modelVariant == ModelVariant::WithinItem) {
		this->pBetween = 0;
	}
	else if (modelConfig.modelVariant == ModelVariant::ZL) {
		this->pBetween = 1;
		this->pContBetween = 1;
	}

}

void MPMP_guessing::setFromLatent(const ModelConfiguration& modelConfig, double pCatGuess_, const vector<double>& catSD_, const vector<double>& pCatUsedToGuess_) {

	this->pCatGuess = CatCont::paramTransform_probability(pCatGuess_);

	size_t nCatActive = catSD_.size();
	this->catSD.resize(nCatActive);
	this->pCatUsedToGuess.resize(nCatActive);

	for (size_t i = 0; i < nCatActive; i++) {
		this->catSD[i] = CatCont::paramTransform_sd(catSD_.at(i), modelConfig.sdRanges, modelConfig.dataType);
		this->pCatUsedToGuess[i] = CatCont::paramTransform_probability(pCatUsedToGuess_.at(i));
	}

	// model variant specific constraints
	if (modelConfig.modelVariant == ModelVariant::ZL) {
		this->pCatGuess = 0;
	}
}

  ///////////////////
 // MPMP_complete //
///////////////////



void MPMP_complete::setFromLatent(const ModelConfiguration& modelConfig, const CategoryParam& partCat, const StandardParam& partStd, const StandardParam& condStd, 
	const CatTypeEffects& catEff, const SpecialParam& special)
{

	// Combine standard parameters for part and cond
	StandardParam standard_ij = partStd + condStd;

	// Set categorization param
	this->categorization.setFromLatent(modelConfig, partCat, standard_ij.catSelectivity);

	vector<size_t> activeCatType = this->categorization.cats.getActiveCatType();
	size_t nCatActive = activeCatType.size(); // this->categorization.cats.nCatActive();

	
	size_t nCatTypes = catEff.catTypeEffects.size();

	// Make place to store guessing parameters (latent)
	vector<double> catSD_t(nCatTypes);

	this->memory.resize(nCatTypes);
	for (size_t t = 0; t < nCatTypes; t++) {

		const SingleCategoryEffect& thisCatEff_t = catEff.catTypeEffects.at(t);

		// Combine standard ij with cat effects
		StandardParam standard_ijt = standard_ij + thisCatEff_t; // i'm the operator with my pocket calculator

		// Copy values for guessing vectors
		catSD_t[t] = standard_ijt.catSD; // + special.catSD_guessShift;

		// Copy to MPMP_memory
		this->memory[t].setFromLatent(modelConfig, standard_ijt);
	}

	// guessing
	vector<double> guessing_catSD_ak(nCatActive);
	vector<double> guessing_pCatUsedToGuess_ak(nCatActive);
	for (size_t ak = 0; ak < nCatActive; ak++) {
		size_t t = activeCatType.at(ak);
		guessing_catSD_ak[ak] = catSD_t[t] + special.catSD_guessShift;

		guessing_pCatUsedToGuess_ak[ak] = special.pCatUsedToGuess.at(ak);
	}

	/*
	if (nCatActive == 0) {
		// For 0 active cats, resize memory to 1 and set based on part/cond parameters.
		this->memory.resize(1);

		this->memory.front().setFromLatent(modelConfig, standard_ij);

		// Special parameters are only related to categories
	}
	else {
	
		// 1 or more active categories
		this->memory.resize(nCatActive);
		for (size_t ak = 0; ak < nCatActive; ak++) {

			size_t catTypeIndex = this->categorization.cats.catType.at(ak);

			const SingleCategoryEffect& thisCatEff = catEff.catTypeEffects.at(catTypeIndex);

			// Combine standard ij with cat effects
			StandardParam standard_ijk = standard_ij + thisCatEff; // i'm the operator with my pocket calculator


			// Copy values for guessing vectors
			guessing_catSD[ak] = standard_ijk.catSD + special.catSD_guessShift;
			guessing_pCatUsedToGuess[ak] = special.pCatUsedToGuess.at(ak);


			// Copy to MPMP_memory
			this->memory[ak].setFromLatent(modelConfig, standard_ijk);

		}
	}
	*/

	this->guessing.setFromLatent(modelConfig, standard_ij.pCatGuess, guessing_catSD_ak, guessing_pCatUsedToGuess_ak);

}

/*
void MPMP_complete::setFromList(const ModelConfiguration& modelConfig, const ParameterList& param, string pnum, string cond) {

	CategoryParam partCat;
	partCat.setFromList(param, pnum, modelConfig.maxCategories);

	StandardParam partStd;
	partStd.setFromList(param, "_part[" + pnum + "]");

	StandardParam condStd;
	condStd.setFromList(param, "_cond[" + cond + "]");

	SpecialParam special;
	CompleteCategoryEffects catEff;
	if (partCat.nCatActive() > 0) {
		special.setFromList(param, pnum, this->categorization.cats.activeCatK);
		catEff.setFromList(param, this->categorization.cats.activeCatK);
	}

	this->setFromLatent(modelConfig, partCat, partStd, condStd, catEff, special);
}
*/

void MPMP_complete::setFromContainer(const ModelConfiguration& modelConfig, const ParamExtractor& pec, const ParamContainer& pc, size_t partInd, size_t condInd) {

	CategoryParam partCat = pec.getParticipantCategory(pc, partInd);
	vector<size_t> activeCatK = partCat.getActiveCatK();

	StandardParam partStd = pec.getParticipantStandard(pc, partInd);

	StandardParam condStd = pec.getConditionStandard(pc, condInd);

	SpecialParam special;
	//CompleteCategoryEffects catEff;
	CatTypeEffects catTypeEff;
	if (partCat.nCatActive() > 0) {
		special = pec.getParticipantSpecial(pc, partInd, activeCatK);
		//catEff = pec.getCategoryEffects(pc, partCat.activeCatK);
		//catEff = pec.getCategoryTypeEffects(pc, partInd, activeCatK);
		catTypeEff = pec.getCategoryTypeEffects_v2(pc);
	}

	this->setFromLatent(modelConfig, partCat, partStd, condStd, catTypeEff, special);

}

std::string MPMP_complete::printValues(void) const {

	std::stringstream ss;

	// Categorization
	ss << "catSelectivity = " << this->categorization.catSelectivity << std::endl;

	const CategoryParam& catPar = this->categorization.cats;

	ss << "catActive[k]   catType[k]   catMu[k]" << std::endl;
	for (size_t k = 0; k < this->categorization.cats.catMu.size(); k++) {
		ss << catPar.catActive[k] << "   " << catPar.catType[k] << "   " << catPar.catMu[k] << std::endl;
	}
	
	ActiveCategoryParam acp = this->categorization.cats.getActiveParam();

	ss << "active: ak   catType[ak]   catMu[ak]" << std::endl;
	for (size_t ak = 0; ak < acp.k.size(); ak++) {
		ss << ak << "   " << acp.catType[ak] << "   " << acp.catMu[ak] << std::endl;
	}

	// Memory
	for (size_t ak = 0; ak < this->memory.size(); ak++) {
		ss << "memory[" << ak << "]" << std::endl;

		const MPMP_memory& mem = this->memory[ak];

		ss << "pMem: " << mem.pMem << std::endl;
		ss << "pBetween: " << mem.pBetween << std::endl;
		ss << "pContBetween: " << mem.pContBetween << std::endl;
		ss << "pContWithin: " << mem.pContWithin << std::endl;
		ss << "contSD: " << mem.contSD << std::endl;
		ss << "catSD: " << mem.catSD << std::endl;
		//ss << "catSelectivity: " << mem.catSelectivity << std::endl;

	}

	// Guessing
	ss << "guessing" << std::endl;
	const MPMP_guessing& gs = this->guessing;

	ss << "pCatGuess: " << gs.pCatGuess << std::endl;
	for (size_t ak = 0; ak < gs.catSD.size(); ak++) {
		ss << "catSD: " << gs.catSD[ak] << std::endl;
		ss << "pCatUsedToGuess: " << gs.pCatUsedToGuess[ak] << std::endl;
	}

	return ss.str();
}

  //////////////////////////////
 // ParamExtractor //
//////////////////////////////

// Set the proxy indices based on pc
void ParamExtractor::setup(const CF_CompositeConfig& cfConfig, const DataCollection& data, const ParamContainer& pc) {

	this->_cfConfig = cfConfig;
	const unsigned int& maxCat = cfConfig.modelConfig.maxCategories;

	//CatCont::message("MaxCat=" + CatCont::toString(maxCat), "ParamExtractor::setup"); // DEBUG

	// Participants
	this->part.resize(data.participants.size());
	for (size_t i = 0; i < data.participants.size(); i++) {
		PI_Participant& pref = this->part[i];

		pref.pnum = data.participants[i].pnum;

		// Standard
		string partIS = "_part[" + pref.pnum + "]";
		pref.standard = _setStandardProxies(pc, partIS);

		// Category
		pref.cat.catMu.resize(maxCat);
		pref.cat.catActive.resize(maxCat);
		pref.cat.catType.resize(maxCat);

		// and Special
		pref.special.catSD_guessShift = _makeProxy(pc, "catSD_guessShift" + partIS);
		pref.special.pCatUsedToGuess.resize(maxCat);

		for (size_t k = 0; k < maxCat; k++) {
			string catIS = "_part[" + pref.pnum + "," + CatCont::catIndexString(k) + "]";

			pref.cat.catMu[k] = _makeProxy(pc, "catMu" + catIS);
			pref.cat.catActive[k] = _makeProxy(pc, "catActive" + catIS);
			pref.cat.catType[k] = _makeProxy(pc, "catType" + catIS);

			pref.special.pCatUsedToGuess[k] = _makeProxy(pc, "pCatUsedToGuess" + catIS);
		}
	}

	// Conditions
	this->cond.resize(data.conditionNames.size());
	for (size_t j = 0; j < data.conditionNames.size(); j++) {
		this->cond[j].condName = data.conditionNames[j];

		string condIS = "_cond[" + data.conditionNames[j] + "]";
		this->cond[j].standard = _setStandardProxies(pc, condIS);
	}

	// CatEffects
	/*
	this->catEffects.resize(modelConfig.maxCategories);
	for (size_t k = 0; k < modelConfig.maxCategories; k++) {

		PI_CatEffect& picck = this->catEffects[k];
		picck.catIndexStr = CatCont::catIndexString(k);

		string ccIndexStr = "_cat[" + picck.catIndexStr + "]";

		picck.pMem = _makeProxy(pc, "pMem" + ccIndexStr);
		picck.pBetween = _makeProxy(pc, "pBetween" + ccIndexStr);
		picck.pContBetween = _makeProxy(pc, "pContBetween" + ccIndexStr);
		picck.pContWithin = _makeProxy(pc, "pContWithin" + ccIndexStr);
		picck.contSD = _makeProxy(pc, "contSD" + ccIndexStr);
		picck.catSD = _makeProxy(pc, "catSD" + ccIndexStr);
	}
	*/

	// catTypes
	const vector<string> catTypeNames = cfConfig.cfConfig.catTypes.at("type");
	this->catTypes.resize(catTypeNames.size());
	for (size_t ct = 0; ct < catTypeNames.size(); ct++) {

		PI_CatType& pict = this->catTypes[ct];
		pict.typeNameStr = catTypeNames.at(ct);

		string ctIndexStr = "_cat[" + catTypeNames.at(ct) + "]";

		pict.pMem = _makeProxy(pc, "pMem" + ctIndexStr);
		pict.pBetween = _makeProxy(pc, "pBetween" + ctIndexStr);
		pict.pContBetween = _makeProxy(pc, "pContBetween" + ctIndexStr);
		pict.pContWithin = _makeProxy(pc, "pContWithin" + ctIndexStr);
		pict.contSD = _makeProxy(pc, "contSD" + ctIndexStr);
		pict.catSD = _makeProxy(pc, "catSD" + ctIndexStr);
	}
}

// This function puts all of the logic of extraction and transformation into one function.
/*
MPMP_complete ParamExtractor::getCompleteMPMP(const ModelConfiguration& modelConfig, const ParamContainer& param, size_t partInd, size_t condInd) const {

	MPMP_complete rval;

	const PI_Standard& ps = this->part.at(partInd).standard;

	const PI_Condition& pc = this->cond.at(condInd).standard;

	double pMem_ij =           param.get(ps.pMem) +           param.get(pc.pMem);
	double pBetween_ij =       param.get(ps.pBetween) +       param.get(pc.pBetween);
	double pContBetween_ij =   param.get(ps.pContBetween) +   param.get(pc.pContBetween);
	double pContWithin_ij =    param.get(ps.pContWithin) +    param.get(pc.pContWithin);
	double pCatGuess_ij =      param.get(ps.pCatGuess) +      param.get(pc.pCatGuess);

	double contSD_ij =         param.get(ps.contSD) +         param.get(pc.contSD);
	double catSD_ij =          param.get(ps.catSD) +          param.get(pc.catSD);
	double catSelectivity_ij = param.get(ps.catSelectivity) + param.get(pc.catSelectivity);

	// TODO


}
*/

// This function is just a little wrapper
MPMP_complete ParamExtractor::getCompleteMPMP(const ParamContainer& param, size_t partInd, size_t condInd) const {

	MPMP_complete rval;

	rval.setFromContainer(_cfConfig.modelConfig, *this, param, partInd, condInd);

	return rval;
}

StandardParam ParamExtractor::getParticipantStandard(const ParamContainer& pc, size_t partInd) const {

	const PI_Standard& ps = this->part.at(partInd).standard;

	return _getStandard(pc, ps);
}

StandardParam ParamExtractor::getConditionStandard(const ParamContainer& pc, size_t condInd) const {
	const PI_Standard& ps = this->cond.at(condInd).standard;

	return _getStandard(pc, ps);
}

CategoryParam ParamExtractor::getParticipantCategory(const ParamContainer& pc, size_t partInd) const {

	const PI_Category& catInd = this->part[partInd].cat;

	size_t maxCat = this->_cfConfig.modelConfig.maxCategories;

	//CatCont::message("MaxCat=" + CatCont::toString(maxCat), "ParamExtractor::getParticipantCategory"); // DEBUG

	vector<double> catMu(maxCat);
	vector<double> catActive(maxCat);
	vector<double> catType(maxCat);

	for (size_t k = 0; k < maxCat; k++) {
		catMu[k] = pc.get(catInd.catMu[k]);
		catActive[k] = pc.get(catInd.catActive[k]);
		catType[k] = pc.get(catInd.catType[k]);
	}

	CategoryParam rval;
	rval.setFromVectors(catMu, catActive, catType);
	return rval;
}

SpecialParam ParamExtractor::getParticipantSpecial(const ParamContainer& pc, size_t partInd, const vector<size_t>& activeCatK) const {
	SpecialParam rval;
	const PI_Special& specialInd = this->part[partInd].special;

	rval.catSD_guessShift = pc.get(specialInd.catSD_guessShift);

	rval.pCatUsedToGuess.resize(activeCatK.size());
	for (size_t ak = 0; ak < activeCatK.size(); ak++) {
		rval.pCatUsedToGuess[ak] = pc.get(specialInd.pCatUsedToGuess[activeCatK[ak]]);
	}

	return rval;
}


/*
CompleteCategoryEffects ParamExtractor::getCategoryEffects(const ParamContainer& pc, const vector<size_t>& activeCatK) const {
	CompleteCategoryEffects rval;

	rval.paramPerAK.resize(activeCatK.size());
	for (size_t ak = 0; ak < activeCatK.size(); ak++) {

		size_t k = activeCatK[ak];

		const PI_CatEffect& picck = this->catEffects[k];
		SingleCategoryEffect& ccak = rval.paramPerAK[ak];

		ccak.pMem = pc.get(picck.pMem);
		ccak.pBetween = pc.get(picck.pBetween);
		ccak.pContBetween = pc.get(picck.pContBetween);
		ccak.pContWithin = pc.get(picck.pContWithin);

		ccak.contSD = pc.get(picck.contSD);
		ccak.catSD = pc.get(picck.catSD);
	}

	return rval;
}
*/

/*
// activeCatK is not required, it could be looked up based on partInd
CompleteCategoryEffects ParamExtractor::getCategoryTypeEffects(const ParamContainer& pc, size_t partInd, const vector<size_t>& activeCatK) const {
	// For this participant, get active categories.
	// For each active category, get its type.
	// Use type to look up matching category effects for that active category.

	// Alias this participant indices for the catTypes.
	const PI_Participant& pip = this->part.at(partInd);

	// Prepare the rval
	CompleteCategoryEffects rval;
	rval.paramPerAK.resize(activeCatK.size());

	for (size_t ak = 0; ak < activeCatK.size(); ak++) {

		// Use the index of this active category (k) to get the index of the category type for this participant for category k (catTypeInd).
		size_t k = activeCatK.at(ak);
		unsigned int catTypeInd = pc.get(pip.cat.catType.at(k));

		// DEBUG
		//CatCont::message("catType index: " + CatCont::toString(catTypeInd));

		// Alias the catType indices
		const PI_CatType& pict = this->catTypes.at(catTypeInd);

		// Alias the category effects for this active category
		SingleCategoryEffect& cek = rval.paramPerAK[ak];

		// Fill these category effects from these category type indices
		cek.pMem = pc.get(pict.pMem);
		cek.pBetween = pc.get(pict.pBetween);
		cek.pContBetween = pc.get(pict.pContBetween);
		cek.pContWithin = pc.get(pict.pContWithin);

		cek.contSD = pc.get(pict.contSD);
		cek.catSD = pc.get(pict.catSD);
	}

	return rval;
}
*/

CatTypeEffects ParamExtractor::getCategoryTypeEffects_v2(const ParamContainer& pc) const {

	size_t nTypes = this->catTypes.size();

	CatTypeEffects rval;
	rval.typeNames.resize(nTypes);
	rval.catTypeEffects.resize(nTypes);

	for (size_t i = 0; i < nTypes; i++) {

		// Alias the catType indices
		const PI_CatType& pict = this->catTypes.at(i);

		rval.typeNames[i] = pict.typeNameStr;

		// Reference the effect to fill it
		SingleCategoryEffect& sce = rval.catTypeEffects[i];

		// Fill these category effects from these category type indices
		sce.pMem = pc.get(pict.pMem);
		sce.pBetween = pc.get(pict.pBetween);
		sce.pContBetween = pc.get(pict.pContBetween);
		sce.pContWithin = pc.get(pict.pContWithin);

		sce.contSD = pc.get(pict.contSD);
		sce.catSD = pc.get(pict.catSD);
	}

	return rval;
}

ParamContainer::ParamProxy ParamExtractor::_makeProxy(const ParamContainer& pc, const std::string& name) {
	ParamContainer::ParamProxy proxy;
	if (!proxy.setup(pc, name)) {
		CatCont::stop("Error getting index for parameter " + proxy.name, "ParamExtractor::setup");
	}
	return proxy;
}

PI_Standard ParamExtractor::_setStandardProxies(const ParamContainer& pc, const string& indexStr) const {
	PI_Standard rval;

	rval.pMem = _makeProxy(pc, "pMem" + indexStr);
	rval.pBetween = _makeProxy(pc, "pBetween" + indexStr);
	rval.pContBetween = _makeProxy(pc, "pContBetween" + indexStr);
	rval.pContWithin = _makeProxy(pc, "pContWithin" + indexStr);
	rval.pCatGuess = _makeProxy(pc, "pCatGuess" + indexStr);

	rval.contSD = _makeProxy(pc, "contSD" + indexStr);
	rval.catSD = _makeProxy(pc, "catSD" + indexStr);
	rval.catSelectivity = _makeProxy(pc, "catSelectivity" + indexStr);

	return rval;
}

StandardParam ParamExtractor::_getStandard(const ParamContainer& pc, const PI_Standard& ps) const {
	StandardParam rval;

	rval.pMem = pc.get(ps.pMem);
	rval.pBetween = pc.get(ps.pBetween);
	rval.pContBetween = pc.get(ps.pContBetween);
	rval.pContWithin = pc.get(ps.pContWithin);
	rval.pCatGuess = pc.get(ps.pCatGuess);

	rval.contSD = pc.get(ps.contSD);
	rval.catSD = pc.get(ps.catSD);
	rval.catSelectivity = pc.get(ps.catSelectivity);

	// Interface idea:
	//rval.catSelectivity = pc[ps.catSelectivity]; // get only? or return double&?

	return rval;
}

////////////////////////////
// CategorizationParam_v2 //
////////////////////////////

/*
bool CategorizationParam_v2::setFromList(const ParameterList& param, string pnum, string cond, size_t maxCategories) {

	this->catSelectivity = extractParam_standard(param, "catSelectivity", pnum, cond, -1);

	this->activeCatMu.clear();

	string catIStrBase = ".part[" + pnum + ",";

	for (size_t k = 0; k < maxCategories; k++) {
		string catIstr = catIStrBase + CatCont::catIndexString(k) + "]";

		int catActive = param.at("catActive" + catIstr);

		if (catActive == 1) {
			this->activeCatMu.push_back(param.at("catMu" + catIstr));
		}
	}

}

bool CategorizationParam_v2::set(double catSelectivity_, vector<double> activeCatMu_) {
	this->catSelectivity = catSelectivity_;
	this->activeCatMu = activeCatMu_;
}
*/

} // namespace CatFirst
} // namespace CatCont