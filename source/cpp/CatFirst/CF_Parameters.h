#pragma once

//#include "Rcpp.h
//#include "./R_Interface.h"

#include "CCM_Main.h"
#include "CCM_Data.h"

#include "CCM_Util.h"
#include "GibbsSampler.h" // for ParameterList

#include "CF_ModelConfig.h"

namespace CatCont {
namespace CatFirst {

struct ParamExtractor;
struct SingleCategoryEffect;

double extractParam_standard(const ParameterList& param, string paramBaseName, string pnum, string cond, int cat);
double extractParam_category(const ParameterList& param, string paramBaseName, string pnum, int cat);

// A single parameter than can have participant and condition effects
struct SingleStandardParam {
	double value = 0;

	void setFromList(const ParameterList& param, string paramBaseName, string pnum, string cond, int cat);

	void toManifestProb(void);
	void toManifestSD(const SDRanges& sdRanges, DataType dataType);
};





struct ActiveCategoryParam {
	vector<size_t> k; // index of this active category
	vector<double> catMu;
	vector<size_t> catType;
};


// Categorization parameters can be shared by all participants or individual to each participant.
// They do not have cond or cat effects.
struct CategoryParam {

	// TODO: These are not all categorization param (catSelectivity). Rename to categoryParam?

	// 1) Track catMu and catActive
	vector<double> catMu;
	vector<unsigned int> catActive; // 0 or 1
	vector<size_t> catType; // index

	// 2) Track active categories and their indices
	//vector<double> activeCatMu;
	//vector<size_t> activeCatK; // Track index (k) of active cats
	// activeCatType? 
	//vector<unsigned int> getActiveCatType(void) const;

	// 3) Getters
	vector<double> getActiveCatMu(void) const; // This is used a lot
	vector<size_t> getActiveCatK(void) const;
	vector<size_t> getActiveCatType(void) const;
	ActiveCategoryParam getActiveParam(void) const;

	// General functions

	int nCatActive(int version = 1) const;

	void clear(void);

	//bool setFromList(const ParameterList& param, string pnum, size_t maxCategories);
	bool setFromVectors(vector<double> catMu_, vector<double> catActive_, vector<double> catType_);
};




// Standard parameters can vary by participant, condition, and/or category.
struct StandardParam {
	double pMem = 0;
	double pBetween = 0; // TODO: Keep or drop? Could change pBetween, pContBetween, and pContWithin into pCont. modelVariant selects likelihood.
	double pContBetween = 0;
	double pContWithin = 0;
	double pCatGuess = 0; // Can't vary with category.

	double contSD = 0;
	double catSD = 0;
	double catSelectivity = 0; // Can't vary with category.

	// indexStr can be like _pnum[123] or _cond[C2] or _cat[5]
	bool setFromList(const ParameterList& param, string indexStr);

	void transformToManifest(const SDRanges& sdRanges, DataType dataType);

	StandardParam operator+(const StandardParam& rhs) const;
	StandardParam operator+(const SingleCategoryEffect& rhs) const;

	map<string, double> toMap(void) const;
};

// Parameters that can vary with category
struct SingleCategoryEffect {
	double pMem = 0;
	double pBetween = 0;
	double pContBetween = 0;
	double pContWithin = 0;
	//double pCatGuess = 0;

	double contSD = 0;
	double catSD = 0;
	//double catSelectivity = 0;

	// This is not a shift. Maybe should not be in cat effect param, which are shifts. Moved to special.
	//double pCatUsedToGuess = 1;

	//void setFromList(const ParameterList& param, size_t k);
};



struct CompleteCategoryEffects {

	vector<SingleCategoryEffect> paramPerAK;
	
	//void setFromList(const ParameterList& param, const vector<size_t>& activeCatK);
};

// Filled with all types, not just active types.
// This is a candidate for being stored and only updated when needed. GibbsSampler feature like parameter update callback?
struct CatTypeEffects {
	vector<string> typeNames; // not really needed
	vector<SingleCategoryEffect> catTypeEffects;
};

// Strange things. Are these per participant?
struct SpecialParam {
	vector<double> pCatUsedToGuess; // One per k. TODO: Per ak?

	// When guessing, catSD is shifted (by the same amount for each category?)
	double catSD_guessShift = 0;

	//bool setFromList(const ParameterList& param, string pnum, size_t maxCategories);
	//void setFromList(const ParameterList& param, string pnum, const vector<size_t>& aks);
};




struct MPMP_memory {
	
	double pMem = -1;
	double pBetween = -1;
	double pContBetween = -1;
	double pContWithin = -1;

	double contSD = -1;
	double catSD = -1;
	//double catSelectivity = -1; // Moved to MPMP_categorization.

	void setFromLatent(const ModelConfiguration& modelConfig, StandardParam standard);
};

// Parameters needed in the guessing state
struct MPMP_guessing {

	double pCatGuess;

	// Per ak. This is because pCatUsedToGuess is not really a catType thing.
	vector<double> catSD;
	vector<double> pCatUsedToGuess; // From special

	void setFromLatent(const ModelConfiguration& modelConfig, double pCatGuess_, const vector<double>& catSD_, const vector<double>& pCatUsedToGuess_);

};

struct MPMP_categorization {
	CategoryParam cats;
	double catSelectivity = -1;

	void setFromLatent(const ModelConfiguration& modelConfig, const CategoryParam& catPar, double catSelectivity_);
};



struct MPMP_complete {
	MPMP_categorization categorization;

	// Per category type. Cornerstone type is index 0.
	vector<MPMP_memory> memory;

	MPMP_guessing guessing;

	//void setFromList(const ModelConfiguration& modelConfig, const ParameterList& param, string pnum, string cond);
	void setFromContainer(const ModelConfiguration& modelConfig, const ParamExtractor& pec, const ParamContainer& pc, size_t partInd, size_t condInd);

	void setFromLatent(const ModelConfiguration& modelConfig, const CategoryParam& partCat, const StandardParam& partStd, 
		const StandardParam& condStd, const CatTypeEffects& catEff, const SpecialParam& special);

	std::string printValues(void) const;
};






//////////////////////////////////////////////
// ParamExtractor stuff

struct PI_Standard {
	ParamContainer::ParamProxy pMem;
	ParamContainer::ParamProxy pBetween;
	ParamContainer::ParamProxy pContBetween;
	ParamContainer::ParamProxy pContWithin;
	ParamContainer::ParamProxy pCatGuess;

	ParamContainer::ParamProxy contSD;
	ParamContainer::ParamProxy catSD;
	ParamContainer::ParamProxy catSelectivity;
};

struct PI_Category {
	vector<ParamContainer::ParamProxy> catMu;
	vector<ParamContainer::ParamProxy> catActive;
	vector<ParamContainer::ParamProxy> catType;
};

struct PI_Special {
	ParamContainer::ParamProxy catSD_guessShift;
	//vector<ParamContainer::ParamProxy> catSD_guessShift; // technically could vary with participant*category
	vector<ParamContainer::ParamProxy> pCatUsedToGuess;
};

struct PI_Participant {
	string pnum;

	PI_Standard standard;

	PI_Category cat;
	//vector<PI_Category> cat; // without vectors in the struct. maybe

	PI_Special special;
};

struct PI_Condition {
	string condName;

	PI_Standard standard;
};

/*
struct PI_CatEffect {
	string catIndexStr;

	ParamContainer::ParamProxy pMem;
	ParamContainer::ParamProxy pBetween;
	ParamContainer::ParamProxy pContBetween;
	ParamContainer::ParamProxy pContWithin;

	ParamContainer::ParamProxy contSD;
	ParamContainer::ParamProxy catSD;
};
*/

// Each catType has a name and indexed parameter values
struct PI_CatType {
	string typeNameStr;

	ParamContainer::ParamProxy pMem;
	ParamContainer::ParamProxy pBetween;
	ParamContainer::ParamProxy pContBetween;
	ParamContainer::ParamProxy pContWithin;

	ParamContainer::ParamProxy contSD;
	ParamContainer::ParamProxy catSD;
};

// Stores indices of parameters in a ParamContainer
// ParamExtractor_Index
struct ParamExtractor {

	vector<PI_Participant> part;
	vector<PI_Condition> cond;
	vector<PI_CatType> catTypes;

	// Set the proxy indices based on pc
	void setup(const CF_CompositeConfig& cfConfig, const DataCollection& data, const ParamContainer& pc);

	// The confusing thing about the param extractor is that there are the index classes (PI_), the intermediate classes (StandardParam), then the MPMP.
	// It might be nice to go straight from PI_ to MPMP_ with a function like this:
	//MPMP_complete getCompleteMPMP(const ModelConfiguration& modelConfig, const ParamContainer& pc, size_t partInd, size_t condInd) const;
	MPMP_complete getCompleteMPMP(const ParamContainer& pc, size_t partInd, size_t condInd) const;

	// Strange interface
	StandardParam getParticipantStandard(const ParamContainer& pc, size_t partInd) const;
	CategoryParam getParticipantCategory(const ParamContainer& pc, size_t partInd) const;
	SpecialParam getParticipantSpecial(const ParamContainer& pc, size_t partInd, const vector<size_t>& activeCatK) const;

	StandardParam getConditionStandard(const ParamContainer& pc, size_t condInd) const;

	//CompleteCategoryEffects getCategoryEffects(const ParamContainer& pc, const vector<size_t>& activeCatK) const;
	//CompleteCategoryEffects getCategoryTypeEffects(const ParamContainer& pc, size_t partInd, const vector<size_t>& activeCatK) const;
	CatTypeEffects getCategoryTypeEffects_v2(const ParamContainer& pc) const;

private:

	CF_CompositeConfig _cfConfig;

	static ParamContainer::ParamProxy _makeProxy(const ParamContainer& pc, const std::string& name);

	PI_Standard _setStandardProxies(const ParamContainer& pc, const string& indexStr) const;
	StandardParam _getStandard(const ParamContainer& pc, const PI_Standard& ps) const;
};




////////////////////////////////////////////////////////////////////////////////////////
// Depreciated

/*
struct CategorizationParam_v2 {

	double catSelectivity;
	vector<double> activeCatMu;

	bool setFromList(const ParameterList& param, string pnum, string cond, size_t maxCategories);
	bool set(double catSelectivity_, vector<double> activeCatMu_);
};

struct MPMP_memory_v2 {
	double pMem;
	double pBetween;
	double pContBetween;
	double pContWithin;

	double catSD;
	double contSD;
};

struct MPMP_guessing_v2 {
	double pCatGuess; // one per part/cond

	// Per ak
	vector<double> catSD;
	vector<double> pCatUsedToGuess; // From special

	void setFromList(const ParameterList& param, string pnum, string cond, vector<size_t> activeCatInd) {
		SingleStandardParam temp;
		SingleStandardParam _catSD_mem;
		SingleStandardParam _catSD_guess;

		temp.setFromList(param, "pCatGuess", pnum, cond, -1);
		this->pCatGuess = temp.value;

		for (size_t ak : activeCatInd) {
			_catSD_mem.setFromList(param, "catSD", pnum, cond, ak);
			_catSD_guess.setFromList(param, "catSD_guessShift", pnum, "", -1); // no cond or cat.

			this->catSD.push_back(_catSD_mem.value + _catSD_guess.value);

			temp.setFromList(param, "pCatUsedToGuess", pnum, "", -1); // no cond or cat.
			this->pCatUsedToGuess.push_back(temp.value);
		}
	}

	void toManifest(const ModelConfiguration& modCfg) {

	}
};

struct MPMP_complete_v2 {
	CategorizationParam_v2 cat;

	vector<MPMP_memory_v2> memory; // When 0 categories active, memory has 1 element.
	MPMP_guessing_v2 guessing;

	void setFromList(const ModelConfiguration& modelConfig, const ParameterList& param, string pnum, string cond);

	// Convenience wrapper?
	void setFromRList(const ModelConfiguration& modelConfig, Rcpp::List paramList);

	void toManifest(const ModelConfiguration& modelConfig);
};
*/

struct CompleteParam {

	StandardParam standard;

	CategoryParam cat;

	SpecialParam special;
};

// An option so getParticipantParam returns both. or struct ParticipantParam = CompleteParam.
struct ParticipantParam {
	StandardParam standard;

	CategoryParam cat;

	SpecialParam special;
};

struct ParamPerK {
	vector<size_t> k; // Which k the standard are for (not all calculated)
	vector<StandardParam> standardPerK;

	// Not per K
	CategoryParam cat;

	// umm
	SpecialParam special;

	void set(bool activeOnly, const ModelConfiguration& modelConfig, const StandardParam& part, const StandardParam& cond, const CompleteCategoryEffects& catClass);
};

// Parameters used in the process model calculations for a single i,j,k combination
struct ManifestProcessModelParameters {

	// Track source?
	struct {
		size_t i;
		size_t j;
		size_t k;
	} indices;

	// Standard parameters that can vary with i, j, and k.
	double pMem = 0;
	double pBetween = 0;
	double pContBetween = 0;
	double pContWithin = 0;

	double contSD = 0;

	// catSD usually varies with i, sometimes j and k.
	double catSD_mem = 0; // For memory only
	double catSD_guess = 0; // For guessing only


	// Can't vary with k
	double pCatGuess = 0;
	double catSelectivity = 0; // Used to calculate weights. Varies by i and sometimes j.

	// Varies with k, but does not use study weights.
	double pCatUsedToGuess = 1;

	vector<double> activeCatMu; // Duplicated (can be same in all i, j, and k.)
	// size_t thisK; ?
	// double thisCatMu; ?
};


struct CatClassMap {

	// All vectors have the same length, K.

	// For each parameter, one element per k. Maps to category class parameters
	vector<string> pMem;

	// These 3 or just pCont
	vector<string> pBetween;
	vector<string> pContBetween;
	vector<string> pContWithin;
	//vector<string> pCont;

	vector<string> contSD;

	vector<string> catSD;

	// pCatGuess doesn't connect to the logic of this.
	// catSelectivity might connect to the logic, but is hard to estimate, so probably not worth doing.

	// R only
	/*
	bool setFromList(Rcpp::List ccList, size_t maxCategories) {
		if (ccList.containsElementNamed("pMem")) {
			pMem = ccList["pMem"];
		}
		else {
			pMem.resize(maxCategories, "constant=0"); // lol wtf
		}
	}
	*/
};

} // namespace CatFirst
} // namespace CatCont