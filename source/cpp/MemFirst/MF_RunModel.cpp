#include "MF_RunModel.h"

namespace CatCont {
namespace MemFirst {

// Decorrelating stuff didn't work for me. Maybe code is bad, maybe concept is bad, dunno.
bool usingDecorrelatingSteps = false; // Disabled, not possible to enable from R.

double MemFirstModelRunner::_llFunction(const CombinedParameters& par, const ConditionData& data) const {
	double rval = 0;
	if (this->config.dataType == DataType::Circular) {
		rval = Circular::betweenAndWithinLL(par, data, config);
	} else if (this->config.dataType == DataType::Linear) {
		rval = Linear::betweenAndWithinLL(par, data, config);
	}
	
	return rval;
}

double MemFirstModelRunner::singleParticipant_ll(const ParamContainer& param, unsigned int pIndex) const {

	ParticipantParameters partPar = _getParticipantParameters(param, pIndex);

	const ParticipantData& pdata = data.participants.at(pIndex);

	double llSum = 0;

	for (unsigned int condIndex = 0; condIndex < pdata.condData.size(); condIndex++) {

		ConditionParameters conditionPar = _getConditionParameters(param, condIndex);

		CombinedParameters combinedParameters = _combineParameters(partPar, conditionPar);

		llSum += _llFunction(combinedParameters, pdata.condData[condIndex]);

	}

	return llSum;
}

//double MemFirstModelRunner::_singleParticipantLLSamplingFunction(const ParamContainer& param, unsigned int pIndex) const {
//	return singleParticipant_ll(param, pIndex);
//}


double MemFirstModelRunner::singleCond_ll(const ParamContainer& param, unsigned int condIndex) const {

	ConditionParameters conditionPar = _getConditionParameters(param, condIndex);

	double llSum = 0;

	for (size_t pIndex = 0; pIndex < data.participants.size(); pIndex++) {
			
		ParticipantParameters partPar = _getParticipantParameters(param, pIndex);

		CombinedParameters combinedParameters = _combineParameters(partPar, conditionPar);

		const ConditionData& cd = data.participants.at(pIndex).condData.at(condIndex);
		llSum += _llFunction(combinedParameters, cd);
	}

	return llSum;
}

double MemFirstModelRunner::_allData_ll(const ParamContainer& param) const {

	double llSum = 0;

	for (size_t i = 0; i < this->data.participants.size(); i++) {
		llSum += this->singleParticipant_ll(param, i);
	}

	return llSum;
}

/*
double MemFirstModelRunner::genericConditionParameter_ll(double thisParam, const ParamContainer& param, unsigned int condIndex, string baseName) const {
	double ll = singleCond_ll(param, condIndex);
	double prior = cauchyLL(thisParam, priors.at(baseName + "_cond.loc"), priors.at(baseName + "_cond.scale"));
	return ll + prior;
}
*/

double MemFirstModelRunner::multiConditionParameter_ll(double thisParam, const ParamContainer& param, vector<unsigned int> condIndices, string baseName) const {

	double llSum = 0;
	for (unsigned int cond : condIndices) {
		llSum += singleCond_ll(param, cond);
	}

	//There is only 1 parameter, so only 1 prior.
	double prior = cauchyLL(thisParam, priors.at(baseName + "_cond.loc"), priors.at(baseName + "_cond.scale"));

	return llSum + prior;
}

double MemFirstModelRunner::normalPriorParameter_ll(double thisParam, const ParamContainer& param, unsigned int pIndex, string paramName) const {
	double ll = singleParticipant_ll(param, pIndex);
	double prior = normalLL(thisParam, param.at(paramName + "_part.mu"), param.at(paramName + "_part.var"));
	return ll + prior;
}



double MemFirstModelRunner::_distancePrior_catMu(const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const {
	string catIStrBase = "[" + data.participants.at(pIndex).pnum + ",";

	vector<double> mus(this->config.maxCategories);
	vector<unsigned int> catActives(this->config.maxCategories);

	for (unsigned int k = 0; k < this->config.maxCategories; k++) {
		string catIStr = catIStrBase + CatCont::catIndexString(k) + "]";

		mus[k] = param.at("catMu_part" + catIStr);
		catActives[k] = (unsigned int)param.at("catActive_part" + catIStr);
	}

	return _catPriorCalc.distancePrior_catMu(mus, catActives, catIndex, true); // log = true
}

double MemFirstModelRunner::_distancePrior_catActive(const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const {
	string catIStrBase = "[" + data.participants.at(pIndex).pnum + ",";

	vector<double> mus(this->config.maxCategories);
	vector<unsigned int> catActives(this->config.maxCategories);

	for (unsigned int k = 0; k < this->config.maxCategories; k++) {
		string catIStr = catIStrBase + CatCont::catIndexString(k) + "]";

		mus[k] = param.at("catMu_part" + catIStr);
		catActives[k] = (unsigned int)param.at("catActive_part" + catIStr);
	}

	return _catPriorCalc.distancePrior_catActive(mus, catActives, catIndex, true); // log = true
}


double MemFirstModelRunner::_catActive_ll(double catActive, const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const {

	bool thisCatActive = (int)catActive == 1;

	double ll = singleParticipant_ll(param, pIndex);

	double priorLL = 0;

	if (config.catActiveDistancePrior) {
		priorLL += _distancePrior_catActive(param, pIndex, catIndex);
	}

	// Used when catMu is shared, catActive is not shared, and there is a heirarchical prior on catActive.
	// This component of the prior on catActive relates to the other participants' catActives for this category.
	if (config.catActiveHierPrior) {

		string catIndexSuffix = "," + CatCont::catIndexString(catIndex) + "]";

		double catActiveSum = 0.0; // Sum across participants of catActive for this category.

		for (size_t i = 0; i < data.participants.size(); i++) {
			if (i == pIndex) {
				catActiveSum += 1; // For this participant, always treat the category as active to prevent the heir prior from going to 0.
			} else {
				string catActiveName = "catActive_part[" + data.participants.at(i).pnum + catIndexSuffix;
				catActiveSum += param.at(catActiveName);
			}
		}

		double pCatsActive = catActiveSum / data.participants.size();

		// Make prior weak, swings from 0.3 to 0.7. (0.4 is interval around 0.5 to swing)
		pCatsActive = 0.5 + (pCatsActive - 0.5) * 0.4;

		if (thisCatActive) {
			priorLL += std::log(pCatsActive);
		} else {
			priorLL += std::log(1 - pCatsActive);
		}

		//priorLL += std::log(thisCatActive ? pCatsActive : 1 - pCatsActive);

	}

	// Fixed prior prob (default 0.5 for active/inactive).
	if (this->priors.find("catActivePriorProb") != this->priors.end()) {
		double caPriorP = this->priors.at("catActivePriorProb");
		if (thisCatActive) {
			priorLL += std::log(caPriorP);
		} else {
			priorLL += std::log(1 - caPriorP);
		}
	}

	return ll + priorLL;
}

double MemFirstModelRunner::_catMu_ll(double catMu, const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const {

	double ll = 0;

	string catIStrBase = "[" + data.participants.at(pIndex).pnum + "," + CatCont::catIndexString(catIndex) + "]";
	string catActiveName = "catActive_part" + catIStrBase;
	if (param.at(catActiveName) == 1.0) {
		ll = singleParticipant_ll(param, pIndex); //This only needs to be calculated for active categories.
	}

	//double prior = _catMuPenalityPrior(param, pIndex, catIndex);
	//double prior = _catPriorCalc.distancePrior_catMu(param, data.participants.at(pIndex).pnum, catIndex);
	double prior = _distancePrior_catMu(param, pIndex, catIndex);

	return ll + prior;
}

// catActive can only be shared if catMu is also shared.
double MemFirstModelRunner::_catActive_shared_ll(double catActive, const ParamContainer& param, unsigned int catIndex) const {

	double dataLL = 0;
	double priorLL = 0;

	for (size_t i = 0; i < this->data.participants.size(); i++) {
		dataLL += singleParticipant_ll(param, i);
	}

	if (config.catActiveDistancePrior) {
		//priorLL += _catActivePenaltyPrior(param, 0, catIndex); 
		priorLL += _distancePrior_catActive(param, 0, catIndex); // pInd=0 to just calculate for shared holder.
	}

	if (this->priors.find("catActivePriorProb") != this->priors.end()) {
		double caPriorP = this->priors.at("catActivePriorProb");
		if ((int)catActive == 1) {
			priorLL += std::log(caPriorP);
		} else {
			priorLL += std::log(1 - caPriorP);
		}
	}

	return dataLL + priorLL;
}

double MemFirstModelRunner::_catMu_shared_ll(double catMu, const ParamContainer& param, unsigned int catIndex) const {

	// likelihood for all participants and conditions
	double ll = this->_allData_ll(param);

	double priorLL = 0;
	if (this->_catActiveShared) {
		// If catActive is shared, then catMu and catActive are the same for all participants.
		// Only calculate the prior once, for the holder participant (index 0).
		priorLL += _distancePrior_catMu(param, 0, catIndex);
	}
	else {
		// If catActive is not shared, the prior is different for each participant.
		// TODO: Some thought is needed. Average likelihood? Use the union of catActive from all participants?
		for (size_t i = 0; i < this->data.participants.size(); i++) {
			//priorLL += _catMuPenalityPrior(param, i, catIndex);

			priorLL += _distancePrior_catMu(param, i, catIndex);
		}

	}

	//double priorLL = _catPriorCalc.distancePrior_catMu(param, 0, catIndex); // pInd=0 to just calculate for shared holder.

	return ll + priorLL;
}


double MemFirstModelRunner::postMuSample(const ParamContainer& param, string baseParamName, double mu0, double var0) const {
	string param_part = baseParamName + "_part";
	vector<double> y = this->gibbs.getCurrentGroupValues(param_part);
	return normal_muPostSample(y, param.at(param_part + ".var"), mu0, var0);
}

double MemFirstModelRunner::postVarSample(const ParamContainer& param, string baseParamName, double a0, double b0) const {
	string param_part = baseParamName + "_part";
	vector<double> y = this->gibbs.getCurrentGroupValues(param_part);
	return normal_varPostSample(y, param.at(param_part + ".mu"), a0, b0);
}


void MemFirstModelRunner::createParameters(void) {

	using namespace std::placeholders;

	this->_setMhTuning();
	this->_setPriors();

	gibbs.clear();

	gibbs.sectionTracker.sectionStart("Total");

	//////////////////////////////////
	// Primary participant parameters
	gibbs.sectionTracker.sectionStart("Participant Primary");
	for (size_t i = 0; i < data.participants.size(); i++) {

		vector<string> probParam = { "pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess" };
		for (string pp : probParam) {

			string pp_part = pp + "_part";

			MH_Parameter par;
			par.name = pp_part + "[" + data.participants[i].pnum + "]";
			par.group = pp_part;

			par.deviateFunction = bind(normalDeviate, _1, mhTuningSd.at(pp_part));
			par.llFunction = bind(&MemFirstModelRunner::normalPriorParameter_ll, this, _1, _2, i, pp);

			gibbs.addParameter(par, uniformDeviate(-2, 2));
		}
		
		vector<string> sdParam = { "contSD", "catSelectivity", "catSD" };
		for (string sdp : sdParam) {

			string sdp_part = sdp + "_part";

			MH_Parameter par;
			par.name = sdp_part + "[" + data.participants[i].pnum + "]";
			par.group = sdp_part;

			par.deviateFunction = bind(normalDeviate, _1, mhTuningSd.at(sdp_part));
			par.llFunction = bind(&MemFirstModelRunner::normalPriorParameter_ll, this, _1, _2, i, sdp);

			gibbs.addParameter(par, uniformDeviate(config.sdRanges.minSd + 2, 40));
		}
	}
	gibbs.sectionTracker.sectionEnd("Participant Primary");

	//////////////////////////////////
	// Category parameters
	gibbs.sectionTracker.sectionStart("Category Parameters");

	// The start value for catMu are a grid within the response range.
	// Note that this uses the data response range rather than the user-settable config response range.
	// This works in the same way for both linear and circular. Imagine a circular design with data on only part of the circle.
	std::vector<double> catMuStartGrid(config.maxCategories);
	double catMuGridStepSize = (data.responseRange.upper - data.responseRange.lower) / (double)config.maxCategories;
	for (unsigned int k = 0; k < config.maxCategories; k++) {
		catMuStartGrid[k] = (k + 0.5) * catMuGridStepSize + data.responseRange.lower;
	}
	// Start grid uses data range, but data are in radians and catMu param are in degrees.
	if (config.dataType == DataType::Circular) {
		catMuStartGrid = CatCont::Circular::radiansToDegrees(catMuStartGrid);
	}

	// If a category parameter is shared between participants, 
	// use the first participant as the holder for the shared parameters.
	// Copy the shared values to all other participants with DependentParameter.
	const unsigned int sharedHolderIndex = 0;
	bool catMuShared = config.sharedParameters.count("catMu") > 0;
	this->_catActiveShared = config.sharedParameters.count("catActive") > 0;
	const bool shouldModuloCatMu = config.dataType == DataType::Circular;

	for (size_t i = 0; i < data.participants.size(); i++) {
		for (unsigned int k = 0; k < config.maxCategories; k++) {

			string catIstr = "[" + data.participants.at(i).pnum + "," + CatCont::catIndexString(k) + "]";

			// ---------------------
			// catMu
			if (catMuShared) {

				if (i == sharedHolderIndex) {
					MH_Parameter sharedCatMu;
					sharedCatMu.name = "catMu_part" + catIstr;
					sharedCatMu.group = "catMu_part";

					//double catMuCandidateSd = mhTuningSd.at("catMu_part"); //linear
					//if (config.dataType == DataType::Circular) {
					//	catMuCandidateSd = CatCont::Circular::degreesToRadians(catMuCandidateSd);
					//}
					//sharedCatMu.deviateFunction = bind(normalDeviate, _1, catMuCandidateSd);

					// catMu estimated in degrees or linear units
					sharedCatMu.deviateFunction = bind(CatCont::MemFirst::catMuDeviateFunction, _1, mhTuningSd.at("catMu_part"), shouldModuloCatMu);

					sharedCatMu.llFunction = bind(&MemFirstModelRunner::_catMu_shared_ll, this, _1, _2, k);

					gibbs.addParameter(sharedCatMu, catMuStartGrid[k]);
				} else {
					// If not the holder, catMu and catActive are dependent parameters
					string catMuSource = "catMu_part[" + data.participants.at(sharedHolderIndex).pnum + "," + CatCont::catIndexString(k) + "]";
					string catMuName = "catMu_part" + catIstr;

					DependentParameter catMuDep(catMuName, "catMu_part", catMuSource);
					gibbs.addParameter(catMuDep);
				}

			} else {
				// catMu individual
				MH_Parameter catMu;
				catMu.name = "catMu_part" + catIstr;
				catMu.group = "catMu_part";

				//double catMuCandidateSd = mhTuningSd.at("catMu_part"); //linear
				//if (config.dataType == DataType::Circular) {
				//	catMuCandidateSd = CatCont::Circular::degreesToRadians(catMuCandidateSd);
				//}
				//catMu.deviateFunction = bind(normalDeviate, _1, catMuCandidateSd);

				// catMu estimated in degrees or linear units
				catMu.deviateFunction = bind(CatCont::MemFirst::catMuDeviateFunction, _1, mhTuningSd.at("catMu_part"), shouldModuloCatMu);

				catMu.llFunction = bind(&MemFirstModelRunner::_catMu_ll, this, _1, _2, i, k);

				gibbs.addParameter(catMu, catMuStartGrid[k]);
			}

			// ---------------------
			// catActive
			if (this->_catActiveShared) {
				if (i == sharedHolderIndex) {
					MH_Parameter sharedCatActive;
					sharedCatActive.name = "catActive_part" + catIstr;
					sharedCatActive.group = "catActive_part";

					sharedCatActive.deviateFunction = CatCont::MemFirst::catActiveDeviateFunction;
					sharedCatActive.llFunction = bind(&MemFirstModelRunner::_catActive_shared_ll, this, _1, _2, k);

					gibbs.addParameter(sharedCatActive, 0.0); // Start all categories inactive. TODO
				} else {
					string catActiveSource = "catActive_part[" + data.participants.at(sharedHolderIndex).pnum + "," + CatCont::catIndexString(k) + "]";
					string catActiveName = "catActive_part" + catIstr;

					DependentParameter catActiveDep(catActiveName, "catActive_part", catActiveSource);
					gibbs.addParameter(catActiveDep);
				}
			} else {
				// catActive individual
				MH_Parameter catActive;
				catActive.name = "catActive_part" + catIstr;
				catActive.group = "catActive_part";

				catActive.deviateFunction = CatCont::MemFirst::catActiveDeviateFunction;
				catActive.llFunction = bind(&MemFirstModelRunner::_catActive_ll, this, _1, _2, i, k);

				gibbs.addParameter(catActive, 0.0); //Start all categories inactive. TODO
			}

		}
	}

	gibbs.sectionTracker.sectionEnd("Category Parameters");


	/////////////////////////////////////////////
	// hierarchical priors on participant parameters
	gibbs.sectionTracker.sectionStart("Hierarchical Priors on Participant", nullptr);

	vector<string> paramWithHierarchicalPrior = config.getParamWithHierachicalPriors();

	for (const string& php : paramWithHierarchicalPrior) {

		ConjugateParameter popMu;

		popMu.name = php + "_part.mu";
		popMu.group = "populationParam";

		popMu.samplingFunction = std::bind(&MemFirstModelRunner::postMuSample, this, _1, php, priors.at(php + "_part.mu.mu"), priors.at(php + "_part.mu.var"));

		gibbs.addParameter(popMu, 0);


		ConjugateParameter popVar;

		popVar.name = php + "_part.var";
		popVar.group = "populationParam";

		popVar.samplingFunction = std::bind(&MemFirstModelRunner::postVarSample, this, _1, php, priors.at(php + "_part.var.a"), priors.at(php + "_part.var.b"));

		gibbs.addParameter(popVar, 30);

	}
	gibbs.sectionTracker.sectionEnd("Hierarchical Priors on Participant");
	// end hierarchical priors
	/////////////////////////////////////////////


	/////////////////////////////////////////////
	// Condition effect parameters
	gibbs.sectionTracker.sectionStart("Condition Effects");
	
	EqualityConstraints eq;
	
	bool eqSetupSuccess = eq.setup(config.overrides.equalityConstraints, data.conditionNames, vector<string>(0));
	if (!eqSetupSuccess) {
		CatCont::message("Error while setting up equality constraints for condition effects.", "createParameters");
		return;
	}

	vector<string> conditionEffectsToCreate = config.getParamWithAndWithoutConditionEffects();
		
	for (size_t condIndex = 0; condIndex < data.conditionNames.size(); condIndex++) {

		string cstr = "[" + data.conditionNames[condIndex] + "]";

		for (const string& parName : conditionEffectsToCreate) {

			string name = parName + "_cond" + cstr;
			string group = parName + "_cond";

			if (data.conditionNames[condIndex] == config.cornerstoneConditionName) {
			//if (condIndex == config.cornerstoneConditionIndex) {

				gibbs.createConstantParameter(name, group, 0);

			} else {

				string source = eq.getSourceParameter(name);

				if (source == EqualityConstraints::FreeParameter) {
					//set up free parameter with multiConditionParameter_ll

					//unsigned int sourceIndex = condIndex;
					vector<unsigned int> conditionIndices = eq.getEqualConditionIndices(parName, condIndex);

					MH_Parameter par;
					par.name = name;
					par.group = group;

					par.deviateFunction = std::bind(normalDeviate, _1, mhTuningSd.at(par.group));
					par.llFunction = std::bind(&MemFirstModelRunner::multiConditionParameter_ll, this, _1, _2, conditionIndices, parName);

					gibbs.addParameter(par, 0);

				} else {

					//Not a free parameter: Set up dependent parameter with a source.
					DependentParameter dp;
					dp.name = name;
					dp.group = group;

					dp.sourceParameter = source;

					gibbs.addParameter(dp);
				}

			}
		} // parName

	} //condIndex

	gibbs.sectionTracker.sectionEnd("Condition Effects");

	// end condition parameters
	/////////////////////////////////////////////


	/////////////////////////////////////////////
	// decorrelating steps
	
	if (usingDecorrelatingSteps) {
		gibbs.sectionTracker.sectionStart("Decorrelating Steps");

		vector<string> decorrelatedParameters = config.getParamWithAndWithoutConditionEffects();

		for (const string& s : decorrelatedParameters) {

			DecorrelatingStep dcs;

			dcs.name = s + "_deCor";
			dcs.group = "deCor";

			dcs.currentValuesFunction = bind(&MemFirstModelRunner::decorrelatingCurrent, this, _1, s);
			dcs.candidateValuesFunction = bind(&MemFirstModelRunner::decorrelatingCandidate, this, _1, s, mhTuningSd.at(dcs.name));

			dcs.llFunction = bind(&MemFirstModelRunner::decorrelating_ll, this, _1, _2, s);

			gibbs.addParameter(dcs, 0);

		}

		gibbs.sectionTracker.sectionEnd("Decorrelating Steps");
	}

	// end decorrelating steps
	/////////////////////////////////////////////


	///////////////////
	// Then go back and zero out conditions (and decorrelating steps) that don't have condition parameters.
	// This seems stupid, but it allows for easier extensibility.
	vector<string> condParamToZero = config.getParamWithoutConditionEffects();

	for (unsigned int i = 0; i < condParamToZero.size(); i++) {
		//TODO: This is not currently needed, because of how the condition effects make equality constraints.
		//The equality constraints make parameters without condition effects all be equal to the cornerstone, which is 0.
		gibbs.setParameterGroupToConstantValue(condParamToZero[i] + "_cond", 0);

		if (usingDecorrelatingSteps) {
			gibbs.replaceParameter(condParamToZero[i] + "_deCor", ConstantParameter(0, condParamToZero[i] + "_deCor", "deCor"), 0);
		}
	}


	//between-item variant: set pBetween and pContWithin to 1 (and also their condition effects to 0)
	if (config.modelVariant == ModelVariant::BetweenItem) {

		vector<string> pToSet = { "pBetween", "pContWithin" };

		for (string& s : pToSet) {
			gibbs.setParameterGroupToConstantValue(s + "_part", 100); //100 being essentially 1 once transformed, 
															//but it actually doesn't matter because the between likelihood function is used.
			gibbs.setParameterGroupToConstantValue(s + "_cond", 0);

			//Get rid of population parameters
			gibbs.replaceParameter(s + "_part.mu", ConstantParameter(0, s + "_part.mu", "populationParam"), 0);
			gibbs.replaceParameter(s + "_part.var", ConstantParameter(1, s + "_part.var", "populationParam"), 1);

			if (usingDecorrelatingSteps) {
				gibbs.replaceParameter(s + "_deCor", ConstantParameter(0, s + "_deCor", "deCor"), 0);
			}
		}
	}

	//within-item variant: set pBetween and pContWithin to 0 (and also their condition effects to 0)
	if (config.modelVariant == ModelVariant::WithinItem) {
		vector<string> pToSet = { "pBetween", "pContBetween" };

		for (string& s : pToSet) {
			gibbs.setParameterGroupToConstantValue(s + "_part", -100); //-100 being essentially 0 once transformed, 
															 //but it actually doesn't matter because the within likelihood function is used.
			gibbs.setParameterGroupToConstantValue(s + "_cond", 0);

			//Get rid of population parameters
			gibbs.replaceParameter(s + "_part.mu", ConstantParameter(0, s + "_part.mu", "populationParam"), 0);
			gibbs.replaceParameter(s + "_part.var", ConstantParameter(1, s + "_part.var", "populationParam"), 1);

			if (usingDecorrelatingSteps) {
				gibbs.replaceParameter(s + "_deCor", ConstantParameter(0, s + "_deCor", "deCor"), 0);
			}
		}
	}


	//ZL model variant: pBetween = 1, pContBetween = 1, pContWithin = 1, catActive = 0, maxCategories = 0
	if (config.modelVariant == ModelVariant::ZL) {
			
		//100 is tranformed to 1, -100 is transformed to 0
		map<string, double> paramValues;
		paramValues["pBetween"] = 100; //whatever: doesn't have to be 1
		paramValues["pContWithin"] = 100;
		paramValues["pContBetween"] = 100;

		paramValues["pCatGuess"] = -100;

		paramValues["catSelectivity"] = 2; //whatever: doesn't have to be 2
		paramValues["catSD"] = 2;

		for (auto it = paramValues.begin(); it != paramValues.end(); it++) {

			gibbs.setParameterGroupToConstantValue(it->first + "_part", it->second);

			gibbs.setParameterGroupToConstantValue(it->first + "_cond", 0);

			//Get rid of population parameters
			gibbs.replaceParameter(it->first + "_part.mu", ConstantParameter(0, it->first + "_part.mu", "populationParam"), 0);
			gibbs.replaceParameter(it->first + "_part.var", ConstantParameter(1, it->first + "_part.var", "populationParam"), 1);

			if (usingDecorrelatingSteps) {
				gibbs.replaceParameter(it->first + "_deCor", ConstantParameter(0, it->first + "_deCor", "deCor"), 0);
			}

		}

		//Also set to 0 for all catActive
		gibbs.setParameterGroupToConstantValue("catActive_part", 0); // this seems unecessary
	}


	if (config.calculateParticipantLikelihoods) {
		gibbs.sectionTracker.sectionStart("Participant Likelihood");

		for (size_t i = 0; i < data.participants.size(); i++) {

			string istr = "[" + data.participants.at(i).pnum + "]";

			CalculatedParameter spLL;
			spLL.name = "participantLL" + istr;
			spLL.group = "participantLL";

			spLL.updateContinuously = false; // Only update once per iteration

			//spLL.samplingFunction = std::bind(&MemFirstModelRunner::_singleParticipantLLSamplingFunction, this, _1, i);
			spLL.samplingFunction = std::bind(&MemFirstModelRunner::singleParticipant_ll, this, _1, i);

			gibbs.addParameter(spLL);
		}

		gibbs.sectionTracker.sectionEnd("Participant Likelihood");
	}

	gibbs.sectionTracker.sectionEnd("Total");
	
	if (!this->runConfig.profileParameterTypes) {
		gibbs.sectionTracker.clear();
	}

	this->_doStartingValueOverrides();

	this->_doConstantParameterOverrides();
}






// Gets the latent parameters associated with a single participant. NO TRANSFORMATIONS ARE APPLIED HERE.
ParticipantParameters MemFirstModelRunner::_getParticipantParameters(const ParamContainer& param, unsigned int pIndex) const {
	return CatCont::MemFirst::getParticipantParameters(param, this->data.participants.at(pIndex).pnum, this->config.maxCategories);
}

// Gets the latent parameters associated with a single condition
ConditionParameters MemFirstModelRunner::_getConditionParameters(const ParamContainer& param, unsigned int condIndex) const {
	return CatCont::MemFirst::getConditionParameters(param, this->data.conditionNames[condIndex]);
}

// Combines the latent participant and condition parameters, transforming them to the mainfest space
CombinedParameters MemFirstModelRunner::_combineParameters(const ParticipantParameters& part, const ConditionParameters& cond) const {
	return CatCont::MemFirst::combineParameters(part, cond, this->config.sdRanges, this->config.dataType);
}


/*
void MemFirstModelRunner::setData(vector<ParticipantData> partData) {

	data = CatCont::DataCollection(); //reset the data struct

	data.participants = partData;

	//Get the condition names and cornerstone condition
	//config.cornerstoneConditionIndex = -1;

	for (size_t j = 0; j < data.participants.front().condData.size(); j++) {
		string conditionName = data.participants.front().condData[j].condition;
		data.conditionNames.push_back(conditionName);

		//if (conditionName == this->config.cornerstoneConditionName) {
		//	config.cornerstoneConditionIndex = j;
		//}
	}

	//Get the data ranges
	for (size_t i = 0; i < data.participants.size(); i++) {
		for (size_t j = 0; j < data.conditionNames.size(); j++) {

			const ConditionData& cd = data.participants[i].condData[j];

			for (size_t obs = 0; obs < cd.study.size(); obs++) {

				data.studyRange.lower = min(data.studyRange.lower, cd.study[obs]);
				data.studyRange.upper = max(data.studyRange.upper, cd.study[obs]);

				data.responseRange.lower = min(data.responseRange.lower, cd.response[obs]);
				data.responseRange.upper = max(data.responseRange.upper, cd.response[obs]);
			} // obs
		} // j
	} // i

}
*/



map<string, double> MemFirstModelRunner::getDefaultPriors(void) {

	map<string, double> defPriors;

	vector<string> probParams = { "pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess" };

	for (const string& pp : probParams) {

		// Hierarchical priors on participant parameters
		defPriors[pp + "_part.mu.mu"] = 0;
		defPriors[pp + "_part.mu.var"] = pow(1.2, 2); //sd = 1.2, var = 1.44
		defPriors[pp + "_part.var.a"] = 1.5;
		defPriors[pp + "_part.var.b"] = 1.5;

		// Priors on condition effects
		defPriors[pp + "_cond.loc"] = 0;
		defPriors[pp + "_cond.scale"] = 0.3;

	}

	// TODO: Apply linear range to SD params and catMu
	//Argument: std::vector<double> linearRange = { 0, 100 };
	//double sdWidth = 360;
	//if (linearRange.size() == 2) {
	//	sdWidth = abs(linearRange[1] - linearRange[2]);
	//}
	//double sdScale = 1;
	//if (linearRange.size() == 2) {
	//	sdScale = abs(linearRange[1] - linearRange[2]) / 360.0;
	//}

	vector<string> sdParams = { "contSD", "catSD", "catSelectivity" };

	for (const string& sp : sdParams) {

		// Hierarchical priors on participant parameters
		defPriors[sp + "_part.mu.mu"] = 35;
		defPriors[sp + "_part.mu.var"] = pow(15, 2); //sd = 15
		defPriors[sp + "_part.var.a"] = 0.75;
		defPriors[sp + "_part.var.b"] = 0.75;

		// Priors on condition effects
		defPriors[sp + "_cond.loc"] = 0;
		defPriors[sp + "_cond.scale"] = 3;

	}

	defPriors["catMuPriorSD"] = 12;
	defPriors["catActivePriorProb"] = 0.5;

	return defPriors;
}

map<string, double> MemFirstModelRunner::getDefaultMHTuning(void) {

	// TODO: Apply linear range to SD params and catMu for both participant and condition params.

	// These are standard deviations of normal candidate distributions.
	map<string, double> defTuningSD;

	///////////////////////////////////////////
	// Participant parameters. All are latent.
	defTuningSD["pMem_part"] = 0.3;
	defTuningSD["pBetween_part"] = 0.9;
	defTuningSD["pContBetween_part"] = 0.4;
	defTuningSD["pContWithin_part"] = 0.7;
	defTuningSD["pCatGuess_part"] = 0.7;

	// SD parameters and catMu are in degrees
	defTuningSD["contSD_part"] = 2;
	defTuningSD["catSD_part"] = 1;
	defTuningSD["catSelectivity_part"] = 2;

	defTuningSD["catMu_part"] = 5;

	///////////////////////////////////////////
	// Condition effects. All are latent.
	defTuningSD["pMem_cond"] = 0.2;
	defTuningSD["pBetween_cond"] = 0.4;
	defTuningSD["pContBetween_cond"] = 0.2;
	defTuningSD["pContWithin_cond"] = 0.2;
	defTuningSD["pCatGuess_cond"] = 0.3;

	defTuningSD["contSD_cond"] = 1;
	defTuningSD["catSD_cond"] = 0.5;
	defTuningSD["catSelectivity_cond"] = 0.8;


	if (usingDecorrelatingSteps) {

		// When using decorrelating steps, the normal parameters need smaller step sizes by about a half.
		for (auto& sd : defTuningSD) {
			sd.second /= 2;
		}
		// Unadjust params without decorrelating steps.
		defTuningSD["catMu_part"] *= 2;


		//decorrelating parameters
		defTuningSD["pMem_deCor"] = 0.1;
		defTuningSD["pBetween_deCor"] = 0.2;
		defTuningSD["pContBetween_deCor"] = 0.2;
		defTuningSD["pContWithin_deCor"] = 0.3;
		defTuningSD["pCatGuess_deCor"] = 0.1;

		defTuningSD["contSD_deCor"] = 2;
		defTuningSD["catSD_deCor"] = 1;
		defTuningSD["catSelectivity_deCor"] = 1;
	}

	return defTuningSD;
}


void MemFirstModelRunner::_setPriors(void) {

	this->priors = getDefaultPriors();

	//Once all of the defaults are in, do the overrides
	_doPriorOverrides();

	//After the overrides are in, calculate these things:
	_catPriorCalc.setup(this->config, this->priors.at("catMuPriorSD"));
}

void MemFirstModelRunner::_setMhTuning(void) {
	this->mhTuningSd = getDefaultMHTuning();

	// Once defaults have been set, override defaults
	_doMhOverrides();
}



/*
ParticipantParameters MemFirstModelRunner::getParticipantParameters(const ParamContainer& param, string pnum, unsigned int maxCategories) {
	string istr = "[" + pnum + "]";

	ParticipantParameters part;
	part.pMem = param.at("pMem" + istr);
	part.pBetween = param.at("pBetween" + istr);
	part.pContBetween = param.at("pContBetween" + istr);
	part.pContWithin = param.at("pContWithin" + istr);
	part.pCatGuess = param.at("pCatGuess" + istr);

	part.contSD = param.at("contSD" + istr);

	string catIStrBase = "[" + pnum + ",";

	for (unsigned int k = 0; k < maxCategories; k++) {
		string catIstr = catIStrBase + CatCont::catIndexString(k) + "]";

		double active = param.at("catActive" + catIstr);

		if (active == 1.0) {
			part.cat.mu.push_back(param.at("catMu" + catIstr));
		}
		//if not active, don't add it to the list
	}

	part.cat.selectivity = param.at("catSelectivity" + istr);
	part.cat.SD = param.at("catSD" + istr);

	return part;
}

ConditionParameters MemFirstModelRunner::getConditionParameters(const ParamContainer& param, string condName) {

	ConditionParameters cond;

	string istr = "[" + condName + "]";

	cond.pMem = param.at("pMem_cond" + istr);
	cond.pBetween = param.at("pBetween_cond" + istr);
	cond.pContBetween = param.at("pContBetween_cond" + istr);
	cond.pContWithin = param.at("pContWithin_cond" + istr);
	cond.pCatGuess = param.at("pCatGuess_cond" + istr);

	cond.contSD = param.at("contSD_cond" + istr);
	cond.cat.selectivity = param.at("catSelectivity_cond" + istr);
	cond.cat.SD = param.at("catSD_cond" + istr);

	return cond;
}

CombinedParameters MemFirstModelRunner::combineParameters(const ParticipantParameters& part, const ConditionParameters& cond, const SDRanges& sdRanges, DataType dataType) {
	CombinedParameters rval;

	rval.pMem = CatCont::paramTransform_probability(part.pMem + cond.pMem);
	rval.pBetween = CatCont::paramTransform_probability(part.pBetween + cond.pBetween);
	rval.pContBetween = CatCont::paramTransform_probability(part.pContBetween + cond.pContBetween);
	rval.pContWithin = CatCont::paramTransform_probability(part.pContWithin + cond.pContWithin);
	rval.pCatGuess = CatCont::paramTransform_probability(part.pCatGuess + cond.pCatGuess);

	rval.contSD = CatCont::paramTransform_sd(part.contSD + cond.contSD, sdRanges, dataType);
	rval.cat.selectivity = CatCont::paramTransform_sd(part.cat.selectivity + cond.cat.selectivity, sdRanges, dataType);
	rval.cat.SD = CatCont::paramTransform_sd(part.cat.SD + cond.cat.SD, sdRanges, dataType);

	// catMu are estimated in radians. No condition parameters or transformations are applied.
	// The only possible transformation would be fmod(mu, 2*PI), but it is unnecessary, so skipped to increase estimation speed.
	rval.cat.mu = part.cat.mu;

	return rval;
}





//This is a combination of the g and h functions from HVR17
double MemFirstModelRunner::_catMuPenaltyDensity(const vector<double>& mus, const vector<unsigned int>& catActives, unsigned int catIndex) const {

	double density = 1;

	for (unsigned int k = 0; k < mus.size(); k++) {

		if (k != catIndex) {
			if (catActives[catIndex] == 1 && catActives[k] == 1) {

				double like;
				if (config.dataType == DataType::Circular) {
					like = vmLut.dVonMises(mus[catIndex], mus[k], this->_catMuPriorData.kappa);
				}
				else {
					like = Linear::normalPDF(mus[catIndex], mus[k], _catMuPriorData.sd); //NOT truncated
				}

				double ratio = 1 - like / _catMuPriorData.maxLikelihood;

				density *= (ratio * ratio); //square the ratio to make flatter bottoms on the notches
			}
		}
	}

	return density; //return non-log density
}

// This is an approximation of the S function from HVR17
// mus are copied not referenced because they are modified
double MemFirstModelRunner::_estimateCatMuScaleFactor(vector<double> mus, const vector<unsigned int>& catActives, unsigned int catIndex, unsigned int steps) const {
	//circular
	double stepSize = 2 * PI / steps;
	double startPoint = 0;

	if (config.dataType == DataType::Linear) {
		stepSize = (config.linearConfiguration.catMu.upper - config.linearConfiguration.catMu.lower) / steps;
		startPoint = config.linearConfiguration.catMu.lower;
	}

	double densSum = 0;

	for (unsigned int i = 0; i < steps; i++) {
		mus[catIndex] = startPoint + i * stepSize;

		densSum += _catMuPenaltyDensity(mus, catActives, catIndex);
	}

	return 1.0 / (densSum * stepSize);
}

//This is the f function from HVR17
double MemFirstModelRunner::_scaledCatMuDensity(const vector<double>& mus, const vector<unsigned int>& catActives, unsigned int catIndex, unsigned int steps) const {

	if (config.dataType == DataType::Linear) {
		// The range of catMu defines a uniform distribution that is multiplied by the rest of the prior.
		// If mu[catIndex] is out of range, it has 0 density.
		if (mus[catIndex] < config.linearConfiguration.catMu.lower || mus[catIndex] > config.linearConfiguration.catMu.upper) {
			return 0;
		}
	}

	//If this category is inactive, the density is just the height of a uniform distribution
	if (catActives[catIndex] == 0) {
		if (config.dataType == DataType::Linear) {
			return 1.0 / (config.linearConfiguration.catMu.upper - config.linearConfiguration.catMu.lower);
		}
		else {
			//circular
			return 1.0 / (2 * PI);
		}
	}

	double unscaled = _catMuPenaltyDensity(mus, catActives, catIndex);
	double scaleFactor = _estimateCatMuScaleFactor(mus, catActives, catIndex, steps);

	return unscaled * scaleFactor;
}

double MemFirstModelRunner::_catMuPenalityPrior(const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const {
	string catIStrBase = "[" + data.participants.at(pIndex).pnum + ",";

	vector<double> mus(config.maxCategories);
	vector<unsigned int> catActives(config.maxCategories);

	for (unsigned int k = 0; k < config.maxCategories; k++) {
		string catIStr = catIStrBase + CatCont::catIndexString(k) + "]";

		mus[k] = param.at("catMu" + catIStr);
		catActives[k] = (unsigned int)param.at("catActive" + catIStr);
	}


	double density = _scaledCatMuDensity(mus, catActives, catIndex, config.catMuPriorApproximationPrecision);

	return std::log(density);
}

double MemFirstModelRunner::_catActivePenaltyPrior(const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const {
	string catIStrBase = "[" + data.participants.at(pIndex).pnum + ",";

	vector<double> mus(config.maxCategories);
	vector<unsigned int> catActives(config.maxCategories);

	for (unsigned int k = 0; k < config.maxCategories; k++) {
		string catIStr = catIStrBase + CatCont::catIndexString(k) + "]";

		mus[k] = param.at("catMu" + catIStr);
		catActives[k] = (unsigned int)param.at("catActive" + catIStr);
	}

	vector<unsigned int> activeCatActives = catActives;
	vector<unsigned int> inactiveCatActives = catActives;

	activeCatActives[catIndex] = 1;
	inactiveCatActives[catIndex] = 0;

	double activeDens = _scaledCatMuDensity(mus, activeCatActives, catIndex, config.catMuPriorApproximationPrecision);
	double inactiveDens = _scaledCatMuDensity(mus, inactiveCatActives, catIndex, config.catMuPriorApproximationPrecision);

	double numDens = (catActives[catIndex] == 1) ? activeDens : inactiveDens;
	double denDens = activeDens + inactiveDens;

	return std::log(numDens / denDens);
}

*/

} // namespace MemFirst
} // namespace CatCont
