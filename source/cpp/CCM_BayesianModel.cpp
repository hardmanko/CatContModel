#include "CCM_BayesianModel.h"


namespace CatCont {

	bool usingDecorrelatingSteps = false;
	const unsigned int catMuScaleFactorApproximationSteps = 60;

	Bayesian::Bayesian() {
		_catMuPriorData.catMuScaleFactorApproximationSteps = catMuScaleFactorApproximationSteps;
		config.calculateParticipantLikelihoods = false;
	}


	double Bayesian::_llFunction(const CombinedParameters& par, const ConditionData& data) const {
		double rval = 0;
		if (this->config.dataType == DataType::Circular) {
			rval = Circular::betweenAndWithinLL(par, data, config.modelVariant);
		} else if (this->config.dataType == DataType::Linear) {
			rval = Linear::betweenAndWithinLL(par, data, config.linearConfiguration);
		}
		
		return rval;
	}

	double Bayesian::singleParticipant_ll(const ParameterList& param, unsigned int pIndex) const {

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

	double Bayesian::_singleParticipantLLSamplingFunction(const ParameterList& param, unsigned int i) const {
		return singleParticipant_ll(param, i);
	}


	double Bayesian::singleCond_ll(const ParameterList& param, unsigned int condIndex) const {

		ConditionParameters conditionPar = _getConditionParameters(param, condIndex);

		double llSum = 0;

		for (unsigned int pIndex = 0; pIndex < data.participants.size(); pIndex++) {
			
			ParticipantParameters partPar = _getParticipantParameters(param, pIndex);

			CombinedParameters combinedParameters = _combineParameters(partPar, conditionPar);

			const ConditionData& cd = data.participants.at(pIndex).condData.at(condIndex);
			llSum += _llFunction(combinedParameters, cd);
		}

		return llSum;
	}

	double Bayesian::genericConditionParameter_ll(double thisParam, const ParameterList& param, unsigned int condIndex, string baseName) const {
		double ll = singleCond_ll(param, condIndex);
		double prior = cauchyLL(thisParam, priors.at(baseName + "_cond.loc"), priors.at(baseName + "_cond.scale"));
		return ll + prior;
	}

	double Bayesian::normalPriorParameter_ll(double thisParam, const ParameterList& param, unsigned int pIndex, string paramName) const {
		double ll = singleParticipant_ll(param, pIndex);
		double prior = normalLL(thisParam, param.at(paramName + ".mu"), param.at(paramName + ".var"));

		/*
		if (ll != ll) {
			Rcpp::Rcout << "Death for param: " << paramName << " & pnum " << pIndex << " with value " << thisParam << endl;
		}
		if (ll == std::numeric_limits<double>::infinity() || ll == -std::numeric_limits<double>::infinity()) {
			Rcpp::Rcout << "Infinite Death for param: " << paramName << " & pnum " << pIndex << " with value " << thisParam << endl;
		}
		*/

		return ll + prior;
	}


	//This is a combination of the g and h functions from the article
	double Bayesian::_catMuPenaltyDensity(const vector<double>& mus, const vector<unsigned int>& catActives, unsigned int k) const {

		double density = 1;

		for (unsigned int kp = 0; kp < mus.size(); kp++) {
			if (kp != k) {
				if (catActives[k] == 1 && catActives[kp] == 1) {

					double like;
					if (config.dataType == DataType::Circular) {
						like = vmLut.dVonMises(mus[k], mus[kp], this->_catMuPriorData.kappa);
					} else {
						//TODO: It really seems like just a normal would be fine, not truncated.
						like = Linear::dnorm(mus[k], mus[kp], _catMuPriorData.sd);
					}

					double ratio = 1 - like / _catMuPriorData.maxLikelihood;

					density *= (ratio * ratio); //square the ratio to make flatter bottoms on the notches
				}
			}
		}

		return density; //return non-log density
	}

	//This is an approximation of the S function from the article
	//mus are not reference because they are modified
	double Bayesian::_estimateCatMuScaleFactor(vector<double> mus, const vector<unsigned int>& catActives, unsigned int k, unsigned int steps) const {
		//circular
		double stepSize = 2 * PI / steps; 
		double startPoint = 0;

		if (config.dataType == DataType::Linear) {
			stepSize = (config.linearConfiguration.catMu.upper - config.linearConfiguration.catMu.lower) / steps; //TODO: Is this ok???
			startPoint = config.linearConfiguration.catMu.lower;
		}

		double densSum = 0;

		for (unsigned int i = 0; i < steps; i++) {
			mus[k] = startPoint + i * stepSize;

			densSum += _catMuPenaltyDensity(mus, catActives, k);
		}

		return 1.0 / (densSum * stepSize);
	}

	//This is the f function from the article
	double Bayesian::_scaledCatMuDensity(const vector<double>& mus, const vector<unsigned int>& catActives, unsigned int k, unsigned int steps) const {

		if (config.dataType == DataType::Linear) {
			if (mus[k] < config.linearConfiguration.catMu.lower || mus[k] > config.linearConfiguration.catMu.upper) {
				//Rcpp::Rcout << "mu[k] rejected: " << mus[k] << endl;
				return 0;
			}
		}

		//If this category is inactive, the density is just the height of a uniform distribution
		if (catActives[k] == 0) {
			if (config.dataType == DataType::Linear) {
				return 1.0 / (config.linearConfiguration.catMu.upper - config.linearConfiguration.catMu.lower); //TODO: Is this correct?
			} else {
				//circular
				return 1.0 / (2 * PI); 
			}
		}

		double unscaled = _catMuPenaltyDensity(mus, catActives, k);
		double scaleFactor = _estimateCatMuScaleFactor(mus, catActives, k, steps);

		return unscaled * scaleFactor;
	}

	double Bayesian::_catMuPenalityPrior(const ParameterList& param, unsigned int pIndex, unsigned int catIndex) const {
		string catIStrBase = "[" + data.participants.at(pIndex).pnum + ",";

		vector<double> mus(config.maxCategories);
		vector<unsigned int> catActives(config.maxCategories);

		for (unsigned int j = 0; j < config.maxCategories; j++) {
			string catIStr = catIStrBase + _convertToString(j) + "]";

			mus[j] = param.at("catMu" + catIStr);
			catActives[j] = (unsigned int)param.at("catActive" + catIStr);
		}

		double density = _scaledCatMuDensity(mus, catActives, catIndex, _catMuPriorData.catMuScaleFactorApproximationSteps);

		return std::log(density);
	}

	double Bayesian::_catActivePenaltyPrior(const ParameterList& param, unsigned int pIndex, unsigned int catIndex) const {
		string catIStrBase = "[" + data.participants.at(pIndex).pnum + ",";

		vector<double> mus(config.maxCategories);
		vector<unsigned int> catActives(config.maxCategories);

		for (unsigned int j = 0; j < config.maxCategories; j++) {
			string catIStr = catIStrBase + _convertToString(j) + "]";

			mus[j] = param.at("catMu" + catIStr);
			catActives[j] = (unsigned int)param.at("catActive" + catIStr);
		}

		vector<unsigned int> activeCatActives = catActives;
		vector<unsigned int> inactiveCatActives = catActives;

		activeCatActives[catIndex] = 1;
		inactiveCatActives[catIndex] = 0;

		double activeDens = _scaledCatMuDensity(mus, activeCatActives, catIndex, _catMuPriorData.catMuScaleFactorApproximationSteps);
		double inactiveDens = _scaledCatMuDensity(mus, inactiveCatActives, catIndex, _catMuPriorData.catMuScaleFactorApproximationSteps);

		double numDens = (catActives[catIndex] == 1) ? activeDens : inactiveDens;
		double denDens = activeDens + inactiveDens;

		return std::log(numDens / denDens);
	}

	double Bayesian::catActive_ll(double catActive, const ParameterList& param, unsigned int pIndex, unsigned int catIndex) const {

		double ll = singleParticipant_ll(param, pIndex);
		double prior = _catActivePenaltyPrior(param, pIndex, catIndex);

		return ll + prior;
	}

	double Bayesian::catMu_ll(double catMu, const ParameterList& param, unsigned int pIndex, unsigned int catIndex) const {

		double ll = 0;

		string catIStrBase = "[" + data.participants.at(pIndex).pnum + "," + _convertToString(catIndex) + "]";
		string catActiveName = "catActive" + catIStrBase;
		if (param.at(catActiveName) == 1.0) {
			ll = singleParticipant_ll(param, pIndex); //This only needs to be calculated for active categories.
		}

		double prior = _catMuPenalityPrior(param, pIndex, catIndex);

		return ll + prior;
	}




	double Bayesian::postMuSample(const ParameterList& param, string paramSetName, double mu0, double var0) const {
		vector<double> y = this->gibbs.getCurrentGroupValues(paramSetName);
		return normal_muPostSample(y, param.at(paramSetName + ".var"), mu0, var0);
	}

	double Bayesian::postVarSample(const ParameterList& param, string paramSetName, double a0, double b0) const {
		vector<double> y = this->gibbs.getCurrentGroupValues(paramSetName);
		return normal_varPostSample(y, param.at(paramSetName + ".mu"), a0, b0);
	}


	//This function produces deviates for the MH step.
	double Bayesian::_catActiveDeviateFunction(double active) {
		if (active == 1.0) {
			return 0.0;
		}
		return 1.0;
	}

	double Bayesian::_sdParameterTransformation(double sd, const SDRanges& ranges, DataType dataType) {

		if (dataType == DataType::Circular) {
			//transform from infinite space to bounded space
			sd = clamp(sd, ranges.minSd, ranges.maxSd);

			//convert to precision
			double kappa = Circular::sdDeg_to_precRad(sd);

			//clamp again. TODO: This is still not needed.
			double clampedKappa = clamp(kappa, ranges.minPrecision, ranges.maxPrecision);

			return clampedKappa;
		} else if (dataType == DataType::Linear) {
			sd = clamp(sd, ranges.minSd, ranges.maxSd); //This should be treated differently for normal. The minSd can be small.
			return sd;
		}


		return 0;
	}

	void Bayesian::createParameters(void) {

		using namespace std::placeholders;

		gibbs.clear();

		/////////////////////////////////////////////
		// participant parameters
		for (unsigned int i = 0; i < data.participants.size(); i++) {

			string istr = "[" + data.participants.at(i).pnum + "]";

			vector<string> probParam;
			probParam.push_back("pMem");
			probParam.push_back("pBetween");
			probParam.push_back("pContBetween");
			probParam.push_back("pContWithin");
			probParam.push_back("pCatGuess");

			for (string param : probParam) {
				MH_Parameter par;
				par.name = param + istr;
				par.group = param;

				par.deviateFunction = bind(normalDeviate, _1, mhTuningSd.at(param));
				par.llFunction = bind(&Bayesian::normalPriorParameter_ll, this, _1, _2, i, param);

				gibbs.addParameter(par, uniformDeviate(-2, 2));
			}

			vector<string> sdParam;
			sdParam.push_back("contSD");
			sdParam.push_back("catSelectivity");
			sdParam.push_back("catSD");

			for (string param : sdParam) {
				MH_Parameter par;
				par.name = param + istr;
				par.group = param;

				par.deviateFunction = bind(normalDeviate, _1, mhTuningSd.at(param));
				par.llFunction = bind(&Bayesian::normalPriorParameter_ll, this, _1, _2, i, param);

				gibbs.addParameter(par, uniformDeviate(config.ranges.minSd + 2, 40));
			}



			for (unsigned int j = 0; j < config.maxCategories; j++) {

				string catIstr = "[" + _convertToString(data.participants.at(i).pnum) + "," + _convertToString(j) + "]";

				{
					MH_Parameter catActive;
					catActive.name = "catActive" + catIstr;
					catActive.group = "catActive";

					catActive.deviateFunction = &Bayesian::_catActiveDeviateFunction;
					catActive.llFunction = bind(&Bayesian::catActive_ll, this, _1, _2, i, j);

					gibbs.addParameter(catActive, 1.0); //Start all categories active.
				}

				{
					MH_Parameter catMu;
					catMu.name = "catMu" + catIstr;
					catMu.group = "catMu";


					double catMuCandidateSd = mhTuningSd.at("catMu"); //linear
					if (config.dataType == DataType::Circular) {
						catMuCandidateSd = Circular::degreesToRadians(catMuCandidateSd);
					}

					catMu.deviateFunction = bind(normalDeviate, _1, catMuCandidateSd); //TODO: Is this a good candidate distribution??? VM?
					catMu.llFunction = bind(&Bayesian::catMu_ll, this, _1, _2, i, j);

					//The start value is in a grid within the response range
					//TODO: Maybe only for linear?
					double stepSize = (data.responseRange.upper - data.responseRange.lower) / (double)config.maxCategories;
					double startValue = (j + 0.5) * stepSize + data.responseRange.lower;

					gibbs.addParameter(catMu, startValue);
				}

			} //end category parameters

		} 
		//end participant parameters
		/////////////////////////////////////////////


		/////////////////////////////////////////////
		// hierarchical priors on participant parameters
		vector<string> paramWithHierarchicalPrior = config.getParamWithHierachicalPriors();

		for (const string& param : paramWithHierarchicalPrior) {

			ConjugateParameter popMu;

			popMu.name = param + ".mu";
			popMu.group = "populationParam";

			popMu.samplingFunction = std::bind(&Bayesian::postMuSample, this, _1, param, priors.at(param + ".mu.mu"), priors.at(param + ".mu.var"));

			gibbs.addParameter(popMu, 0);


			ConjugateParameter popVar;

			popVar.name = param + ".var";
			popVar.group = "populationParam";

			popVar.samplingFunction = std::bind(&Bayesian::postVarSample, this, _1, param, priors.at(param + ".var.a"), priors.at(param + ".var.b"));

			gibbs.addParameter(popVar, 30);

		}
		// end hierarchical priors
		/////////////////////////////////////////////


		/////////////////////////////////////////////
		//condition parameters
		vector<string> conditionEffectsToCreate = config.getParamWithAndWithoutConditionEffects();
		
		for (unsigned int condIndex = 0; condIndex < data.conditionNames.size(); condIndex++) {

			string cstr = "[" + data.conditionNames[condIndex] + "]";

			for (const string& s : conditionEffectsToCreate) {

				if (condIndex == config.cornerstoneConditionIndex) { 

					DependentParameter par;

					par.name = s + "_cond" + cstr;
					par.group = s + "_cond";

					par.sourceParameter = par.name; //Make itself its source: That will just result in it keeping the same value each time it updates.

					gibbs.addParameter(par, 0);

				} else {

					MH_Parameter par;
					par.name = s + "_cond" + cstr;
					par.group = s + "_cond";

					par.deviateFunction = bind(normalDeviate, _1, mhTuningSd.at(par.group));
					par.llFunction = std::bind(&Bayesian::genericConditionParameter_ll, this, _1, _2, condIndex, s);

					gibbs.addParameter(par, 0);

				}
			}

		}
		// end condition parameters
		/////////////////////////////////////////////


		/////////////////////////////////////////////
		// decorrelating steps
		if (usingDecorrelatingSteps) {
			vector<string> decorrelatedParameters = config.getParamWithAndWithoutConditionEffects();

			for (const string& s : decorrelatedParameters) {

				DecorrelatingStep dcs;

				dcs.name = s + "_deCor";
				dcs.group = "deCor";

				dcs.currentValuesFunction = bind(&Bayesian::decorrelatingCurrent, this, _1, s);
				dcs.candidateValuesFunction = bind(&Bayesian::decorrelatingCandidate, this, _1, s, mhTuningSd.at(dcs.name));

				dcs.llFunction = bind(&Bayesian::decorrelating_ll, this, _1, _2, s);

				gibbs.addParameter(dcs, 0);

			}
		}
		// end decorrelating steps
		/////////////////////////////////////////////


		///////////////////
		//Then go back and zero out conditions (and decorrelating steps) that don't have condition parameters.
		//This seems stupid, but it allows for easier extensibility.
		vector<string> condParamToZero = config.getParamWithoutConditionEffects();

		for (unsigned int i = 0; i < condParamToZero.size(); i++) {
			gibbs.setParameterGroupToConstantValue(condParamToZero[i] + "_cond", 0);

			if (usingDecorrelatingSteps) {
				gibbs.replaceParameter(condParamToZero[i] + "_deCor", ConstantParameter(0, condParamToZero[i] + "_deCor", "deCor"), 0);
			}
		}


		//between only variant: set pBetween and pContWithin to 1 (and also their condition effects to 0)
		if (config.modelVariant == ModelVariant::BetweenItem) {

			vector<string> pToSet;
			pToSet.push_back("pBetween");
			pToSet.push_back("pContWithin");

			for (string& s : pToSet) {
				gibbs.setParameterGroupToConstantValue(s, 100); //100 being essentially 1 once transformed, 
																   //but it actually doesn't matter because the between likelihood function is used.
				gibbs.setParameterGroupToConstantValue(s + "_cond", 0);

				//Get rid of population parameters
				gibbs.replaceParameter(s + ".mu", ConstantParameter(0, s + ".mu", "populationParam"), 0);
				gibbs.replaceParameter(s + ".var", ConstantParameter(1, s + ".var", "populationParam"), 1);

				if (usingDecorrelatingSteps) {
					gibbs.replaceParameter(s + "_deCor", ConstantParameter(0, s + "_deCor", "deCor"), 0);
				}
			}
		}

		//within only variant: set pBetween and pContWithin to 0 (and also their condition effects to 0)
		if (config.modelVariant == ModelVariant::WithinItem) {
			vector<string> pToSet;
			pToSet.push_back("pBetween");
			pToSet.push_back("pContBetween");

			for (string& s : pToSet) {
				gibbs.setParameterGroupToConstantValue(s, -100); //-100 being essentially 0 once transformed, 
																//but it actually doesn't matter because the within likelihood function is used.
				gibbs.setParameterGroupToConstantValue(s + "_cond", 0);

				//Get rid of population parameters
				gibbs.replaceParameter(s + ".mu", ConstantParameter(0, s + ".mu", "populationParam"), 0);
				gibbs.replaceParameter(s + ".var", ConstantParameter(1, s + ".var", "populationParam"), 1);

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

				gibbs.setParameterGroupToConstantValue(it->first, it->second);

				gibbs.setParameterGroupToConstantValue(it->first + "_cond", 0);

				//Get rid of population parameters
				gibbs.replaceParameter(it->first + ".mu", ConstantParameter(0, it->first + ".mu", "populationParam"), 0);
				gibbs.replaceParameter(it->first + ".var", ConstantParameter(1, it->first + ".var", "populationParam"), 1);

				if (usingDecorrelatingSteps) {
					gibbs.replaceParameter(it->first + "_deCor", ConstantParameter(0, it->first + "_deCor", "deCor"), 0);
				}

			}

			//Also set to 0 for all catActive
			gibbs.setParameterGroupToConstantValue("catActive", 0);
		}


		if (config.calculateParticipantLikelihoods) {
			for (unsigned int i = 0; i < data.participants.size(); i++) {

				string istr = "[" + data.participants.at(i).pnum + "]";

				DependentParameter spLL;
				spLL.name = "participantLL" + istr;
				spLL.group = "participantLL";

				spLL.samplingFunction = std::bind(&Bayesian::_singleParticipantLLSamplingFunction, this, _1, i);

				gibbs.addParameter(spLL);
			}
		}
		

		this->_doStartingValueOverrides();

		this->_doConstantParameterOverrides();
	}



	//Participant level parameters with estimated priors on them. 
	//Basically, every parameter except the catMu and catActive.
	vector<string> Bayesian::Configuration::getParamWithHierachicalPriors(void) const {
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

	vector<string> Bayesian::Configuration::getParamWithAndWithoutConditionEffects(void) const {
		return this->getParamWithHierachicalPriors();
	}

	vector<string> Bayesian::Configuration::getParamWithConditionEffects(void) const {
		return this->paramWithConditionEffects;
	}

	vector<string> Bayesian::Configuration::getParamWithoutConditionEffects(void) const {

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





	ParticipantParameters Bayesian::getParticipantParameters(const ParameterList& param, string pnum, unsigned int maxCategories) {
		string istr = "[" + pnum + "]";

		ParticipantParameters part;
		part.pMem = param.at("pMem" + istr);
		part.pBetween = param.at("pBetween" + istr);
		part.pContBetween = param.at("pContBetween" + istr);
		part.pContWithin = param.at("pContWithin" + istr);
		part.pCatGuess = param.at("pCatGuess" + istr);

		part.contSD = param.at("contSD" + istr);

		string catIStrBase = "[" + pnum + ",";

		for (unsigned int j = 0; j < maxCategories; j++) {
			string catIstr = catIStrBase + _convertToString(j) + "]";

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

	ConditionParameters Bayesian::getConditionParameters(const ParameterList& param, string condName) {

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

	CombinedParameters Bayesian::combineParameters(const ParticipantParameters& part, const ConditionParameters& cond, const SDRanges& ranges, DataType dataType) {
		CombinedParameters rval;

		rval.pMem = logitInverse(part.pMem + cond.pMem);
		rval.pBetween = logitInverse(part.pBetween + cond.pBetween);
		rval.pContBetween = logitInverse(part.pContBetween + cond.pContBetween);
		rval.pContWithin = logitInverse(part.pContWithin + cond.pContWithin);
		rval.pCatGuess = logitInverse(part.pCatGuess + cond.pCatGuess);

		rval.contSD = _sdParameterTransformation(part.contSD + cond.contSD, ranges, dataType);
		rval.cat.selectivity = _sdParameterTransformation(part.cat.selectivity + cond.cat.selectivity, ranges, dataType);
		rval.cat.SD = _sdParameterTransformation(part.cat.SD + cond.cat.SD, ranges, dataType);

		rval.cat.mu = part.cat.mu; //no condition parameters or transformations. the only possible transformation would be fmod(mu, 2*PI), but it is unnecessary

		return rval;
	}

	//Gets the untransformed parameters associated with a single participant. NO TRANSFORMATIONS ARE APPLIED HERE.
	ParticipantParameters Bayesian::_getParticipantParameters(const ParameterList& param, unsigned int pIndex) const {
		return getParticipantParameters(param, this->data.participants.at(pIndex).pnum, this->config.maxCategories);
	}

	//Gets the untransformed parameters associated with a single condition
	ConditionParameters Bayesian::_getConditionParameters(const ParameterList& param, unsigned int condIndex) const {
		return getConditionParameters(param, this->data.conditionNames[condIndex]);
	}

	//Combines the untransformed participant and condition parameters, tranforming them to the mainfest space
	CombinedParameters Bayesian::_combineParameters(const ParticipantParameters& part, const ConditionParameters& cond) const {
		return combineParameters(part, cond, this->config.ranges, this->config.dataType);
	}



	void Bayesian::setData(vector<ParticipantData> partData) {

		data = CatCont::Data(); //reset the data struct

		data.participants = partData;

		//Get the condition names and cornerstone condition
		config.cornerstoneConditionIndex = -1;

		for (unsigned int j = 0; j < data.participants.front().condData.size(); j++) {
			string conditionName = data.participants.front().condData[j].condition;
			data.conditionNames.push_back(conditionName);

			if (conditionName == this->config.cornerstoneConditionName) {
				config.cornerstoneConditionIndex = j;
			}
		}

		//Get the data ranges
		for (unsigned int i = 0; i < data.participants.size(); i++) {
			for (unsigned int j = 0; j < data.conditionNames.size(); j++) {

				const ConditionData& cd = data.participants[i].condData[j];

				for (unsigned int k = 0; k < cd.study.size(); k++) {

					data.studyRange.lower = min(data.studyRange.lower, cd.study[k]);
					data.studyRange.upper = max(data.studyRange.upper, cd.study[k]);

					data.responseRange.lower = min(data.responseRange.lower, cd.response[k]);
					data.responseRange.upper = max(data.responseRange.upper, cd.response[k]);
				} // k
			} // j
		} // i

	}

	//This function must be called after createParameters().
	void Bayesian::setParameterStartingValues(map<string, double> vals) {

		vector<GibbsParameter*> parameters = gibbs.getParameters();
		for (unsigned int i = 0; i < parameters.size(); i++) {
			GibbsParameter* par = parameters[i];

			if (vals.find(par->name) != vals.end()) {
				//Found
				vector<double>& samples = par->getSamples();
				samples.clear();
				samples.push_back( vals.at(par->name) );
			} else {
				std::stringstream ss;
				ss << "Note: No starting value found for \"" << par->name << "\".";
				logMessage("setParameterStartingValues", ss.str());
			}
		}

		//Also check which values are provided but for which there is no parameter.
		for (map<string, double>::iterator it = vals.begin(); it != vals.end(); it++) {
			if (!gibbs.hasParameter(it->first)) {
				std::stringstream ss;
				ss << "Note: Starting value provided for \"" << it->first << "\", but there is no parameter by that name.";
				logMessage("setParameterStartingValues", ss.str());
			}
		}

	}

	void Bayesian::setPriors(void) {

		priors.clear();

		vector<string> probParams;
		probParams.push_back("pMem");
		probParams.push_back("pBetween");
		probParams.push_back("pContBetween");
		probParams.push_back("pContWithin");
		probParams.push_back("pCatGuess");


		for (const string& s : probParams) {

			priors[s + ".mu.mu"] = 0;
			priors[s + ".mu.var"] = pow(3, 2); //sd = 3
			priors[s + ".var.a"] = 0.1;
			priors[s + ".var.b"] = 0.1;

			priors[s + "_cond.loc"] = 0;
			priors[s + "_cond.scale"] = 2;

		}

		vector<string> sdParams;
		sdParams.push_back("contSD");
		sdParams.push_back("catSelectivity");
		sdParams.push_back("catSD");

		for (const string& s : sdParams) {

			priors[s + ".mu.mu"] = 20;
			priors[s + ".mu.var"] = pow(20, 2); //sd = 20
			priors[s + ".var.a"] = 0.1;
			priors[s + ".var.b"] = 0.1;

			priors[s + "_cond.loc"] = 0;
			priors[s + "_cond.scale"] = 5;

		}

		priors["catMuPriorSD"] = 12;


		//Once all of the defaults are in, do the overrides
		_doPriorOverrides();


		//After the overrides are in, calculate these things.

		_catMuPriorData.sd = priors["catMuPriorSD"];
		if (config.dataType == DataType::Circular) {

			_catMuPriorData.kappa = Circular::sdDeg_to_precRad(_catMuPriorData.sd);
			_catMuPriorData.maxLikelihood = vmLut.dVonMises(0, 0, _catMuPriorData.kappa);

		} else if (config.dataType == DataType::Linear) {

			//TODO
			_catMuPriorData.kappa = std::numeric_limits<double>::infinity();
			_catMuPriorData.maxLikelihood = Linear::dnorm(0, 0, _catMuPriorData.sd); //TODO: Not truncated???

		}
	}



	void Bayesian::setMhTuning(void) {
		mhTuningSd.clear();

		//participant parameters
		mhTuningSd["pMem"] = 0.3;
		mhTuningSd["pBetween"] = 0.9;
		mhTuningSd["pContBetween"] = 0.4;
		mhTuningSd["pContWithin"] = 0.7;
		mhTuningSd["pCatGuess"] = 0.7;

		mhTuningSd["contSD"] = 2;
		mhTuningSd["catSelectivity"] = 2;
		mhTuningSd["catSD"] = 1;

		mhTuningSd["catMu"] = 5; //degrees


		//condition effects
		mhTuningSd["pMem_cond"] = 0.2;
		mhTuningSd["pBetween_cond"] = 0.4;
		mhTuningSd["pContBetween_cond"] = 0.2;
		mhTuningSd["pContWithin_cond"] = 0.2;
		mhTuningSd["pCatGuess_cond"] = 0.3;

		mhTuningSd["contSD_cond"] = 1;
		mhTuningSd["catSelectivity_cond"] = 0.8;
		mhTuningSd["catSD_cond"] = 0.5;


		
		if (usingDecorrelatingSteps) {

			for (auto& sd : mhTuningSd) {
				//To account for the decorrelating steps, the step sizes have to be smaller for the rest of the parameters.
				sd.second /= 2;
			}
			//Then unadjust things without condition effects
			mhTuningSd["catMu"] *= 2;


			//decorrelating parameters
			mhTuningSd["pMem_deCor"] = 0.1;
			mhTuningSd["pBetween_deCor"] = 0.2;
			mhTuningSd["pContBetween_deCor"] = 0.2;
			mhTuningSd["pContWithin_deCor"] = 0.3;
			mhTuningSd["pCatGuess_deCor"] = 0.1;

			mhTuningSd["contSD_deCor"] = 2;

			mhTuningSd["catSelectivity_deCor"] = 1;
			mhTuningSd["catSD_deCor"] = 1;
		}


		_doMhOverrides();

	}


	map<string, map<string, double>> calculateWAIC(const vector<ParticipantData>& allData,
		const Bayesian::Configuration& modelConfig,
		const vector< ParameterList >& posteriorIterations)
	{

		double LPPD = 0;
		double P_1 = 0;
		double P_2 = 0;

		vector<string> pnums(allData.size());
		vector<double> LPPD_part(allData.size(), 0);
		vector<double> P_1_part = LPPD_part;
		vector<double> P_2_part = LPPD_part;
		vector<double> WAIC_1_part = LPPD_part;
		vector<double> WAIC_2_part = LPPD_part;


		for (unsigned int pIndex = 0; pIndex < allData.size(); pIndex++) {

			const ParticipantData& pData = allData[pIndex];
			pnums[pIndex] = pData.pnum;

			for (unsigned int condIndex = 0; condIndex < pData.condData.size(); condIndex++) {

				const ConditionData& condData = pData.condData[condIndex];
				
				
				vector< double > currentLikelihoodSum(condData.study.size(), 0); //Initialize to 0
				vector< double > currentLogLikelihoodSum = currentLikelihoodSum; //Also initialize to 0

				//First index is the observation, second is the iteration.
				vector< vector< double > > currentLogLikelihoods;

				currentLogLikelihoods.resize(condData.study.size());
				for (unsigned int i = 0; i < currentLogLikelihoods.size(); i++) {
					currentLogLikelihoods[i].resize(posteriorIterations.size());
				}


				for (unsigned int iteration = 0; iteration < posteriorIterations.size(); iteration++) {

					const ParameterList& param = posteriorIterations[iteration];

					ParticipantParameters partParam = Bayesian::getParticipantParameters(param, pData.pnum, modelConfig.maxCategories);
					ConditionParameters condParam = Bayesian::getConditionParameters(param, condData.condition);

					CombinedParameters combinedParam = Bayesian::combineParameters(partParam, condParam, modelConfig.ranges, modelConfig.dataType);


					//For each iteration, get the likelihoods for all observation for this participant/condition pair
					vector<double> likelihoods;
					if (modelConfig.dataType == DataType::Circular) {
						likelihoods = Circular::betweenAndWithinLikelihood(combinedParam, condData, modelConfig.modelVariant);
					} else if (modelConfig.dataType == DataType::Linear) {
						likelihoods = Linear::betweenAndWithinLikelihood(combinedParam, condData, modelConfig.linearConfiguration);
					}

					//Store the likelihoods
					for (unsigned int obs = 0; obs < condData.study.size(); obs++) {

						currentLikelihoodSum[obs] += likelihoods[obs];

						double ll = log(likelihoods[obs]);

						currentLogLikelihoodSum[obs] += ll;

						currentLogLikelihoods[obs][iteration] = ll;
					}

				}

				//Process the likelihoods
				for (unsigned int obs = 0; obs < condData.study.size(); obs++) {


					//Convert from sum to mean
					double logOfMeanL = log(currentLikelihoodSum[obs] / posteriorIterations.size());
					double meanOfLogL = currentLogLikelihoodSum[obs] / posteriorIterations.size();


					//Add on to LPPD
					LPPD += logOfMeanL;
					LPPD_part[pIndex] += logOfMeanL;

					//Add on to P_1
					double a = logOfMeanL;
					double b = meanOfLogL;

					P_1 += 2 * (a - b);
					P_1_part[pIndex] += 2 * (a - b);

					//Add on to P_2: sample variance of LL
					double varOfLLAccum = 0;

					for (unsigned int iteration = 0; iteration < posteriorIterations.size(); iteration++) {
						double dOfLL = currentLogLikelihoods[obs][iteration] - meanOfLogL;
						varOfLLAccum += dOfLL * dOfLL;
					}

					double varOfLL = varOfLLAccum / (posteriorIterations.size() - 1);
					P_2 += varOfLL;
					P_2_part[pIndex] += varOfLL;


				} //obs


			} //condIndex

			WAIC_1_part[pIndex] = -2 * (LPPD_part[pIndex] - P_1_part[pIndex]);
			WAIC_2_part[pIndex] = -2 * (LPPD_part[pIndex] - P_2_part[pIndex]);

		} //pIndex

		double WAIC_1 = -2 * (LPPD - P_1);
		double WAIC_2 = -2 * (LPPD - P_2);

		//These insertions do nothing because the sort order is changed when put into the map. They could just be push_back.
		pnums.insert(pnums.begin(), "Total");
		WAIC_1_part.insert(WAIC_1_part.begin(), WAIC_1);
		WAIC_2_part.insert(WAIC_2_part.begin(), WAIC_2);

		LPPD_part.insert(LPPD_part.begin(), LPPD);
		P_1_part.insert(P_1_part.begin(), P_1);
		P_2_part.insert(P_2_part.begin(), P_2);

		map<string, map<string, double> > rval;

		for (unsigned int i = 0; i < pnums.size(); i++) {
			map<string, double> temp;
			temp["WAIC_1"] = WAIC_1_part[i];
			temp["WAIC_2"] = WAIC_2_part[i];
			temp["P_1"] = P_1_part[i];
			temp["P_2"] = P_2_part[i];
			temp["LPPD"] = LPPD_part[i];
			rval[pnums[i]] = temp;
		}

		return rval;
	}


	//The decorrelating stuff should work in theory, but it doesn't seem to work well in practice, at least
	//not all of the time. I don't recommend using it.
	ParameterList Bayesian::_getDecorrelatingValues(const ParameterList& param, string paramSetName, double deviate) const {
		ParameterList rval;

		for (unsigned int i = 0; i < data.participants.size(); i++) {
			string istr = paramSetName + "[" + data.participants.at(i).pnum + "]";
			rval[istr] = param.at(istr) + deviate;
		}

		for (unsigned int j = 0; j < data.conditionNames.size(); j++) {
			if (j != config.cornerstoneConditionIndex) {
				string istr = paramSetName + "_cond[" + data.conditionNames[j] + "]";
				rval[istr] = param.at(istr) - deviate;
			}
		}

		return rval;
	}

	ParameterList Bayesian::decorrelatingCurrent(const ParameterList& param, string paramSetName) const {
		return _getDecorrelatingValues(param, paramSetName, 0);
	}

	ParameterList Bayesian::decorrelatingCandidate(const ParameterList& param, string paramSetName, double mhSd) const {
		double deviate = normalDeviate(0, mhSd);

		return _getDecorrelatingValues(param, paramSetName, deviate);
	}

	double Bayesian::decorrelating_ll(const ParameterList& theseParam, const ParameterList& allParam, string paramSetName) const {

		vector<ConditionParameters> condPar(data.conditionNames.size());

		double prior = 0;

		for (unsigned int condIndex = 0; condIndex < data.conditionNames.size(); condIndex++) {

			condPar[condIndex] = _getConditionParameters(allParam, condIndex);

			if (condIndex != config.cornerstoneConditionIndex) {
				prior += cauchyLL(theseParam.at(paramSetName + "_cond[" + data.conditionNames[condIndex] + "]"),
					priors.at(paramSetName + "_cond.loc"), priors.at(paramSetName + "_cond.scale"));
			}
		}

		double llSum = 0;
		for (unsigned int pIndex = 0; pIndex < data.participants.size(); pIndex++) {


			ParticipantParameters partPar = _getParticipantParameters(allParam, pIndex);

			for (unsigned int condIndex = 0; condIndex < data.conditionNames.size(); condIndex++) {
				CombinedParameters combinedParameters = _combineParameters(partPar, condPar[condIndex]);

				const ConditionData& cd = data.participants.at(pIndex).condData.at(condIndex);
				llSum += _llFunction(combinedParameters, cd);
			}

			prior += normalLL(theseParam.at(paramSetName + "[" + data.participants.at(pIndex).pnum + "]"),
				allParam.at(paramSetName + ".mu"), allParam.at(paramSetName + ".var"));
		}

		return llSum + prior;
	}


} // namespace CatCont

