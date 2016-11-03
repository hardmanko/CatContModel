#include "CCM_BayesianModel.h"


namespace CatCont {

	void Bayesian::_doMhOverrides(void) {

		stringstream ss;
		ss << "Number of MH overrides: " << overrides.mhTunings.size();
		logMessage("", ss.str());

		for (auto it = overrides.mhTunings.begin(); it != overrides.mhTunings.end(); it++) {

			string name = it->first;

			if (mhTuningSd.find(name) != mhTuningSd.end()) {

				mhTuningSd[name] = it->second;
			} else {
				logMessage("_doMhOverrides", "Warning: Invalid MH tuning override key: \"" + name + "\" ignored.");
			}
		}
	}

	void Bayesian::_doPriorOverrides(void) {

		stringstream ss;
		ss << "Number of prior overrides: " << overrides.priors.size();
		logMessage("", ss.str());

		for (auto it = overrides.priors.begin(); it != overrides.priors.end(); it++) {

			string name = it->first;

			if (priors.find(name) != priors.end()) {
				priors[name] = it->second;
			} else {
				logMessage("_doPriorOverrides", "Warning: Invalid prior override key: \"" + name + "\" ignored.");
			}

		}
	}

	void Bayesian::_doConstantParameterOverrides(void) {

		stringstream ss;
		ss << "Number of constant parameter value overrides: " << overrides.constantValues.size();
		logMessage("", ss.str());


		vector<GibbsParameter*> parameters = gibbs.getParameters();
		for (unsigned int i = 0; i < parameters.size(); i++) {
			GibbsParameter* par = parameters[i];

			map<string, double>::iterator it = overrides.constantValues.find(par->name);

			if (it != overrides.constantValues.end()) {
				gibbs.replaceParameter(par->name, ConstantParameter(it->second, par->name, par->group));
			}

		}

		//Also check which values are provided but for which there is no parameter.
		for (map<string, double>::iterator it = overrides.constantValues.begin(); it != overrides.constantValues.end(); it++) {
			if (!gibbs.hasParameter(it->first)) {
				std::stringstream ss;
				ss << "Note: Constant value provided for \"" << it->first << "\", but there is no parameter by that name.";
				logMessage("setParameterStartingValues", ss.str());
			}
		}

	}

	//This function must be called after createParameters().
	void Bayesian::_doStartingValueOverrides(void) {

		stringstream ss;
		ss << "Number of starting parameter value overrides: " << overrides.startingValues.size();
		logMessage("", ss.str());

		vector<GibbsParameter*> parameters = gibbs.getParameters();
		for (unsigned int i = 0; i < parameters.size(); i++) {
			GibbsParameter* par = parameters[i];

			if (overrides.startingValues.find(par->name) != overrides.startingValues.end()) {
				//Found
				vector<double>& samples = par->getSamples();
				samples.clear();
				samples.push_back(overrides.startingValues.at(par->name));
			} 
			/*
			else {
				std::stringstream ss;
				ss << "Note: No starting value found for \"" << par->name << "\".";
				logMessage("setParameterStartingValues", ss.str());
			}
			*/
		}

		//Also check which values are provided but for which there is no parameter.
		for (map<string, double>::iterator it = overrides.startingValues.begin(); it != overrides.startingValues.end(); it++) {
			if (!gibbs.hasParameter(it->first)) {
				std::stringstream ss;
				ss << "Note: Starting value provided for \"" << it->first << "\", but there is no parameter by that name.";
				logMessage("setParameterStartingValues", ss.str());
			}
		}

	}



	//TODO: This is unfinished and maybe wrongheaded
	void Bayesian::_doParameterEqualityConstraints(void) {

		//assume that you have this
		map<string, string> mapping;


		using namespace std::placeholders;

		vector<GibbsParameter*> param = this->gibbs.getParameters();
		vector<string> paramNames(param.size());
		for (unsigned int i = 0; i < param.size(); i++) {
			paramNames[i] = param[i]->name;
		}

		EqualityConstraints eq;
		eq.setup(mapping, paramNames, this->data.conditionNames);

		for (GibbsParameter* target : param) {
			string source = eq.getSourceParameter(target->name);
			if (source == "FREE_PARAMETER") {


				//If a condition effect, get the condition group indices
				if (target->name.find("_cond[") != string::npos) {

					string baseName = extractBaseParameterName(target->name);
					string condName = extractIndex(target->name);

					if (condName == config.cornerstoneConditionName) {
						//Don't change anything for the cornerstone condition
						continue;
					}

					//Get the index of this condition. Consider putting this into a function...
					unsigned int condIndex = -1;
					for (unsigned int i = 0; i < data.conditionNames.size(); i++) {
						if (data.conditionNames[i] == condName) {
							condIndex = i;
							break;
						}
					}

					vector<unsigned int> conditionIndices = eq.getEqualConditionIndices(baseName, condIndex);

					MH_Parameter par;
					par.name = target->name;
					par.group = target->group;

					par.deviateFunction = bind(normalDeviate, _1, mhTuningSd.at(par.group));
					par.llFunction = std::bind(&Bayesian::multiConditionParameter_ll, this, _1, _2, conditionIndices, baseName);

					gibbs.replaceParameter(target->name, par, 0);

				}

			} else {
				DependentParameter dp;
				dp.name = target->name;
				dp.group = target->group;
				dp.sourceParameter = source; //set this parameter to use the source
				gibbs.replaceParameter(target->name, dp);
			}
		}


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




	/////////////////////////////////
	// EqualityConstraints
	/////////////////////////////////

	string Bayesian::EqualityConstraints::FreeParameter = "FREE_PARAMETER";

	bool Bayesian::EqualityConstraints::setup(const IndividualMappings& mappings, const vector<string>& conditionNames, const vector<string>& allParameterNames) {
		individualMappings = mappings;

		if (allParameterNames.size() > 0) {
			vector<string> keysToErase;

			for (auto it = individualMappings.begin(); it != individualMappings.end(); it++) {

				bool firstNotFound = find(allParameterNames.begin(), allParameterNames.end(), it->first) == allParameterNames.end();
				bool secondNotFound = find(allParameterNames.begin(), allParameterNames.end(), it->second) == allParameterNames.end();

				if (it->first == FreeParameter) {
					firstNotFound = false;
				}
				if (it->second == FreeParameter) {
					secondNotFound = false;
				}

				if (firstNotFound) {
					logMessage("EqualityConstraints", "Parameter " + it->first + " does not exist and will be ignored.");
				}

				if (secondNotFound) {
					logMessage("EqualityConstraints", "Parameter " + it->second + " does not exist and will be ignored.");
				}

				if (firstNotFound || secondNotFound) {
					keysToErase.push_back(it->first);
				}
			}

			for (string s : keysToErase) {
				individualMappings.erase(s);
			}
		}

		bool simplifySuccess = simplifyIndividualMappings(&individualMappings);
		if (!simplifySuccess) {
			return false;
		}

		groupMappings = calculateEqualConditionIndices(individualMappings, conditionNames);
		return true;
	}



	//This should include the given condition
	const vector<unsigned int>& Bayesian::EqualityConstraints::getEqualConditionIndices(string param, unsigned int cond) const {
		return groupMappings.at(param).at(cond);
	}

	//Basically, read out of individualMappings
	string Bayesian::EqualityConstraints::getSourceParameter(string parameter) const {
		if (individualMappings.find(parameter) == individualMappings.end()) {
			return "PARAMETER_NOT_FOUND";
		}

		return individualMappings.at(parameter);
	}

	//Account for parameters not yet mentioned(?)
	Bayesian::EqualityConstraints::IndividualMappings Bayesian::EqualityConstraints::incorporateAdditionalParameters(IndividualMappings mappings, const vector<string>& allNames) const {
		for (string s : allNames) {
			if (mappings.find(s) == mappings.end()) {
				mappings[s] = FreeParameter;
			}
		}
		return mappings;
	}

	//Reduce to simplest form
	bool Bayesian::EqualityConstraints::simplifyIndividualMappings(IndividualMappings* mapping) const {

		//get parameter names
		vector<string> names;
		for (auto it = mapping->begin(); it != mapping->end(); it++) {
			names.push_back(it->first);
		}

		bool mappingChanged = false;
		for (const string& currentTarget : names) {

			string currentSource = mapping->at(currentTarget);

			if (currentSource == FreeParameter) {
				continue;
			}

			set<string> cycleHistory;
			cycleHistory.insert(currentSource);

			string mostSimpleSource = currentSource;
			while (true) {
				string nextPossibleSource = mapping->at(mostSimpleSource);

				if (cycleHistory.find(nextPossibleSource) != cycleHistory.end()) {
					logMessage("EqualityConstraints::simplifyIndividualMappings()", "Error: Equality constraint cycle detected.");
					return false;
				} else {
					cycleHistory.insert(nextPossibleSource);
				}

				if (nextPossibleSource != FreeParameter) {
					mostSimpleSource = nextPossibleSource;
				} else {
					break;
				}
			}
			if (mostSimpleSource != currentSource) {
				mapping->at(currentTarget) = mostSimpleSource;
				mappingChanged = true;
			}

		}

		if (mappingChanged) {
			return simplifyIndividualMappings(mapping);
		}

		return true;
	}



	//conditionNames comes from the Data struct.
	Bayesian::EqualityConstraints::GroupMappings Bayesian::EqualityConstraints::calculateEqualConditionIndices(const IndividualMappings& mapping, const vector<string>& conditionNames) const {

		GroupMappings rval;

		map<string, unsigned int> cnameToIndex;
		for (unsigned int i = 0; i < conditionNames.size(); i++) {
			cnameToIndex[conditionNames[i]] = i;
		}

		vector<string> names;
		for (auto it = mapping.begin(); it != mapping.end(); it++) {
			if (it->first.find("_cond[") != string::npos) {
				names.push_back(it->first);
			}
		}

		for (const string& possibleSource : names) {

			string baseParameterName = extractBaseParameterName(possibleSource);

			string sourceCondName = extractIndex(possibleSource);

			vector<string> cNames;
			vector<unsigned int> cInds;

			//If this is a free parameter, add it to its own group
			if (mapping.at(possibleSource) == FreeParameter) {
				string sourceCondName = extractIndex(possibleSource);

				cNames.push_back(sourceCondName);
				cInds.push_back(cnameToIndex.at(sourceCondName));
			}

			for (const string& possibleTarget : names) {

				if (mapping.at(possibleTarget) == possibleSource) {
					string targetCondName = extractIndex(possibleTarget);

					cNames.push_back(targetCondName);
					cInds.push_back(cnameToIndex.at(targetCondName));
				}

			}

			unsigned int sourceConditionIndex = cnameToIndex.at(sourceCondName);
			rval[baseParameterName][sourceConditionIndex] = cInds;
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
}