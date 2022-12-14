#include "CCM_BayesianModel.h"


namespace CatCont {

	void Bayesian::_doMhOverrides(void) {

		if (this->runConfig.verbose) {
			stringstream ss;
			ss << "Number of MH overrides: " << config.overrides.mhTunings.size();
			logMessage("", ss.str());
		}

		for (auto it = config.overrides.mhTunings.begin(); it != config.overrides.mhTunings.end(); it++) {

			string name = it->first;

			if (mhTuningSd.find(name) != mhTuningSd.end()) {

				mhTuningSd[name] = it->second;
			} else {
				logMessage("_doMhOverrides", "Warning: Invalid MH tuning override key: \"" + name + "\" ignored.");
			}
		}
	}

	void Bayesian::_doPriorOverrides(void) {

		if (this->runConfig.verbose) {
			stringstream ss;
			ss << "Number of prior overrides: " << config.overrides.priors.size();
			logMessage("", ss.str());
		}

		for (auto it = config.overrides.priors.begin(); it != config.overrides.priors.end(); it++) {

			string name = it->first;

			if (this->priors.find(name) != this->priors.end()) {
				this->priors[name] = it->second;
			} else {
				logMessage("_doPriorOverrides", "Warning: Invalid prior override key: \"" + name + "\" ignored.");
			}

		}
	}

	void Bayesian::_doConstantParameterOverrides(void) {

		if (this->runConfig.verbose) {
			stringstream ss;
			ss << "Number of constant parameter value overrides: " << config.overrides.constantValues.size();
			logMessage("", ss.str());
		}

		vector<GibbsParameter*> parameters = gibbs.getParameters();
		for (unsigned int i = 0; i < parameters.size(); i++) {
			GibbsParameter* par = parameters[i];

			map<string, double>::iterator it = config.overrides.constantValues.find(par->name);

			if (it != config.overrides.constantValues.end()) {
				gibbs.replaceParameter(par->name, ConstantParameter(it->second, par->name, par->group));
			}

		}

		//Also check which values are provided but for which there is no parameter.
		for (map<string, double>::iterator it = config.overrides.constantValues.begin(); it != config.overrides.constantValues.end(); it++) {
			if (!gibbs.hasParameter(it->first)) {
				std::stringstream ss;
				ss << "Note: Constant value provided for \"" << it->first << "\", but there is no parameter by that name.";
				logMessage("_doConstantParameterOverrides", ss.str());
			}
		}

	}

	//This function must be called after createParameters().
	void Bayesian::_doStartingValueOverrides(void) {

		if (this->runConfig.verbose) {
			stringstream ss;
			ss << "Number of starting parameter value overrides: " << config.overrides.startingValues.size();
			logMessage("", ss.str());
		}

		vector<GibbsParameter*> parameters = gibbs.getParameters();
		for (unsigned int i = 0; i < parameters.size(); i++) {
			GibbsParameter* par = parameters[i];

			if (config.overrides.startingValues.find(par->name) != config.overrides.startingValues.end()) {
				//Found
				vector<double>& samples = par->getSamples();
				samples.clear();
				samples.push_back(config.overrides.startingValues.at(par->name));
			}
		}

		//Also check which values are provided but for which there is no parameter.
		for (map<string, double>::iterator it = config.overrides.startingValues.begin(); it != config.overrides.startingValues.end(); it++) {
			if (!gibbs.hasParameter(it->first)) {
				std::stringstream ss;
				ss << "Note: Starting value provided for \"" << it->first << "\", but there is no parameter by that name.";
				logMessage("_doStartingValueOverrides", ss.str());
			}
		}

	}


	// Note that this function does all calculations at the level of the participants. 
	// This is because the data is at the level of the participants.
	// In R, the participant level WAICs are aggregated to get a whole model WAIC.
	map<string, map<string, double>> calculateWAIC(const vector<ParticipantData>& allData,
		const ModelConfiguration& modelConfig,
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

				//Initialize to 0
				vector< double > currentLikelihoodSum(condData.study.size(), 0); 
				vector< double > currentLogLikelihoodSum = currentLikelihoodSum;

				//First index is the observation, second is the iteration.
				vector< vector< double > > currentLogLikelihoods;

				currentLogLikelihoods.resize(condData.study.size());
				for (unsigned int i = 0; i < currentLogLikelihoods.size(); i++) {
					currentLogLikelihoods[i].resize(posteriorIterations.size());
				}


				for (size_t iteration = 0; iteration < posteriorIterations.size(); iteration++) {

					const ParameterList& param = posteriorIterations[iteration];

					ParticipantParameters partParam = Bayesian::getParticipantParameters(param, pData.pnum, modelConfig.maxCategories);
					ConditionParameters condParam = Bayesian::getConditionParameters(param, condData.condition);

					CombinedParameters combinedParam = Bayesian::combineParameters(partParam, condParam, modelConfig.ranges, modelConfig.dataType);


					//For each iteration, get the likelihoods for all observation for this participant/condition pair
					vector<double> likelihoods;
					if (modelConfig.dataType == DataType::Circular) {
						likelihoods = Circular::betweenAndWithinLikelihood(combinedParam, condData, modelConfig.modelVariant);
					} else if (modelConfig.dataType == DataType::Linear) {
						likelihoods = Linear::betweenAndWithinLikelihood(combinedParam, condData, modelConfig);
					}

					//Store the likelihoods
					for (size_t obs = 0; obs < condData.study.size(); obs++) {

						currentLikelihoodSum[obs] += likelihoods[obs];

						double ll = log(likelihoods[obs]);

						currentLogLikelihoodSum[obs] += ll;

						currentLogLikelihoods[obs][iteration] = ll;
					}
				}

				//Process the likelihoods
				for (size_t obs = 0; obs < condData.study.size(); obs++) {

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

					for (size_t iteration = 0; iteration < posteriorIterations.size(); iteration++) {
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

		//These insertions do nothing special because the sort order is changed when put into the map. They could just be push_back.
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



	
	///////////////////////////////////////////
	// The decorrelating stuff should work in theory, but it doesn't seem to work well in practice, at least
	// not all of the time. I don't recommend using it.
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